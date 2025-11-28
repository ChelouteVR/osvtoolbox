/**
 * OSV Toolbox - DJI OSV File Extractor/Recomposer
 *
 * Cross-platform tool for extracting and recomposing DJI OSV files (MP4 variant)
 *
 * Compilation:
 *   g++ -std=c++17 -O2 -o osvtoolbox osvtoolbox.cpp
 *
 * On Windows (MSVC):
 *   cl /EHsc /O2 osvtoolbox.cpp /Fe:osvtoolbox.exe
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <map>
#include <memory>

// Cross-platform byte order handling
#ifdef _MSC_VER
    #include <stdlib.h>
    #define bswap_32(x) _byteswap_ulong(x)
    #define bswap_64(x) _byteswap_uint64(x)
#elif defined(__APPLE__)
    #include <libkern/OSByteOrder.h>
    #define bswap_32(x) OSSwapInt32(x)
    #define bswap_64(x) OSSwapInt64(x)
#else
    #include <byteswap.h>
#endif

// Read big-endian values
inline uint32_t read_be32(const uint8_t* p) {
    return (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) |
           (uint32_t(p[2]) << 8) | uint32_t(p[3]);
}

inline uint64_t read_be64(const uint8_t* p) {
    return (uint64_t(read_be32(p)) << 32) | uint64_t(read_be32(p + 4));
}

inline void write_be32(uint8_t* p, uint32_t v) {
    p[0] = (v >> 24) & 0xFF;
    p[1] = (v >> 16) & 0xFF;
    p[2] = (v >> 8) & 0xFF;
    p[3] = v & 0xFF;
}

inline void write_be64(uint8_t* p, uint64_t v) {
    write_be32(p, v >> 32);
    write_be32(p + 4, v & 0xFFFFFFFF);
}

// MP4 Box structure
struct Box {
    uint64_t offset;      // Position in file
    uint64_t size;        // Total box size
    uint64_t header_size; // Header size (8 or 16)
    std::string type;     // 4-char box type
    std::vector<uint8_t> data;  // Box content (for small boxes)
    std::vector<std::shared_ptr<Box>> children;  // Child boxes for containers

    uint64_t content_offset() const { return offset + header_size; }
    uint64_t content_size() const { return size - header_size; }
};

// Sample entry for tracks
struct Sample {
    uint64_t offset;
    uint32_t size;
};

// Track information
struct Track {
    uint32_t track_id;
    std::string handler_type;  // vide, soun, meta
    std::string codec;         // hvc1, mp4a, djmd, dbgi
    std::vector<Sample> samples;
    std::vector<uint8_t> stsd_data;  // Sample description
    uint32_t timescale;
    uint64_t duration;
    uint32_t width;
    uint32_t height;
    uint32_t sample_rate;
    uint32_t channels;
};

class OSVParser {
public:
    std::string filename;
    std::ifstream file;
    uint64_t file_size;

    std::vector<std::shared_ptr<Box>> top_boxes;
    std::shared_ptr<Box> moov_box;
    std::shared_ptr<Box> mdat_box;
    std::shared_ptr<Box> udta_box;
    std::shared_ptr<Box> meta_box;
    std::shared_ptr<Box> camd_box;
    std::vector<uint8_t> ftyp_data;
    std::vector<uint8_t> free_index_data;  // The free box with index

    std::vector<Track> tracks;

    OSVParser(const std::string& path) : filename(path) {}

    bool open() {
        file.open(filename, std::ios::binary);
        if (!file) {
            std::cerr << "Error: Cannot open file " << filename << std::endl;
            return false;
        }

        file.seekg(0, std::ios::end);
        file_size = file.tellg();
        file.seekg(0, std::ios::beg);

        return true;
    }

    void close() {
        if (file.is_open()) {
            file.close();
        }
    }

    std::shared_ptr<Box> read_box(uint64_t max_end = 0) {
        uint64_t pos = file.tellg();
        if (max_end == 0) max_end = file_size;
        if (pos >= max_end) return nullptr;

        uint8_t header[8];
        file.read(reinterpret_cast<char*>(header), 8);
        if (file.gcount() < 8) return nullptr;

        auto box = std::make_shared<Box>();
        box->offset = pos;
        box->size = read_be32(header);
        box->type = std::string(reinterpret_cast<char*>(header + 4), 4);
        box->header_size = 8;

        if (box->size == 1) {
            uint8_t ext[8];
            file.read(reinterpret_cast<char*>(ext), 8);
            box->size = read_be64(ext);
            box->header_size = 16;
        } else if (box->size == 0) {
            box->size = max_end - pos;
        }

        return box;
    }

    void read_box_data(std::shared_ptr<Box> box) {
        file.seekg(box->content_offset());
        box->data.resize(box->content_size());
        file.read(reinterpret_cast<char*>(box->data.data()), box->content_size());
    }

    bool is_container(const std::string& type) {
        static const char* containers[] = {
            "moov", "trak", "mdia", "minf", "stbl", "udta", "dinf", "edts", nullptr
        };
        for (const char** p = containers; *p; ++p) {
            if (type == *p) return true;
        }
        return false;
    }

    void parse_container(std::shared_ptr<Box> box) {
        uint64_t end_pos = box->offset + box->size;
        file.seekg(box->content_offset());

        while ((uint64_t)file.tellg() < end_pos) {
            auto child = read_box(end_pos);
            if (!child) break;

            if (is_container(child->type)) {
                parse_container(child);
            } else if (child->type == "meta") {
                // meta box has version/flags before children
                file.seekg(child->content_offset() + 4);
                uint64_t meta_end = child->offset + child->size;
                while ((uint64_t)file.tellg() < meta_end) {
                    auto meta_child = read_box(meta_end);
                    if (!meta_child) break;
                    child->children.push_back(meta_child);
                    file.seekg(meta_child->offset + meta_child->size);
                }
            }

            box->children.push_back(child);
            file.seekg(child->offset + child->size);
        }
    }

    std::shared_ptr<Box> find_box(std::shared_ptr<Box> parent, const std::string& type) {
        for (auto& child : parent->children) {
            if (child->type == type) return child;
        }
        return nullptr;
    }

    std::shared_ptr<Box> find_box_path(std::shared_ptr<Box> parent, const std::vector<std::string>& path) {
        auto current = parent;
        for (const auto& type : path) {
            current = find_box(current, type);
            if (!current) return nullptr;
        }
        return current;
    }

    bool parse_top_level() {
        file.seekg(0);

        while ((uint64_t)file.tellg() < file_size) {
            auto box = read_box();
            if (!box) break;

            std::cout << "Found box: " << box->type << " at " << box->offset
                      << " size " << box->size << std::endl;

            if (box->type == "ftyp") {
                read_box_data(box);
                ftyp_data.resize(box->size);
                file.seekg(box->offset);
                file.read(reinterpret_cast<char*>(ftyp_data.data()), box->size);
            } else if (box->type == "free") {
                // Check if this is the index free box (has covr/snal/camd)
                if (box->size > 100) {
                    read_box_data(box);
                    if (box->data.size() >= 4 &&
                        memcmp(box->data.data(), "covr", 4) == 0) {
                        free_index_data.resize(box->size);
                        file.seekg(box->offset);
                        file.read(reinterpret_cast<char*>(free_index_data.data()), box->size);
                    }
                }
            } else if (box->type == "mdat") {
                mdat_box = box;
            } else if (box->type == "moov") {
                moov_box = box;
                parse_container(box);
            } else if (box->type == "camd") {
                camd_box = box;
            }

            top_boxes.push_back(box);
            file.seekg(box->offset + box->size);
        }

        if (!moov_box) {
            std::cerr << "Error: No moov box found" << std::endl;
            return false;
        }

        // Find udta and meta inside moov
        udta_box = find_box(moov_box, "udta");
        meta_box = find_box(moov_box, "meta");

        return true;
    }

    void parse_stco(std::shared_ptr<Box> stco, std::vector<uint64_t>& offsets) {
        file.seekg(stco->content_offset());
        uint8_t buf[8];
        file.read(reinterpret_cast<char*>(buf), 8);
        uint32_t entry_count = read_be32(buf + 4);

        offsets.resize(entry_count);
        for (uint32_t i = 0; i < entry_count; i++) {
            file.read(reinterpret_cast<char*>(buf), 4);
            offsets[i] = read_be32(buf);
        }
    }

    void parse_co64(std::shared_ptr<Box> co64, std::vector<uint64_t>& offsets) {
        file.seekg(co64->content_offset());
        uint8_t buf[12];
        file.read(reinterpret_cast<char*>(buf), 8);
        uint32_t entry_count = read_be32(buf + 4);

        offsets.resize(entry_count);
        for (uint32_t i = 0; i < entry_count; i++) {
            file.read(reinterpret_cast<char*>(buf), 8);
            offsets[i] = read_be64(buf);
        }
    }

    void parse_stsz(std::shared_ptr<Box> stsz, std::vector<uint32_t>& sizes) {
        file.seekg(stsz->content_offset());
        uint8_t buf[12];
        file.read(reinterpret_cast<char*>(buf), 12);
        uint32_t sample_size = read_be32(buf + 4);
        uint32_t sample_count = read_be32(buf + 8);

        sizes.resize(sample_count);
        if (sample_size != 0) {
            for (uint32_t i = 0; i < sample_count; i++) {
                sizes[i] = sample_size;
            }
        } else {
            for (uint32_t i = 0; i < sample_count; i++) {
                file.read(reinterpret_cast<char*>(buf), 4);
                sizes[i] = read_be32(buf);
            }
        }
    }

    void parse_stsc(std::shared_ptr<Box> stsc,
                    std::vector<uint32_t>& first_chunk,
                    std::vector<uint32_t>& samples_per_chunk) {
        file.seekg(stsc->content_offset());
        uint8_t buf[12];
        file.read(reinterpret_cast<char*>(buf), 8);
        uint32_t entry_count = read_be32(buf + 4);

        first_chunk.resize(entry_count);
        samples_per_chunk.resize(entry_count);

        for (uint32_t i = 0; i < entry_count; i++) {
            file.read(reinterpret_cast<char*>(buf), 12);
            first_chunk[i] = read_be32(buf);
            samples_per_chunk[i] = read_be32(buf + 4);
            // sample_description_index at buf + 8, ignored
        }
    }

    bool parse_tracks() {
        for (auto& child : moov_box->children) {
            if (child->type != "trak") continue;

            Track track;

            // Get track ID from tkhd
            auto tkhd = find_box(child, "tkhd");
            if (tkhd) {
                file.seekg(tkhd->content_offset());
                uint8_t buf[100];
                file.read(reinterpret_cast<char*>(buf), 92);
                uint8_t version = buf[0];
                if (version == 0) {
                    track.track_id = read_be32(buf + 12);
                    track.width = read_be32(buf + 76) >> 16;
                    track.height = read_be32(buf + 80) >> 16;
                } else {
                    track.track_id = read_be32(buf + 20);
                    track.width = read_be32(buf + 84) >> 16;
                    track.height = read_be32(buf + 88) >> 16;
                }
            }

            // Get handler type from mdia/hdlr
            auto mdia = find_box(child, "mdia");
            if (!mdia) continue;

            // Get timescale from mdhd
            auto mdhd = find_box(mdia, "mdhd");
            if (mdhd) {
                file.seekg(mdhd->content_offset());
                uint8_t buf[32];
                file.read(reinterpret_cast<char*>(buf), 32);
                uint8_t version = buf[0];
                if (version == 0) {
                    track.timescale = read_be32(buf + 12);
                    track.duration = read_be32(buf + 16);
                } else {
                    track.timescale = read_be32(buf + 20);
                    track.duration = read_be64(buf + 24);
                }
            }

            auto hdlr = find_box(mdia, "hdlr");
            if (hdlr) {
                file.seekg(hdlr->content_offset() + 8);
                char handler[5] = {0};
                file.read(handler, 4);
                track.handler_type = handler;
            }

            // Get sample info from minf/stbl
            auto minf = find_box(mdia, "minf");
            if (!minf) continue;

            auto stbl = find_box(minf, "stbl");
            if (!stbl) continue;

            // Get codec from stsd
            auto stsd = find_box(stbl, "stsd");
            if (stsd) {
                file.seekg(stsd->content_offset() + 8);
                uint8_t buf[4];
                file.read(reinterpret_cast<char*>(buf), 4);
                uint32_t entry_size = read_be32(buf);
                file.read(reinterpret_cast<char*>(buf), 4);
                track.codec = std::string(reinterpret_cast<char*>(buf), 4);

                // Save entire stsd content for rebuilding
                track.stsd_data.resize(stsd->content_size());
                file.seekg(stsd->content_offset());
                file.read(reinterpret_cast<char*>(track.stsd_data.data()),
                         stsd->content_size());
            }

            // Parse chunk offsets
            std::vector<uint64_t> chunk_offsets;
            auto stco = find_box(stbl, "stco");
            auto co64 = find_box(stbl, "co64");
            if (stco) {
                parse_stco(stco, chunk_offsets);
            } else if (co64) {
                parse_co64(co64, chunk_offsets);
            }

            // Parse sample sizes
            std::vector<uint32_t> sample_sizes;
            auto stsz = find_box(stbl, "stsz");
            if (stsz) {
                parse_stsz(stsz, sample_sizes);
            }

            // Parse sample-to-chunk mapping
            std::vector<uint32_t> first_chunk, samples_per_chunk;
            auto stsc = find_box(stbl, "stsc");
            if (stsc) {
                parse_stsc(stsc, first_chunk, samples_per_chunk);
            }

            // Build sample list
            if (!chunk_offsets.empty() && !sample_sizes.empty()) {
                uint32_t sample_idx = 0;
                uint32_t stsc_idx = 0;

                for (uint32_t chunk_idx = 0; chunk_idx < chunk_offsets.size(); chunk_idx++) {
                    // Find the stsc entry for this chunk
                    while (stsc_idx + 1 < first_chunk.size() &&
                           first_chunk[stsc_idx + 1] <= chunk_idx + 1) {
                        stsc_idx++;
                    }

                    uint32_t samples_in_chunk = samples_per_chunk[stsc_idx];
                    uint64_t offset = chunk_offsets[chunk_idx];

                    for (uint32_t s = 0; s < samples_in_chunk && sample_idx < sample_sizes.size(); s++) {
                        Sample sample;
                        sample.offset = offset;
                        sample.size = sample_sizes[sample_idx];
                        track.samples.push_back(sample);

                        offset += sample.size;
                        sample_idx++;
                    }
                }
            }

            std::cout << "Track " << track.track_id << ": " << track.codec
                      << " (" << track.handler_type << ") - "
                      << track.samples.size() << " samples" << std::endl;

            tracks.push_back(track);
        }

        return !tracks.empty();
    }

    bool extract_track_raw(int track_idx, const std::string& output_path) {
        if (track_idx < 0 || track_idx >= (int)tracks.size()) {
            std::cerr << "Invalid track index: " << track_idx << std::endl;
            return false;
        }

        const Track& track = tracks[track_idx];
        std::ofstream out(output_path, std::ios::binary);
        if (!out) {
            std::cerr << "Cannot create output file: " << output_path << std::endl;
            return false;
        }

        // Also save sample sizes to a .sizes file
        std::string sizes_path = output_path + ".sizes";
        std::ofstream sizes_out(sizes_path, std::ios::binary);

        std::vector<uint8_t> buffer;
        for (const auto& sample : track.samples) {
            buffer.resize(sample.size);
            file.seekg(sample.offset);
            file.read(reinterpret_cast<char*>(buffer.data()), sample.size);
            out.write(reinterpret_cast<char*>(buffer.data()), sample.size);

            // Write sample size as big-endian uint32
            uint8_t size_buf[4];
            write_be32(size_buf, sample.size);
            sizes_out.write(reinterpret_cast<char*>(size_buf), 4);
        }

        out.close();
        sizes_out.close();
        std::cout << "Extracted track " << track_idx << " to " << output_path
                  << " (" << track.samples.size() << " samples)" << std::endl;
        return true;
    }

    // Create a minimal MP4 file with just one track
    bool extract_video_track_mp4(int track_idx, const std::string& output_path) {
        if (track_idx < 0 || track_idx >= (int)tracks.size()) {
            std::cerr << "Invalid track index: " << track_idx << std::endl;
            return false;
        }

        const Track& track = tracks[track_idx];
        if (track.handler_type != "vide") {
            std::cerr << "Track " << track_idx << " is not a video track" << std::endl;
            return false;
        }

        std::ofstream out(output_path, std::ios::binary);
        if (!out) {
            std::cerr << "Cannot create output file: " << output_path << std::endl;
            return false;
        }

        // Calculate mdat size
        uint64_t mdat_content_size = 0;
        for (const auto& sample : track.samples) {
            mdat_content_size += sample.size;
        }

        // Write ftyp
        out.write(reinterpret_cast<const char*>(ftyp_data.data()), ftyp_data.size());

        // Write mdat header
        uint64_t mdat_offset = ftyp_data.size();
        uint64_t mdat_size = 8 + mdat_content_size;
        bool use_64bit_mdat = mdat_size > 0xFFFFFFFF;

        if (use_64bit_mdat) {
            uint8_t mdat_header[16];
            write_be32(mdat_header, 1);
            memcpy(mdat_header + 4, "mdat", 4);
            write_be64(mdat_header + 8, 16 + mdat_content_size);
            out.write(reinterpret_cast<char*>(mdat_header), 16);
            mdat_offset += 16;
        } else {
            uint8_t mdat_header[8];
            write_be32(mdat_header, 8 + mdat_content_size);
            memcpy(mdat_header + 4, "mdat", 4);
            out.write(reinterpret_cast<char*>(mdat_header), 8);
            mdat_offset = ftyp_data.size() + 8;
        }

        // Write mdat content and build chunk offsets
        std::vector<uint64_t> new_offsets;
        std::vector<uint8_t> buffer;
        uint64_t current_offset = mdat_offset;

        for (const auto& sample : track.samples) {
            new_offsets.push_back(current_offset);
            buffer.resize(sample.size);
            file.seekg(sample.offset);
            file.read(reinterpret_cast<char*>(buffer.data()), sample.size);
            out.write(reinterpret_cast<char*>(buffer.data()), sample.size);
            current_offset += sample.size;
        }

        uint64_t moov_offset = out.tellp();

        // Build moov box
        std::vector<uint8_t> moov_data;

        // Helper to append data
        auto append = [&moov_data](const void* data, size_t size) {
            const uint8_t* p = static_cast<const uint8_t*>(data);
            moov_data.insert(moov_data.end(), p, p + size);
        };

        auto append_box_header = [&moov_data](const char* type, uint32_t size) {
            uint8_t header[8];
            write_be32(header, size);
            memcpy(header + 4, type, 4);
            moov_data.insert(moov_data.end(), header, header + 8);
        };

        // We'll build the moov content then write the header
        size_t moov_content_start = 8;
        moov_data.resize(8); // Reserve space for moov header

        // mvhd (movie header)
        {
            uint8_t mvhd[108];
            memset(mvhd, 0, sizeof(mvhd));
            write_be32(mvhd, 108);  // size
            memcpy(mvhd + 4, "mvhd", 4);
            mvhd[8] = 0;  // version
            // timescale at offset 20
            write_be32(mvhd + 20, track.timescale);
            // duration at offset 24
            write_be32(mvhd + 24, track.duration);
            // preferred rate 1.0 at offset 28
            write_be32(mvhd + 28, 0x00010000);
            // preferred volume 1.0 at offset 32
            mvhd[32] = 0x01; mvhd[33] = 0x00;
            // matrix at offset 44
            write_be32(mvhd + 44, 0x00010000);
            write_be32(mvhd + 60, 0x00010000);
            write_be32(mvhd + 84, 0x40000000);
            // next track id at offset 104
            write_be32(mvhd + 104, 2);
            append(mvhd, 108);
        }

        // trak box
        size_t trak_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);  // Reserve for trak header

        // tkhd (track header)
        {
            uint8_t tkhd[92];
            memset(tkhd, 0, sizeof(tkhd));
            write_be32(tkhd, 92);
            memcpy(tkhd + 4, "tkhd", 4);
            tkhd[8] = 0;  // version
            write_be32(tkhd + 8, 0x00000003);  // flags: enabled, in movie
            // track id at offset 12
            write_be32(tkhd + 20, 1);
            // duration at offset 28
            write_be32(tkhd + 28, track.duration);
            // matrix at offset 48
            write_be32(tkhd + 48, 0x00010000);
            write_be32(tkhd + 64, 0x00010000);
            write_be32(tkhd + 88, 0x40000000);
            // width at offset 84
            write_be32(tkhd + 84, track.width << 16);
            // height at offset 88
            write_be32(tkhd + 88, track.height << 16);
            append(tkhd, 92);
        }

        // mdia box
        size_t mdia_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // mdhd
        {
            uint8_t mdhd[32];
            memset(mdhd, 0, sizeof(mdhd));
            write_be32(mdhd, 32);
            memcpy(mdhd + 4, "mdhd", 4);
            mdhd[8] = 0;  // version
            write_be32(mdhd + 12, track.timescale);
            write_be32(mdhd + 16, track.duration);
            mdhd[20] = 0x55; mdhd[21] = 0xC4;  // language "und"
            append(mdhd, 32);
        }

        // hdlr
        {
            const char* hdlr_name = "VideoHandler";
            size_t name_len = strlen(hdlr_name) + 1;
            size_t hdlr_size = 32 + name_len;

            uint8_t hdlr[32];
            memset(hdlr, 0, sizeof(hdlr));
            write_be32(hdlr, hdlr_size);
            memcpy(hdlr + 4, "hdlr", 4);
            hdlr[8] = 0;  // version
            memcpy(hdlr + 16, "vide", 4);
            append(hdlr, 32);
            moov_data.insert(moov_data.end(), hdlr_name, hdlr_name + name_len);
        }

        // minf box
        size_t minf_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // vmhd
        {
            uint8_t vmhd[20];
            memset(vmhd, 0, sizeof(vmhd));
            write_be32(vmhd, 20);
            memcpy(vmhd + 4, "vmhd", 4);
            write_be32(vmhd + 8, 0x00000001);  // version + flags
            append(vmhd, 20);
        }

        // dinf/dref
        {
            uint8_t dinf[36];
            memset(dinf, 0, sizeof(dinf));
            write_be32(dinf, 36);
            memcpy(dinf + 4, "dinf", 4);
            write_be32(dinf + 8, 28);
            memcpy(dinf + 12, "dref", 4);
            write_be32(dinf + 20, 1);  // entry count
            write_be32(dinf + 24, 12);  // url size
            memcpy(dinf + 28, "url ", 4);
            write_be32(dinf + 32, 0x00000001);  // flags: self-contained
            append(dinf, 36);
        }

        // stbl box
        size_t stbl_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // stsd - use original
        {
            size_t stsd_size = 8 + track.stsd_data.size();
            uint8_t stsd_header[8];
            write_be32(stsd_header, stsd_size);
            memcpy(stsd_header + 4, "stsd", 4);
            append(stsd_header, 8);
            append(track.stsd_data.data(), track.stsd_data.size());
        }

        // stts (time to sample)
        {
            uint8_t stts[24];
            memset(stts, 0, sizeof(stts));
            write_be32(stts, 24);
            memcpy(stts + 4, "stts", 4);
            write_be32(stts + 12, 1);  // entry count
            write_be32(stts + 16, track.samples.size());  // sample count
            write_be32(stts + 20, 1001);  // sample delta (for ~30fps at 30000 timescale)
            append(stts, 24);
        }

        // stss (sync samples / keyframes) - find keyframes from original
        std::vector<uint32_t> keyframes;
        for (auto& child : moov_box->children) {
            if (child->type != "trak") continue;
            auto tkhd = find_box(child, "tkhd");
            if (!tkhd) continue;
            file.seekg(tkhd->content_offset() + 12);
            uint8_t buf[4];
            file.read(reinterpret_cast<char*>(buf), 4);
            uint32_t tid = read_be32(buf);
            if (tid != track.track_id) continue;

            auto stbl = find_box_path(child, {"mdia", "minf", "stbl"});
            if (!stbl) continue;
            auto stss = find_box(stbl, "stss");
            if (stss) {
                file.seekg(stss->content_offset() + 4);
                file.read(reinterpret_cast<char*>(buf), 4);
                uint32_t count = read_be32(buf);
                for (uint32_t i = 0; i < count; i++) {
                    file.read(reinterpret_cast<char*>(buf), 4);
                    keyframes.push_back(read_be32(buf));
                }
            }
            break;
        }

        if (!keyframes.empty()) {
            size_t stss_size = 16 + keyframes.size() * 4;
            uint8_t stss_header[16];
            write_be32(stss_header, stss_size);
            memcpy(stss_header + 4, "stss", 4);
            stss_header[8] = stss_header[9] = stss_header[10] = stss_header[11] = 0;
            write_be32(stss_header + 12, keyframes.size());
            append(stss_header, 16);
            for (uint32_t kf : keyframes) {
                uint8_t buf[4];
                write_be32(buf, kf);
                append(buf, 4);
            }
        }

        // stsc (sample to chunk) - one sample per chunk for simplicity
        {
            uint8_t stsc[28];
            memset(stsc, 0, sizeof(stsc));
            write_be32(stsc, 28);
            memcpy(stsc + 4, "stsc", 4);
            write_be32(stsc + 12, 1);  // entry count
            write_be32(stsc + 16, 1);  // first chunk
            write_be32(stsc + 20, 1);  // samples per chunk
            write_be32(stsc + 24, 1);  // sample description index
            append(stsc, 28);
        }

        // stsz (sample sizes)
        {
            size_t stsz_size = 20 + track.samples.size() * 4;
            uint8_t stsz_header[20];
            write_be32(stsz_header, stsz_size);
            memcpy(stsz_header + 4, "stsz", 4);
            memset(stsz_header + 8, 0, 8);  // version/flags and sample_size=0
            write_be32(stsz_header + 16, track.samples.size());
            append(stsz_header, 20);
            for (const auto& sample : track.samples) {
                uint8_t buf[4];
                write_be32(buf, sample.size);
                append(buf, 4);
            }
        }

        // stco/co64 (chunk offsets)
        bool use_co64 = (current_offset > 0xFFFFFFFF);
        if (use_co64) {
            size_t co64_size = 16 + new_offsets.size() * 8;
            uint8_t co64_header[16];
            write_be32(co64_header, co64_size);
            memcpy(co64_header + 4, "co64", 4);
            memset(co64_header + 8, 0, 4);
            write_be32(co64_header + 12, new_offsets.size());
            append(co64_header, 16);
            for (uint64_t off : new_offsets) {
                uint8_t buf[8];
                write_be64(buf, off);
                append(buf, 8);
            }
        } else {
            size_t stco_size = 16 + new_offsets.size() * 4;
            uint8_t stco_header[16];
            write_be32(stco_header, stco_size);
            memcpy(stco_header + 4, "stco", 4);
            memset(stco_header + 8, 0, 4);
            write_be32(stco_header + 12, new_offsets.size());
            append(stco_header, 16);
            for (uint64_t off : new_offsets) {
                uint8_t buf[4];
                write_be32(buf, off);
                append(buf, 4);
            }
        }

        // Fix stbl size
        write_be32(&moov_data[stbl_start], moov_data.size() - stbl_start);
        memcpy(&moov_data[stbl_start + 4], "stbl", 4);

        // Fix minf size
        write_be32(&moov_data[minf_start], moov_data.size() - minf_start);
        memcpy(&moov_data[minf_start + 4], "minf", 4);

        // Fix mdia size
        write_be32(&moov_data[mdia_start], moov_data.size() - mdia_start);
        memcpy(&moov_data[mdia_start + 4], "mdia", 4);

        // Fix trak size
        write_be32(&moov_data[trak_start], moov_data.size() - trak_start);
        memcpy(&moov_data[trak_start + 4], "trak", 4);

        // Fix moov size
        write_be32(&moov_data[0], moov_data.size());
        memcpy(&moov_data[4], "moov", 4);

        out.write(reinterpret_cast<char*>(moov_data.data()), moov_data.size());
        out.close();

        std::cout << "Extracted video track " << track_idx << " to " << output_path << std::endl;
        return true;
    }

    bool extract_audio_track_mp4(int track_idx, const std::string& output_path) {
        if (track_idx < 0 || track_idx >= (int)tracks.size()) {
            std::cerr << "Invalid track index: " << track_idx << std::endl;
            return false;
        }

        const Track& track = tracks[track_idx];
        if (track.handler_type != "soun") {
            std::cerr << "Track " << track_idx << " is not an audio track" << std::endl;
            return false;
        }

        // Similar to video but with audio-specific boxes
        std::ofstream out(output_path, std::ios::binary);
        if (!out) {
            std::cerr << "Cannot create output file: " << output_path << std::endl;
            return false;
        }

        // Calculate mdat size
        uint64_t mdat_content_size = 0;
        for (const auto& sample : track.samples) {
            mdat_content_size += sample.size;
        }

        // Write ftyp
        out.write(reinterpret_cast<const char*>(ftyp_data.data()), ftyp_data.size());

        // Write mdat
        uint8_t mdat_header[8];
        write_be32(mdat_header, 8 + mdat_content_size);
        memcpy(mdat_header + 4, "mdat", 4);
        out.write(reinterpret_cast<char*>(mdat_header), 8);

        uint64_t mdat_offset = ftyp_data.size() + 8;
        std::vector<uint64_t> new_offsets;
        std::vector<uint8_t> buffer;
        uint64_t current_offset = mdat_offset;

        for (const auto& sample : track.samples) {
            new_offsets.push_back(current_offset);
            buffer.resize(sample.size);
            file.seekg(sample.offset);
            file.read(reinterpret_cast<char*>(buffer.data()), sample.size);
            out.write(reinterpret_cast<char*>(buffer.data()), sample.size);
            current_offset += sample.size;
        }

        // Build moov for audio
        std::vector<uint8_t> moov_data;

        auto append = [&moov_data](const void* data, size_t size) {
            const uint8_t* p = static_cast<const uint8_t*>(data);
            moov_data.insert(moov_data.end(), p, p + size);
        };

        moov_data.resize(8);

        // mvhd
        {
            uint8_t mvhd[108];
            memset(mvhd, 0, sizeof(mvhd));
            write_be32(mvhd, 108);
            memcpy(mvhd + 4, "mvhd", 4);
            write_be32(mvhd + 20, track.timescale);
            write_be32(mvhd + 24, track.duration);
            write_be32(mvhd + 28, 0x00010000);
            mvhd[32] = 0x01;
            write_be32(mvhd + 44, 0x00010000);
            write_be32(mvhd + 60, 0x00010000);
            write_be32(mvhd + 84, 0x40000000);
            write_be32(mvhd + 104, 2);
            append(mvhd, 108);
        }

        size_t trak_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // tkhd
        {
            uint8_t tkhd[92];
            memset(tkhd, 0, sizeof(tkhd));
            write_be32(tkhd, 92);
            memcpy(tkhd + 4, "tkhd", 4);
            write_be32(tkhd + 8, 0x00000003);
            write_be32(tkhd + 20, 1);
            write_be32(tkhd + 28, track.duration);
            write_be32(tkhd + 36, 0x01000000);  // volume
            write_be32(tkhd + 48, 0x00010000);
            write_be32(tkhd + 64, 0x00010000);
            write_be32(tkhd + 88, 0x40000000);
            append(tkhd, 92);
        }

        size_t mdia_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // mdhd
        {
            uint8_t mdhd[32];
            memset(mdhd, 0, sizeof(mdhd));
            write_be32(mdhd, 32);
            memcpy(mdhd + 4, "mdhd", 4);
            write_be32(mdhd + 12, track.timescale);
            write_be32(mdhd + 16, track.duration);
            mdhd[20] = 0x55; mdhd[21] = 0xC4;
            append(mdhd, 32);
        }

        // hdlr
        {
            const char* name = "SoundHandler";
            size_t name_len = strlen(name) + 1;
            size_t hdlr_size = 32 + name_len;
            uint8_t hdlr[32];
            memset(hdlr, 0, sizeof(hdlr));
            write_be32(hdlr, hdlr_size);
            memcpy(hdlr + 4, "hdlr", 4);
            memcpy(hdlr + 16, "soun", 4);
            append(hdlr, 32);
            moov_data.insert(moov_data.end(), name, name + name_len);
        }

        size_t minf_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // smhd
        {
            uint8_t smhd[16];
            memset(smhd, 0, sizeof(smhd));
            write_be32(smhd, 16);
            memcpy(smhd + 4, "smhd", 4);
            append(smhd, 16);
        }

        // dinf/dref
        {
            uint8_t dinf[36];
            memset(dinf, 0, sizeof(dinf));
            write_be32(dinf, 36);
            memcpy(dinf + 4, "dinf", 4);
            write_be32(dinf + 8, 28);
            memcpy(dinf + 12, "dref", 4);
            write_be32(dinf + 20, 1);
            write_be32(dinf + 24, 12);
            memcpy(dinf + 28, "url ", 4);
            write_be32(dinf + 32, 0x00000001);
            append(dinf, 36);
        }

        size_t stbl_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // stsd
        {
            size_t stsd_size = 8 + track.stsd_data.size();
            uint8_t stsd_header[8];
            write_be32(stsd_header, stsd_size);
            memcpy(stsd_header + 4, "stsd", 4);
            append(stsd_header, 8);
            append(track.stsd_data.data(), track.stsd_data.size());
        }

        // stts
        {
            uint8_t stts[24];
            memset(stts, 0, sizeof(stts));
            write_be32(stts, 24);
            memcpy(stts + 4, "stts", 4);
            write_be32(stts + 12, 1);
            write_be32(stts + 16, track.samples.size());
            write_be32(stts + 20, 1024);  // typical AAC frame duration
            append(stts, 24);
        }

        // stsc
        {
            uint8_t stsc[28];
            memset(stsc, 0, sizeof(stsc));
            write_be32(stsc, 28);
            memcpy(stsc + 4, "stsc", 4);
            write_be32(stsc + 12, 1);
            write_be32(stsc + 16, 1);
            write_be32(stsc + 20, 1);
            write_be32(stsc + 24, 1);
            append(stsc, 28);
        }

        // stsz
        {
            size_t stsz_size = 20 + track.samples.size() * 4;
            uint8_t stsz_header[20];
            write_be32(stsz_header, stsz_size);
            memcpy(stsz_header + 4, "stsz", 4);
            memset(stsz_header + 8, 0, 8);
            write_be32(stsz_header + 16, track.samples.size());
            append(stsz_header, 20);
            for (const auto& sample : track.samples) {
                uint8_t buf[4];
                write_be32(buf, sample.size);
                append(buf, 4);
            }
        }

        // stco
        {
            size_t stco_size = 16 + new_offsets.size() * 4;
            uint8_t stco_header[16];
            write_be32(stco_header, stco_size);
            memcpy(stco_header + 4, "stco", 4);
            memset(stco_header + 8, 0, 4);
            write_be32(stco_header + 12, new_offsets.size());
            append(stco_header, 16);
            for (uint64_t off : new_offsets) {
                uint8_t buf[4];
                write_be32(buf, off);
                append(buf, 4);
            }
        }

        write_be32(&moov_data[stbl_start], moov_data.size() - stbl_start);
        memcpy(&moov_data[stbl_start + 4], "stbl", 4);

        write_be32(&moov_data[minf_start], moov_data.size() - minf_start);
        memcpy(&moov_data[minf_start + 4], "minf", 4);

        write_be32(&moov_data[mdia_start], moov_data.size() - mdia_start);
        memcpy(&moov_data[mdia_start + 4], "mdia", 4);

        write_be32(&moov_data[trak_start], moov_data.size() - trak_start);
        memcpy(&moov_data[trak_start + 4], "trak", 4);

        write_be32(&moov_data[0], moov_data.size());
        memcpy(&moov_data[4], "moov", 4);

        out.write(reinterpret_cast<char*>(moov_data.data()), moov_data.size());
        out.close();

        std::cout << "Extracted audio track " << track_idx << " to " << output_path << std::endl;
        return true;
    }

    bool save_box_raw(std::shared_ptr<Box> box, const std::string& output_path) {
        if (!box) {
            std::cerr << "Box is null" << std::endl;
            return false;
        }

        std::ofstream out(output_path, std::ios::binary);
        if (!out) {
            std::cerr << "Cannot create output file: " << output_path << std::endl;
            return false;
        }

        std::vector<uint8_t> buffer(box->size);
        file.seekg(box->offset);
        file.read(reinterpret_cast<char*>(buffer.data()), box->size);
        out.write(reinterpret_cast<char*>(buffer.data()), box->size);
        out.close();

        std::cout << "Saved box " << box->type << " to " << output_path << std::endl;
        return true;
    }

    bool save_additional_boxes(const std::string& output_path) {
        std::ofstream out(output_path, std::ios::binary);
        if (!out) {
            std::cerr << "Cannot create output file: " << output_path << std::endl;
            return false;
        }

        // Write a header with box info
        // Format: count (4 bytes), then for each box: type (4), offset (8), size (8)
        // Then the actual box data

        std::vector<std::pair<std::string, std::shared_ptr<Box>>> boxes_to_save;

        if (udta_box) boxes_to_save.push_back({"udta", udta_box});
        if (meta_box) boxes_to_save.push_back({"meta", meta_box});
        if (camd_box) boxes_to_save.push_back({"camd", camd_box});

        // Also save the free index box
        if (!free_index_data.empty()) {
            // Write magic header
            const char* magic = "OSVT";
            out.write(magic, 4);

            // Write version
            uint8_t version[4] = {0, 0, 0, 1};
            out.write(reinterpret_cast<char*>(version), 4);

            // Write ftyp
            uint32_t ftyp_size = ftyp_data.size();
            uint8_t size_buf[4];
            write_be32(size_buf, ftyp_size);
            out.write(reinterpret_cast<char*>(size_buf), 4);
            out.write(reinterpret_cast<const char*>(ftyp_data.data()), ftyp_data.size());

            // Write free index
            uint32_t free_size = free_index_data.size();
            write_be32(size_buf, free_size);
            out.write(reinterpret_cast<char*>(size_buf), 4);
            out.write(reinterpret_cast<const char*>(free_index_data.data()), free_index_data.size());
        }

        // Write box count
        uint8_t count_buf[4];
        write_be32(count_buf, boxes_to_save.size());
        out.write(reinterpret_cast<char*>(count_buf), 4);

        // Write each box
        std::vector<uint8_t> buffer;
        for (const auto& [name, box] : boxes_to_save) {
            // Write box size
            uint8_t size_buf[8];
            write_be64(size_buf, box->size);
            out.write(reinterpret_cast<char*>(size_buf), 8);

            // Write box data
            buffer.resize(box->size);
            file.seekg(box->offset);
            file.read(reinterpret_cast<char*>(buffer.data()), box->size);
            out.write(reinterpret_cast<char*>(buffer.data()), box->size);
        }

        out.close();
        std::cout << "Saved additional boxes to " << output_path << std::endl;
        return true;
    }
};

std::string get_base_name(const std::string& path) {
    size_t pos = path.rfind('.');
    if (pos != std::string::npos) {
        return path.substr(0, pos);
    }
    return path;
}

int extract_data(const std::string& input_file) {
    OSVParser parser(input_file);

    if (!parser.open()) {
        return 1;
    }

    std::cout << "Parsing OSV file: " << input_file << std::endl;
    std::cout << "File size: " << parser.file_size << " bytes" << std::endl;

    if (!parser.parse_top_level()) {
        return 1;
    }

    if (!parser.parse_tracks()) {
        return 1;
    }

    std::string base = get_base_name(input_file);

    // Use ffmpeg to extract video and audio tracks (cleaner output)
    std::cout << "\nExtracting video track 1 with ffmpeg..." << std::endl;
    std::string cmd1 = "ffmpeg -y -v warning -i \"" + input_file + "\" -map 0:0 -c copy \"" + base + "-track1.mp4\"";
    system(cmd1.c_str());

    std::cout << "Extracting video track 2 with ffmpeg..." << std::endl;
    std::string cmd2 = "ffmpeg -y -v warning -i \"" + input_file + "\" -map 0:1 -c copy \"" + base + "-track2.mp4\"";
    system(cmd2.c_str());

    std::cout << "Extracting audio track with ffmpeg..." << std::endl;
    std::string cmd3 = "ffmpeg -y -v warning -i \"" + input_file + "\" -map 0:2 -c copy \"" + base + "-track3.mp4\"";
    system(cmd3.c_str());

    // Extract djmd/dbgi tracks as raw binary
    int djmd_idx = 1;
    int dbgi_idx = 1;

    for (size_t i = 0; i < parser.tracks.size(); i++) {
        const Track& track = parser.tracks[i];

        if (track.codec == "djmd") {
            std::string output = base + "-djmd" + std::to_string(djmd_idx) + ".bin";
            parser.extract_track_raw(i, output);
            djmd_idx++;
        } else if (track.codec == "dbgi") {
            std::string output = base + "-dbgi" + std::to_string(dbgi_idx) + ".bin";
            parser.extract_track_raw(i, output);
            dbgi_idx++;
        }
    }

    // Save additional boxes (udta, meta, camd)
    std::string boxes_output = base + "-additional-boxes.bin";
    parser.save_additional_boxes(boxes_output);

    parser.close();

    std::cout << "\nExtraction complete!" << std::endl;
    return 0;
}

// Helper class to parse input MP4 files for recompose
class MP4Reader {
public:
    std::ifstream file;
    uint64_t file_size;
    std::vector<Sample> samples;
    std::vector<uint8_t> stsd_data;
    uint32_t timescale = 30000;
    uint64_t duration = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    std::string codec;
    std::string handler_type;
    std::vector<uint32_t> keyframes;  // 1-based sample numbers

    bool open(const std::string& path) {
        file.open(path, std::ios::binary);
        if (!file) return false;
        file.seekg(0, std::ios::end);
        file_size = file.tellg();
        file.seekg(0);
        return true;
    }

    void close() { if (file.is_open()) file.close(); }

    std::shared_ptr<Box> read_box(uint64_t max_end) {
        uint64_t pos = file.tellg();
        if (pos >= max_end) return nullptr;

        uint8_t header[8];
        file.read(reinterpret_cast<char*>(header), 8);
        if (file.gcount() < 8) return nullptr;

        auto box = std::make_shared<Box>();
        box->offset = pos;
        box->size = read_be32(header);
        box->type = std::string(reinterpret_cast<char*>(header + 4), 4);
        box->header_size = 8;

        if (box->size == 1) {
            uint8_t ext[8];
            file.read(reinterpret_cast<char*>(ext), 8);
            box->size = read_be64(ext);
            box->header_size = 16;
        } else if (box->size == 0) {
            box->size = max_end - pos;
        }
        return box;
    }

    bool parse() {
        file.seekg(0);
        std::shared_ptr<Box> moov_box = nullptr;
        std::shared_ptr<Box> mdat_box = nullptr;

        // Find moov and mdat
        while ((uint64_t)file.tellg() < file_size) {
            auto box = read_box(file_size);
            if (!box) break;
            if (box->type == "moov") moov_box = box;
            if (box->type == "mdat") mdat_box = box;
            file.seekg(box->offset + box->size);
        }

        if (!moov_box || !mdat_box) return false;

        // Parse moov to find first track
        file.seekg(moov_box->offset + moov_box->header_size);
        while ((uint64_t)file.tellg() < moov_box->offset + moov_box->size) {
            auto box = read_box(moov_box->offset + moov_box->size);
            if (!box) break;
            if (box->type == "trak") {
                parse_trak(box);
                break;  // Only first track
            }
            file.seekg(box->offset + box->size);
        }

        // If duration is invalid, calculate from sample count
        if (!samples.empty() && (duration == 0 || duration > 0x7FFFFFFF)) {
            // For video: 30fps at 30000 timescale = 1001 per frame
            // For audio: 48kHz with 1024 samples per frame
            if (handler_type == "vide") {
                duration = samples.size() * 1001;
            } else if (handler_type == "soun") {
                duration = samples.size() * 1024;
            }
        }
        return !samples.empty();
    }

    void parse_trak(std::shared_ptr<Box> trak) {
        file.seekg(trak->offset + trak->header_size);
        uint64_t trak_end = trak->offset + trak->size;

        while ((uint64_t)file.tellg() < trak_end) {
            auto box = read_box(trak_end);
            if (!box) break;

            if (box->type == "tkhd") {
                file.seekg(box->offset + box->header_size);
                uint8_t buf[92];
                file.read(reinterpret_cast<char*>(buf), 92);
                uint8_t version = buf[0];
                if (version == 0) {
                    width = read_be32(buf + 76) >> 16;
                    height = read_be32(buf + 80) >> 16;
                }
            } else if (box->type == "mdia") {
                parse_mdia(box);
            }
            file.seekg(box->offset + box->size);
        }
    }

    void parse_mdia(std::shared_ptr<Box> mdia) {
        file.seekg(mdia->offset + mdia->header_size);
        uint64_t mdia_end = mdia->offset + mdia->size;

        while ((uint64_t)file.tellg() < mdia_end) {
            auto box = read_box(mdia_end);
            if (!box) break;

            if (box->type == "mdhd") {
                file.seekg(box->offset + box->header_size);
                uint8_t buf[32];
                file.read(reinterpret_cast<char*>(buf), 32);
                uint8_t version = buf[0];
                if (version == 0) {
                    uint32_t ts = read_be32(buf + 12);
                    // Validate timescale - if unreasonable, keep default
                    if (ts > 0 && ts < 1000000) {
                        timescale = ts;
                    }
                    duration = read_be32(buf + 16);
                }
            } else if (box->type == "hdlr") {
                file.seekg(box->offset + box->header_size + 8);
                char h[5] = {0};
                file.read(h, 4);
                handler_type = h;
            } else if (box->type == "minf") {
                parse_minf(box);
            }
            file.seekg(box->offset + box->size);
        }
    }

    void parse_minf(std::shared_ptr<Box> minf) {
        file.seekg(minf->offset + minf->header_size);
        uint64_t minf_end = minf->offset + minf->size;

        while ((uint64_t)file.tellg() < minf_end) {
            auto box = read_box(minf_end);
            if (!box) break;
            if (box->type == "stbl") {
                parse_stbl(box);
            }
            file.seekg(box->offset + box->size);
        }
    }

    void parse_stbl(std::shared_ptr<Box> stbl) {
        file.seekg(stbl->offset + stbl->header_size);
        uint64_t stbl_end = stbl->offset + stbl->size;

        std::vector<uint64_t> chunk_offsets;
        std::vector<uint32_t> sample_sizes;
        std::vector<uint32_t> first_chunk, samples_per_chunk;

        while ((uint64_t)file.tellg() < stbl_end) {
            auto box = read_box(stbl_end);
            if (!box) break;

            if (box->type == "stsd") {
                // Save stsd content
                stsd_data.resize(box->size - box->header_size);
                file.seekg(box->offset + box->header_size);
                file.read(reinterpret_cast<char*>(stsd_data.data()), stsd_data.size());
                // Get codec
                if (stsd_data.size() >= 12) {
                    codec = std::string(reinterpret_cast<char*>(&stsd_data[8]), 4);
                }
            } else if (box->type == "stco") {
                file.seekg(box->offset + box->header_size + 4);
                uint8_t buf[4];
                file.read(reinterpret_cast<char*>(buf), 4);
                uint32_t count = read_be32(buf);
                chunk_offsets.resize(count);
                for (uint32_t i = 0; i < count; i++) {
                    file.read(reinterpret_cast<char*>(buf), 4);
                    chunk_offsets[i] = read_be32(buf);
                }
            } else if (box->type == "co64") {
                file.seekg(box->offset + box->header_size + 4);
                uint8_t buf[8];
                file.read(reinterpret_cast<char*>(buf), 4);
                uint32_t count = read_be32(buf);
                chunk_offsets.resize(count);
                for (uint32_t i = 0; i < count; i++) {
                    file.read(reinterpret_cast<char*>(buf), 8);
                    chunk_offsets[i] = read_be64(buf);
                }
            } else if (box->type == "stsz") {
                file.seekg(box->offset + box->header_size + 4);
                uint8_t buf[8];
                file.read(reinterpret_cast<char*>(buf), 8);
                uint32_t sample_size = read_be32(buf);
                uint32_t count = read_be32(buf + 4);
                sample_sizes.resize(count);
                if (sample_size != 0) {
                    for (uint32_t i = 0; i < count; i++) sample_sizes[i] = sample_size;
                } else {
                    for (uint32_t i = 0; i < count; i++) {
                        file.read(reinterpret_cast<char*>(buf), 4);
                        sample_sizes[i] = read_be32(buf);
                    }
                }
            } else if (box->type == "stsc") {
                file.seekg(box->offset + box->header_size + 4);
                uint8_t buf[12];
                file.read(reinterpret_cast<char*>(buf), 4);
                uint32_t count = read_be32(buf);
                first_chunk.resize(count);
                samples_per_chunk.resize(count);
                for (uint32_t i = 0; i < count; i++) {
                    file.read(reinterpret_cast<char*>(buf), 12);
                    first_chunk[i] = read_be32(buf);
                    samples_per_chunk[i] = read_be32(buf + 4);
                }
            } else if (box->type == "stss") {
                file.seekg(box->offset + box->header_size + 4);
                uint8_t buf[4];
                file.read(reinterpret_cast<char*>(buf), 4);
                uint32_t count = read_be32(buf);
                keyframes.resize(count);
                for (uint32_t i = 0; i < count; i++) {
                    file.read(reinterpret_cast<char*>(buf), 4);
                    keyframes[i] = read_be32(buf);
                }
            }
            file.seekg(box->offset + box->size);
        }

        // Build sample list
        if (!chunk_offsets.empty() && !sample_sizes.empty()) {
            uint32_t sample_idx = 0;
            uint32_t stsc_idx = 0;

            for (uint32_t chunk_idx = 0; chunk_idx < chunk_offsets.size(); chunk_idx++) {
                while (stsc_idx + 1 < first_chunk.size() &&
                       first_chunk[stsc_idx + 1] <= chunk_idx + 1) {
                    stsc_idx++;
                }
                uint32_t samples_in_chunk = samples_per_chunk.empty() ? 1 : samples_per_chunk[stsc_idx];
                uint64_t offset = chunk_offsets[chunk_idx];

                for (uint32_t s = 0; s < samples_in_chunk && sample_idx < sample_sizes.size(); s++) {
                    Sample sample;
                    sample.offset = offset;
                    sample.size = sample_sizes[sample_idx];
                    samples.push_back(sample);
                    offset += sample.size;
                    sample_idx++;
                }
            }
        }
    }
};

// Read binary file into vector
std::vector<uint8_t> read_binary_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) return {};
    size_t size = f.tellg();
    f.seekg(0);
    std::vector<uint8_t> data(size);
    f.read(reinterpret_cast<char*>(data.data()), size);
    return data;
}

// Read sample sizes file (big-endian uint32 per sample)
std::vector<uint32_t> read_sizes_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) return {};
    size_t file_size = f.tellg();
    size_t num_samples = file_size / 4;
    f.seekg(0);
    std::vector<uint32_t> sizes(num_samples);
    for (size_t i = 0; i < num_samples; i++) {
        uint8_t buf[4];
        f.read(reinterpret_cast<char*>(buf), 4);
        sizes[i] = read_be32(buf);
    }
    return sizes;
}

int recompose(const std::string& base_name, const std::string& output_file) {
    std::cout << "Recomposing OSV from: " << base_name << std::endl;

    // Parse input MP4 files
    MP4Reader video1, video2, audio;

    std::string track1_path = base_name + "-track1.mp4";
    std::string track2_path = base_name + "-track2.mp4";
    std::string track3_path = base_name + "-track3.mp4";

    if (!video1.open(track1_path) || !video1.parse()) {
        std::cerr << "Error: Cannot parse " << track1_path << std::endl;
        return 1;
    }
    std::cout << "Video1: " << video1.samples.size() << " samples, "
              << video1.width << "x" << video1.height << std::endl;

    if (!video2.open(track2_path) || !video2.parse()) {
        std::cerr << "Error: Cannot parse " << track2_path << std::endl;
        return 1;
    }
    std::cout << "Video2: " << video2.samples.size() << " samples" << std::endl;

    if (!audio.open(track3_path) || !audio.parse()) {
        std::cerr << "Error: Cannot parse " << track3_path << std::endl;
        return 1;
    }
    std::cout << "Audio: " << audio.samples.size() << " samples" << std::endl;

    // Read metadata binary files
    std::vector<uint8_t> djmd1_data = read_binary_file(base_name + "-djmd1.bin");
    std::vector<uint8_t> djmd2_data = read_binary_file(base_name + "-djmd2.bin");
    std::vector<uint8_t> dbgi1_data = read_binary_file(base_name + "-dbgi1.bin");
    std::vector<uint8_t> dbgi2_data = read_binary_file(base_name + "-dbgi2.bin");

    if (djmd1_data.empty() || djmd2_data.empty() || dbgi1_data.empty() || dbgi2_data.empty()) {
        std::cerr << "Error: Cannot read metadata binary files" << std::endl;
        return 1;
    }

    std::cout << "DJMD1: " << djmd1_data.size() << " bytes" << std::endl;
    std::cout << "DJMD2: " << djmd2_data.size() << " bytes" << std::endl;
    std::cout << "DBGI1: " << dbgi1_data.size() << " bytes" << std::endl;
    std::cout << "DBGI2: " << dbgi2_data.size() << " bytes" << std::endl;

    // Read sample sizes files for metadata tracks
    std::vector<uint32_t> djmd1_sizes = read_sizes_file(base_name + "-djmd1.bin.sizes");
    std::vector<uint32_t> djmd2_sizes = read_sizes_file(base_name + "-djmd2.bin.sizes");
    std::vector<uint32_t> dbgi1_sizes = read_sizes_file(base_name + "-dbgi1.bin.sizes");
    std::vector<uint32_t> dbgi2_sizes = read_sizes_file(base_name + "-dbgi2.bin.sizes");

    if (djmd1_sizes.empty() || djmd2_sizes.empty() || dbgi1_sizes.empty() || dbgi2_sizes.empty()) {
        std::cerr << "Error: Cannot read sample sizes files (.sizes)" << std::endl;
        return 1;
    }

    std::cout << "DJMD1 sizes: " << djmd1_sizes.size() << " samples" << std::endl;
    std::cout << "DJMD2 sizes: " << djmd2_sizes.size() << " samples" << std::endl;
    std::cout << "DBGI1 sizes: " << dbgi1_sizes.size() << " samples" << std::endl;
    std::cout << "DBGI2 sizes: " << dbgi2_sizes.size() << " samples" << std::endl;

    // Read additional boxes
    std::ifstream boxes_file(base_name + "-additional-boxes.bin", std::ios::binary);
    if (!boxes_file) {
        std::cerr << "Error: Cannot open additional-boxes.bin" << std::endl;
        return 1;
    }

    // Parse additional boxes file
    char magic[4];
    boxes_file.read(magic, 4);
    if (memcmp(magic, "OSVT", 4) != 0) {
        std::cerr << "Error: Invalid additional-boxes file format" << std::endl;
        return 1;
    }

    uint8_t buf[8];
    boxes_file.read(reinterpret_cast<char*>(buf), 4);  // version

    // Read ftyp
    boxes_file.read(reinterpret_cast<char*>(buf), 4);
    uint32_t ftyp_size = read_be32(buf);
    std::vector<uint8_t> ftyp_data(ftyp_size);
    boxes_file.read(reinterpret_cast<char*>(ftyp_data.data()), ftyp_size);

    // Read free index
    boxes_file.read(reinterpret_cast<char*>(buf), 4);
    uint32_t free_index_size = read_be32(buf);
    std::vector<uint8_t> free_index_data(free_index_size);
    boxes_file.read(reinterpret_cast<char*>(free_index_data.data()), free_index_size);

    // Read saved boxes (udta, meta, camd)
    boxes_file.read(reinterpret_cast<char*>(buf), 4);
    uint32_t box_count = read_be32(buf);

    std::vector<std::vector<uint8_t>> saved_boxes;
    for (uint32_t i = 0; i < box_count; i++) {
        boxes_file.read(reinterpret_cast<char*>(buf), 8);
        uint64_t box_size = read_be64(buf);
        std::vector<uint8_t> box_data(box_size);
        boxes_file.read(reinterpret_cast<char*>(box_data.data()), box_size);
        saved_boxes.push_back(box_data);
    }
    boxes_file.close();

    std::cout << "Loaded " << saved_boxes.size() << " additional boxes" << std::endl;

    // Calculate total mdat size
    uint64_t video1_total = 0, video2_total = 0, audio_total = 0;
    for (const auto& s : video1.samples) video1_total += s.size;
    for (const auto& s : video2.samples) video2_total += s.size;
    for (const auto& s : audio.samples) audio_total += s.size;

    uint64_t mdat_content_size = video1_total + video2_total + audio_total +
                                  djmd1_data.size() + djmd2_data.size() +
                                  dbgi1_data.size() + dbgi2_data.size();

    std::cout << "Total mdat content: " << mdat_content_size << " bytes" << std::endl;

    // Open output file
    std::ofstream out(output_file, std::ios::binary);
    if (!out) {
        std::cerr << "Error: Cannot create output file" << std::endl;
        return 1;
    }

    // Write ftyp
    out.write(reinterpret_cast<char*>(ftyp_data.data()), ftyp_data.size());

    // Write empty free (8 bytes)
    uint8_t free_8[8] = {0, 0, 0, 8, 'f', 'r', 'e', 'e'};
    out.write(reinterpret_cast<char*>(free_8), 8);

    // Write free index (will update offsets later)
    uint64_t free_index_pos = out.tellp();
    out.write(reinterpret_cast<char*>(free_index_data.data()), free_index_data.size());

    // Calculate mdat position
    uint64_t mdat_start = out.tellp();
    bool use_64bit_mdat = (mdat_content_size + 8) > 0xFFFFFFFF;
    uint64_t mdat_header_size = use_64bit_mdat ? 16 : 8;

    // Write mdat header
    if (use_64bit_mdat) {
        uint8_t mdat_header[16];
        write_be32(mdat_header, 1);
        memcpy(mdat_header + 4, "mdat", 4);
        write_be64(mdat_header + 8, mdat_header_size + mdat_content_size);
        out.write(reinterpret_cast<char*>(mdat_header), 16);
    } else {
        uint8_t mdat_header[8];
        write_be32(mdat_header, 8 + mdat_content_size);
        memcpy(mdat_header + 4, "mdat", 4);
        out.write(reinterpret_cast<char*>(mdat_header), 8);
    }

    uint64_t data_offset = mdat_start + mdat_header_size;

    // Write video1 samples and record offsets
    std::vector<uint64_t> video1_offsets;
    std::vector<uint8_t> copy_buf;
    for (const auto& sample : video1.samples) {
        video1_offsets.push_back(data_offset);
        copy_buf.resize(sample.size);
        video1.file.seekg(sample.offset);
        video1.file.read(reinterpret_cast<char*>(copy_buf.data()), sample.size);
        out.write(reinterpret_cast<char*>(copy_buf.data()), sample.size);
        data_offset += sample.size;
    }

    // Write video2 samples
    std::vector<uint64_t> video2_offsets;
    for (const auto& sample : video2.samples) {
        video2_offsets.push_back(data_offset);
        copy_buf.resize(sample.size);
        video2.file.seekg(sample.offset);
        video2.file.read(reinterpret_cast<char*>(copy_buf.data()), sample.size);
        out.write(reinterpret_cast<char*>(copy_buf.data()), sample.size);
        data_offset += sample.size;
    }

    // Write audio samples
    std::vector<uint64_t> audio_offsets;
    for (const auto& sample : audio.samples) {
        audio_offsets.push_back(data_offset);
        copy_buf.resize(sample.size);
        audio.file.seekg(sample.offset);
        audio.file.read(reinterpret_cast<char*>(copy_buf.data()), sample.size);
        out.write(reinterpret_cast<char*>(copy_buf.data()), sample.size);
        data_offset += sample.size;
    }

    // Write djmd1
    uint64_t djmd1_offset = data_offset;
    out.write(reinterpret_cast<char*>(djmd1_data.data()), djmd1_data.size());
    data_offset += djmd1_data.size();

    // Write djmd2
    uint64_t djmd2_offset = data_offset;
    out.write(reinterpret_cast<char*>(djmd2_data.data()), djmd2_data.size());
    data_offset += djmd2_data.size();

    // Write dbgi1
    uint64_t dbgi1_offset = data_offset;
    out.write(reinterpret_cast<char*>(dbgi1_data.data()), dbgi1_data.size());
    data_offset += dbgi1_data.size();

    // Write dbgi2
    uint64_t dbgi2_offset = data_offset;
    out.write(reinterpret_cast<char*>(dbgi2_data.data()), dbgi2_data.size());
    data_offset += dbgi2_data.size();

    uint64_t moov_start = out.tellp();

    // Build moov
    std::vector<uint8_t> moov_data;
    auto append = [&moov_data](const void* data, size_t size) {
        const uint8_t* p = static_cast<const uint8_t*>(data);
        moov_data.insert(moov_data.end(), p, p + size);
    };

    moov_data.resize(8);  // Reserve for moov header

    // mvhd
    {
        uint8_t mvhd[108];
        memset(mvhd, 0, sizeof(mvhd));
        write_be32(mvhd, 108);
        memcpy(mvhd + 4, "mvhd", 4);
        write_be32(mvhd + 20, video1.timescale);
        write_be32(mvhd + 24, video1.duration);
        write_be32(mvhd + 28, 0x00010000);  // rate 1.0
        mvhd[32] = 0x01;  // volume 1.0
        write_be32(mvhd + 44, 0x00010000);  // matrix
        write_be32(mvhd + 60, 0x00010000);
        write_be32(mvhd + 84, 0x40000000);
        write_be32(mvhd + 104, 8);  // next track id
        append(mvhd, 108);
    }

    // Helper to build a video track
    auto build_video_trak = [&](uint32_t track_id, const MP4Reader& vid,
                                const std::vector<uint64_t>& offsets) {
        size_t trak_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // tkhd
        {
            uint8_t tkhd[92];
            memset(tkhd, 0, sizeof(tkhd));
            write_be32(tkhd, 92);
            memcpy(tkhd + 4, "tkhd", 4);
            write_be32(tkhd + 8, track_id == 1 ? 0x00000003 : 0x00000000);  // flags
            write_be32(tkhd + 20, track_id);
            write_be32(tkhd + 28, vid.duration);
            write_be32(tkhd + 48, 0x00010000);  // matrix
            write_be32(tkhd + 64, 0x00010000);
            write_be32(tkhd + 88, 0x40000000);
            write_be32(tkhd + 84, vid.width << 16);
            write_be32(tkhd + 88, vid.height << 16);
            append(tkhd, 92);
        }

        // edts/elst
        {
            uint8_t edts[36];
            write_be32(edts, 36);
            memcpy(edts + 4, "edts", 4);
            write_be32(edts + 8, 28);
            memcpy(edts + 12, "elst", 4);
            memset(edts + 16, 0, 4);
            write_be32(edts + 20, 1);
            write_be32(edts + 24, vid.duration);
            write_be32(edts + 28, 0);
            write_be32(edts + 32, 0x00010000);
            append(edts, 36);
        }

        // mdia
        size_t mdia_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // mdhd (version 0: 8 header + 4 ver/flags + 4 create + 4 modify + 4 timescale + 4 duration + 2 lang + 2 quality = 32)
        {
            uint8_t mdhd[32];
            memset(mdhd, 0, sizeof(mdhd));
            write_be32(mdhd, 32);
            memcpy(mdhd + 4, "mdhd", 4);
            // version/flags at 8 (already 0)
            // creation_time at 12 (already 0)
            // modification_time at 16 (already 0)
            write_be32(mdhd + 20, vid.timescale);
            write_be32(mdhd + 24, vid.duration);
            mdhd[28] = 0x55; mdhd[29] = 0xC4;  // language "und"
            append(mdhd, 32);
        }

        // hdlr
        {
            const char* name = "VideoHandler";
            size_t name_len = strlen(name) + 1;
            size_t hdlr_size = 32 + name_len;
            uint8_t hdlr[32];
            memset(hdlr, 0, sizeof(hdlr));
            write_be32(hdlr, hdlr_size);
            memcpy(hdlr + 4, "hdlr", 4);
            memcpy(hdlr + 16, "vide", 4);
            append(hdlr, 32);
            moov_data.insert(moov_data.end(), name, name + name_len);
        }

        // minf
        size_t minf_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // vmhd
        {
            uint8_t vmhd[20];
            memset(vmhd, 0, sizeof(vmhd));
            write_be32(vmhd, 20);
            memcpy(vmhd + 4, "vmhd", 4);
            write_be32(vmhd + 8, 0x00000001);
            append(vmhd, 20);
        }

        // dinf/dref
        {
            uint8_t dinf[36];
            memset(dinf, 0, sizeof(dinf));
            write_be32(dinf, 36);
            memcpy(dinf + 4, "dinf", 4);
            write_be32(dinf + 8, 28);
            memcpy(dinf + 12, "dref", 4);
            write_be32(dinf + 20, 1);
            write_be32(dinf + 24, 12);
            memcpy(dinf + 28, "url ", 4);
            write_be32(dinf + 32, 0x00000001);
            append(dinf, 36);
        }

        // stbl
        size_t stbl_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // stsd
        {
            size_t stsd_size = 8 + vid.stsd_data.size();
            uint8_t stsd_header[8];
            write_be32(stsd_header, stsd_size);
            memcpy(stsd_header + 4, "stsd", 4);
            append(stsd_header, 8);
            append(vid.stsd_data.data(), vid.stsd_data.size());
        }

        // stts
        {
            uint8_t stts[24];
            memset(stts, 0, sizeof(stts));
            write_be32(stts, 24);
            memcpy(stts + 4, "stts", 4);
            write_be32(stts + 12, 1);
            write_be32(stts + 16, vid.samples.size());
            write_be32(stts + 20, 1001);  // 30fps at 30000 timescale
            append(stts, 24);
        }

        // stss (keyframes)
        if (!vid.keyframes.empty()) {
            size_t stss_size = 16 + vid.keyframes.size() * 4;
            uint8_t stss_header[16];
            write_be32(stss_header, stss_size);
            memcpy(stss_header + 4, "stss", 4);
            memset(stss_header + 8, 0, 4);
            write_be32(stss_header + 12, vid.keyframes.size());
            append(stss_header, 16);
            for (uint32_t kf : vid.keyframes) {
                uint8_t b[4];
                write_be32(b, kf);
                append(b, 4);
            }
        }

        // stsc
        {
            uint8_t stsc[28];
            memset(stsc, 0, sizeof(stsc));
            write_be32(stsc, 28);
            memcpy(stsc + 4, "stsc", 4);
            write_be32(stsc + 12, 1);
            write_be32(stsc + 16, 1);
            write_be32(stsc + 20, 1);
            write_be32(stsc + 24, 1);
            append(stsc, 28);
        }

        // stsz
        {
            size_t stsz_size = 20 + vid.samples.size() * 4;
            uint8_t stsz_header[20];
            write_be32(stsz_header, stsz_size);
            memcpy(stsz_header + 4, "stsz", 4);
            memset(stsz_header + 8, 0, 8);
            write_be32(stsz_header + 16, vid.samples.size());
            append(stsz_header, 20);
            for (const auto& s : vid.samples) {
                uint8_t b[4];
                write_be32(b, s.size);
                append(b, 4);
            }
        }

        // stco
        {
            size_t stco_size = 16 + offsets.size() * 4;
            uint8_t stco_header[16];
            write_be32(stco_header, stco_size);
            memcpy(stco_header + 4, "stco", 4);
            memset(stco_header + 8, 0, 4);
            write_be32(stco_header + 12, offsets.size());
            append(stco_header, 16);
            for (uint64_t off : offsets) {
                uint8_t b[4];
                write_be32(b, off);
                append(b, 4);
            }
        }

        // Fix sizes
        write_be32(&moov_data[stbl_start], moov_data.size() - stbl_start);
        memcpy(&moov_data[stbl_start + 4], "stbl", 4);
        write_be32(&moov_data[minf_start], moov_data.size() - minf_start);
        memcpy(&moov_data[minf_start + 4], "minf", 4);
        write_be32(&moov_data[mdia_start], moov_data.size() - mdia_start);
        memcpy(&moov_data[mdia_start + 4], "mdia", 4);
        write_be32(&moov_data[trak_start], moov_data.size() - trak_start);
        memcpy(&moov_data[trak_start + 4], "trak", 4);
    };

    // Helper to build audio track
    auto build_audio_trak = [&](uint32_t track_id, const MP4Reader& aud,
                                const std::vector<uint64_t>& offsets) {
        size_t trak_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // tkhd
        {
            uint8_t tkhd[92];
            memset(tkhd, 0, sizeof(tkhd));
            write_be32(tkhd, 92);
            memcpy(tkhd + 4, "tkhd", 4);
            write_be32(tkhd + 8, 0x00000003);
            write_be32(tkhd + 20, track_id);
            write_be32(tkhd + 28, aud.duration);
            write_be32(tkhd + 36, 0x01000000);  // volume
            write_be32(tkhd + 48, 0x00010000);
            write_be32(tkhd + 64, 0x00010000);
            write_be32(tkhd + 88, 0x40000000);
            append(tkhd, 92);
        }

        // edts/elst
        {
            uint8_t edts[36];
            write_be32(edts, 36);
            memcpy(edts + 4, "edts", 4);
            write_be32(edts + 8, 28);
            memcpy(edts + 12, "elst", 4);
            memset(edts + 16, 0, 4);
            write_be32(edts + 20, 1);
            write_be32(edts + 24, aud.duration);
            write_be32(edts + 28, 0);
            write_be32(edts + 32, 0x00010000);
            append(edts, 36);
        }

        // mdia
        size_t mdia_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // mdhd
        {
            uint8_t mdhd[32];
            memset(mdhd, 0, sizeof(mdhd));
            write_be32(mdhd, 32);
            memcpy(mdhd + 4, "mdhd", 4);
            write_be32(mdhd + 20, aud.timescale);
            write_be32(mdhd + 24, aud.duration);
            mdhd[28] = 0x55; mdhd[29] = 0xC4;
            append(mdhd, 32);
        }

        // hdlr
        {
            const char* name = "SoundHandler";
            size_t name_len = strlen(name) + 1;
            size_t hdlr_size = 32 + name_len;
            uint8_t hdlr[32];
            memset(hdlr, 0, sizeof(hdlr));
            write_be32(hdlr, hdlr_size);
            memcpy(hdlr + 4, "hdlr", 4);
            memcpy(hdlr + 16, "soun", 4);
            append(hdlr, 32);
            moov_data.insert(moov_data.end(), name, name + name_len);
        }

        // minf
        size_t minf_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // smhd
        {
            uint8_t smhd[16];
            memset(smhd, 0, sizeof(smhd));
            write_be32(smhd, 16);
            memcpy(smhd + 4, "smhd", 4);
            append(smhd, 16);
        }

        // dinf/dref
        {
            uint8_t dinf[36];
            memset(dinf, 0, sizeof(dinf));
            write_be32(dinf, 36);
            memcpy(dinf + 4, "dinf", 4);
            write_be32(dinf + 8, 28);
            memcpy(dinf + 12, "dref", 4);
            write_be32(dinf + 20, 1);
            write_be32(dinf + 24, 12);
            memcpy(dinf + 28, "url ", 4);
            write_be32(dinf + 32, 0x00000001);
            append(dinf, 36);
        }

        // stbl
        size_t stbl_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // stsd
        {
            size_t stsd_size = 8 + aud.stsd_data.size();
            uint8_t stsd_header[8];
            write_be32(stsd_header, stsd_size);
            memcpy(stsd_header + 4, "stsd", 4);
            append(stsd_header, 8);
            append(aud.stsd_data.data(), aud.stsd_data.size());
        }

        // stts
        {
            uint8_t stts[24];
            memset(stts, 0, sizeof(stts));
            write_be32(stts, 24);
            memcpy(stts + 4, "stts", 4);
            write_be32(stts + 12, 1);
            write_be32(stts + 16, aud.samples.size());
            write_be32(stts + 20, 1024);  // AAC frame size
            append(stts, 24);
        }

        // stsc - audio typically has different chunk layout
        {
            uint8_t stsc[28];
            memset(stsc, 0, sizeof(stsc));
            write_be32(stsc, 28);
            memcpy(stsc + 4, "stsc", 4);
            write_be32(stsc + 12, 1);
            write_be32(stsc + 16, 1);
            write_be32(stsc + 20, 1);
            write_be32(stsc + 24, 1);
            append(stsc, 28);
        }

        // stsz
        {
            size_t stsz_size = 20 + aud.samples.size() * 4;
            uint8_t stsz_header[20];
            write_be32(stsz_header, stsz_size);
            memcpy(stsz_header + 4, "stsz", 4);
            memset(stsz_header + 8, 0, 8);
            write_be32(stsz_header + 16, aud.samples.size());
            append(stsz_header, 20);
            for (const auto& s : aud.samples) {
                uint8_t b[4];
                write_be32(b, s.size);
                append(b, 4);
            }
        }

        // stco
        {
            size_t stco_size = 16 + offsets.size() * 4;
            uint8_t stco_header[16];
            write_be32(stco_header, stco_size);
            memcpy(stco_header + 4, "stco", 4);
            memset(stco_header + 8, 0, 4);
            write_be32(stco_header + 12, offsets.size());
            append(stco_header, 16);
            for (uint64_t off : offsets) {
                uint8_t b[4];
                write_be32(b, off);
                append(b, 4);
            }
        }

        // Fix sizes
        write_be32(&moov_data[stbl_start], moov_data.size() - stbl_start);
        memcpy(&moov_data[stbl_start + 4], "stbl", 4);
        write_be32(&moov_data[minf_start], moov_data.size() - minf_start);
        memcpy(&moov_data[minf_start + 4], "minf", 4);
        write_be32(&moov_data[mdia_start], moov_data.size() - mdia_start);
        memcpy(&moov_data[mdia_start + 4], "mdia", 4);
        write_be32(&moov_data[trak_start], moov_data.size() - trak_start);
        memcpy(&moov_data[trak_start + 4], "trak", 4);
    };

    // Helper to build metadata track (djmd/dbgi)
    auto build_meta_trak = [&](uint32_t track_id, const char* codec_name,
                               const char* handler_name,
                               const std::vector<uint8_t>& data,
                               uint64_t data_offset_in_mdat,
                               const std::vector<uint32_t>& sample_sizes,
                               uint32_t sample_duration) {
        size_t trak_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        uint32_t num_samples = sample_sizes.size();

        // tkhd
        {
            uint8_t tkhd[92];
            memset(tkhd, 0, sizeof(tkhd));
            write_be32(tkhd, 92);
            memcpy(tkhd + 4, "tkhd", 4);
            write_be32(tkhd + 20, track_id);
            write_be32(tkhd + 28, num_samples * sample_duration);
            write_be32(tkhd + 48, 0x00010000);
            write_be32(tkhd + 64, 0x00010000);
            write_be32(tkhd + 88, 0x40000000);
            append(tkhd, 92);
        }

        // edts/elst
        {
            uint8_t edts[36];
            write_be32(edts, 36);
            memcpy(edts + 4, "edts", 4);
            write_be32(edts + 8, 28);
            memcpy(edts + 12, "elst", 4);
            memset(edts + 16, 0, 4);
            write_be32(edts + 20, 1);
            write_be32(edts + 24, num_samples * sample_duration);
            write_be32(edts + 28, 0);
            write_be32(edts + 32, 0x00010000);
            append(edts, 36);
        }

        // mdia
        size_t mdia_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // mdhd
        {
            uint8_t mdhd[32];
            memset(mdhd, 0, sizeof(mdhd));
            write_be32(mdhd, 32);
            memcpy(mdhd + 4, "mdhd", 4);
            write_be32(mdhd + 20, 30000);  // timescale
            write_be32(mdhd + 24, num_samples * sample_duration);
            mdhd[28] = 0x55; mdhd[29] = 0xC4;
            append(mdhd, 32);
        }

        // hdlr
        {
            size_t name_len = strlen(handler_name) + 1;
            size_t hdlr_size = 32 + name_len;
            uint8_t hdlr[32];
            memset(hdlr, 0, sizeof(hdlr));
            write_be32(hdlr, hdlr_size);
            memcpy(hdlr + 4, "hdlr", 4);
            memcpy(hdlr + 16, "meta", 4);
            append(hdlr, 32);
            moov_data.insert(moov_data.end(), handler_name, handler_name + name_len);
        }

        // minf
        size_t minf_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // gmhd (for metadata tracks)
        bool is_dbgi = (strcmp(codec_name, "dbgi") == 0);
        {
            size_t gmhd_size = is_dbgi ? 88 : 84;
            std::vector<uint8_t> gmhd(gmhd_size, 0);
            write_be32(gmhd.data(), gmhd_size);
            memcpy(gmhd.data() + 4, "gmhd", 4);
            // gmin subbox
            write_be32(gmhd.data() + 8, 24);
            memcpy(gmhd.data() + 12, "gmin", 4);
            // text subbox
            size_t text_size = is_dbgi ? 52 : 48;
            write_be32(gmhd.data() + 32, text_size);
            memcpy(gmhd.data() + 36, "text", 4);
            append(gmhd.data(), gmhd_size);
        }

        // dinf/dref
        {
            uint8_t dinf[36];
            memset(dinf, 0, sizeof(dinf));
            write_be32(dinf, 36);
            memcpy(dinf + 4, "dinf", 4);
            write_be32(dinf + 8, 28);
            memcpy(dinf + 12, "dref", 4);
            write_be32(dinf + 20, 1);
            write_be32(dinf + 24, 12);
            memcpy(dinf + 28, "url ", 4);
            write_be32(dinf + 32, 0x00000001);
            append(dinf, 36);
        }

        // stbl
        size_t stbl_start = moov_data.size();
        moov_data.resize(moov_data.size() + 8);

        // stsd with codec-specific entry
        {
            uint8_t stsd[36];
            memset(stsd, 0, sizeof(stsd));
            write_be32(stsd, 36);
            memcpy(stsd + 4, "stsd", 4);
            write_be32(stsd + 12, 1);  // entry count
            write_be32(stsd + 16, 20);  // entry size
            memcpy(stsd + 20, codec_name, 4);
            write_be32(stsd + 32, 1);  // data reference index
            append(stsd, 36);
        }

        // stts
        {
            uint8_t stts[24];
            memset(stts, 0, sizeof(stts));
            write_be32(stts, 24);
            memcpy(stts + 4, "stts", 4);
            write_be32(stts + 12, 1);
            write_be32(stts + 16, num_samples);
            write_be32(stts + 20, sample_duration);
            append(stts, 24);
        }

        // stsc
        {
            uint8_t stsc[28];
            memset(stsc, 0, sizeof(stsc));
            write_be32(stsc, 28);
            memcpy(stsc + 4, "stsc", 4);
            write_be32(stsc + 12, 1);
            write_be32(stsc + 16, 1);
            write_be32(stsc + 20, 1);
            write_be32(stsc + 24, 1);
            append(stsc, 28);
        }

        // stsz - use actual sample sizes from sizes file
        {
            size_t stsz_size = 20 + num_samples * 4;
            uint8_t stsz_header[20];
            write_be32(stsz_header, stsz_size);
            memcpy(stsz_header + 4, "stsz", 4);
            memset(stsz_header + 8, 0, 8);  // sample_size = 0 means variable
            write_be32(stsz_header + 16, num_samples);
            append(stsz_header, 20);
            // Write actual per-sample sizes
            for (uint32_t i = 0; i < num_samples; i++) {
                uint8_t b[4];
                write_be32(b, sample_sizes[i]);
                append(b, 4);
            }
        }

        // stco - use actual sample sizes for offset calculation
        {
            size_t stco_size = 16 + num_samples * 4;
            uint8_t stco_header[16];
            write_be32(stco_header, stco_size);
            memcpy(stco_header + 4, "stco", 4);
            memset(stco_header + 8, 0, 4);
            write_be32(stco_header + 12, num_samples);
            append(stco_header, 16);
            uint64_t off = data_offset_in_mdat;
            for (uint32_t i = 0; i < num_samples; i++) {
                uint8_t b[4];
                write_be32(b, off);
                append(b, 4);
                off += sample_sizes[i];  // Use actual sample size
            }
        }

        // Fix sizes
        write_be32(&moov_data[stbl_start], moov_data.size() - stbl_start);
        memcpy(&moov_data[stbl_start + 4], "stbl", 4);
        write_be32(&moov_data[minf_start], moov_data.size() - minf_start);
        memcpy(&moov_data[minf_start + 4], "minf", 4);
        write_be32(&moov_data[mdia_start], moov_data.size() - mdia_start);
        memcpy(&moov_data[mdia_start + 4], "mdia", 4);
        write_be32(&moov_data[trak_start], moov_data.size() - trak_start);
        memcpy(&moov_data[trak_start + 4], "trak", 4);
    };

    // Build all tracks
    uint32_t num_video_samples = video1.samples.size();
    uint32_t sample_duration = 1001;  // for 30fps at 30000 timescale

    build_video_trak(1, video1, video1_offsets);
    build_video_trak(2, video2, video2_offsets);
    build_audio_trak(3, audio, audio_offsets);
    build_meta_trak(4, "djmd", "CAM meta", djmd1_data, djmd1_offset, djmd1_sizes, sample_duration);
    build_meta_trak(5, "djmd", "CAM meta", djmd2_data, djmd2_offset, djmd2_sizes, sample_duration);
    build_meta_trak(6, "dbgi", "CAM dbgi", dbgi1_data, dbgi1_offset, dbgi1_sizes, sample_duration);
    build_meta_trak(7, "dbgi", "CAM dbgi", dbgi2_data, dbgi2_offset, dbgi2_sizes, sample_duration);

    // Add udta from saved boxes (first box should be udta)
    if (saved_boxes.size() >= 1) {
        append(saved_boxes[0].data(), saved_boxes[0].size());
    }

    // Add meta from saved boxes (second box)
    if (saved_boxes.size() >= 2) {
        append(saved_boxes[1].data(), saved_boxes[1].size());
    }

    // Fix moov size
    write_be32(&moov_data[0], moov_data.size());
    memcpy(&moov_data[4], "moov", 4);

    // Write moov
    out.write(reinterpret_cast<char*>(moov_data.data()), moov_data.size());

    // Write camd (third saved box)
    if (saved_boxes.size() >= 3) {
        out.write(reinterpret_cast<char*>(saved_boxes[2].data()), saved_boxes[2].size());
    }

    out.close();

    video1.close();
    video2.close();
    audio.close();

    std::cout << "\nRecompose complete: " << output_file << std::endl;
    return 0;
}

void print_usage(const char* prog) {
    std::cout << "OSV Toolbox - DJI OSV File Extractor/Recomposer\n\n"
              << "Usage:\n"
              << "  " << prog << " --extract-data <input.OSV>\n"
              << "  " << prog << " --recompose <base-name> <output.OSV>\n\n"
              << "Operations:\n"
              << "  --extract-data  Extract video/audio tracks, djmd/dbgi samples, and metadata\n"
              << "  --recompose     Rebuild OSV file from extracted components\n\n"
              << "Extract outputs:\n"
              << "  <base>-track1.mp4        First video track\n"
              << "  <base>-track2.mp4        Second video track\n"
              << "  <base>-track3.mp4        Audio track (AAC in MP4 container)\n"
              << "  <base>-djmd1.bin         First DJMD metadata samples\n"
              << "  <base>-djmd2.bin         Second DJMD metadata samples\n"
              << "  <base>-dbgi1.bin         First DBGI metadata samples\n"
              << "  <base>-dbgi2.bin         Second DBGI metadata samples\n"
              << "  <base>-additional-boxes.bin  udta, meta, camd boxes\n";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }

    std::string operation = argv[1];

    if (operation == "--extract-data") {
        if (argc < 3) {
            std::cerr << "Error: Missing input file\n";
            print_usage(argv[0]);
            return 1;
        }
        return extract_data(argv[2]);
    } else if (operation == "--recompose") {
        if (argc < 4) {
            std::cerr << "Error: Missing arguments for recompose\n";
            print_usage(argv[0]);
            return 1;
        }
        return recompose(argv[2], argv[3]);
    } else {
        std::cerr << "Unknown operation: " << operation << "\n";
        print_usage(argv[0]);
        return 1;
    }

    return 0;
}
