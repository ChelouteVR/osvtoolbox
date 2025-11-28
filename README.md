# OSVToolbox

A cross-platform command-line tool for extracting and recomposing DJI Osmo 360 video files (`.OSV`).

## Overview

DJI Osmo 360 cameras produce `.OSV` files, which are MP4 containers with a custom structure containing:
- 2 HEVC video tracks (3840x3840 each, for front and back lenses)
- 1 AAC audio track
- 2 `djmd` metadata tracks (DJI camera metadata)
- 2 `dbgi` metadata tracks (DJI debug info)
- A `camd` box with additional camera data
- Embedded thumbnail

OSVToolbox allows you to:
1. **Extract** all components from an OSV file for processing
2. **Recompose** a valid OSV file from modified components

This is useful for video processing workflows where you need to modify the video streams while preserving all DJI metadata.

## Building

### Requirements

- C++17 compatible compiler
- FFmpeg (for video/audio extraction)

### macOS

```bash
g++ -std=c++17 -O2 -o osvtoolbox osvtoolbox.cpp
```

Or with Clang:
```bash
clang++ -std=c++17 -O2 -o osvtoolbox osvtoolbox.cpp
```

### Linux

```bash
g++ -std=c++17 -O2 -o osvtoolbox osvtoolbox.cpp
```

### Windows (MSVC)

```cmd
cl /std:c++17 /O2 /EHsc osvtoolbox.cpp /Fe:osvtoolbox.exe
```

Or with MinGW:
```bash
g++ -std=c++17 -O2 -o osvtoolbox.exe osvtoolbox.cpp
```

## Usage

### Extract Data

Extract all components from an OSV file:

```bash
./osvtoolbox --extract-data <input.OSV>
```

This creates the following files in the same directory as the input:

| File | Description |
|------|-------------|
| `<basename>-track1.mp4` | First video track (HEVC, 3840x3840) |
| `<basename>-track2.mp4` | Second video track (HEVC, 3840x3840) |
| `<basename>-track3.mp4` | Audio track (AAC, 48kHz stereo) |
| `<basename>-djmd1.bin` | First DJI metadata track (raw samples) |
| `<basename>-djmd2.bin` | Second DJI metadata track (raw samples) |
| `<basename>-dbgi1.bin` | First DJI debug track (raw samples) |
| `<basename>-dbgi2.bin` | Second DJI debug track (raw samples) |
| `<basename>-djmd1.bin.sizes` | Sample sizes for djmd1 |
| `<basename>-djmd2.bin.sizes` | Sample sizes for djmd2 |
| `<basename>-dbgi1.bin.sizes` | Sample sizes for dbgi1 |
| `<basename>-dbgi2.bin.sizes` | Sample sizes for dbgi2 |
| `<basename>-additional-boxes.bin` | ftyp, udta, meta, camd boxes |

**Example:**
```bash
./osvtoolbox --extract-data /path/to/VIDEO.OSV
```

### Recompose

Rebuild an OSV file from extracted (and optionally modified) components:

```bash
./osvtoolbox --recompose <basename> <output.OSV>
```

The tool expects all extracted files to be present with the same basename.

**Example:**
```bash
./osvtoolbox --recompose /path/to/VIDEO /path/to/OUTPUT.OSV
```

## Workflow Example

A typical workflow for processing OSV files:

```bash
# 1. Extract all data
./osvtoolbox --extract-data input.OSV

# 2. Process video tracks with your tools (e.g., stabilization, color grading)
#    The video files are standard MP4/HEVC and can be processed with any video tool
ffmpeg -i input-track1.mp4 -vf "your_filter" -c:v libx265 output-track1.mp4
ffmpeg -i input-track2.mp4 -vf "your_filter" -c:v libx265 output-track2.mp4

# 3. Replace the original track files with processed ones
mv output-track1.mp4 input-track1.mp4
mv output-track2.mp4 input-track2.mp4

# 4. Recompose the OSV file
./osvtoolbox --recompose input output.OSV
```

## Technical Details

### OSV File Structure

```
[ftyp]          - File type box (isom/mp41)
[free]          - Empty padding
[free]          - Index with thumbnail offset
[mdat]          - Media data (all samples)
  ├── Video 1 samples
  ├── Video 2 samples
  ├── Audio samples
  ├── DJMD 1 samples
  ├── DJMD 2 samples
  ├── DBGI 1 samples
  └── DBGI 2 samples
[moov]          - Movie metadata
  ├── [mvhd]    - Movie header
  ├── [trak]    - Track 1 (video)
  ├── [trak]    - Track 2 (video)
  ├── [trak]    - Track 3 (audio)
  ├── [trak]    - Track 4 (djmd)
  ├── [trak]    - Track 5 (djmd)
  ├── [trak]    - Track 6 (dbgi)
  ├── [trak]    - Track 7 (dbgi)
  ├── [udta]    - User data (thumbnail)
  └── [meta]    - Metadata
[camd]          - DJI camera metadata box
```

### Sample Sizes

The `djmd` and `dbgi` tracks have variable sample sizes. The `.sizes` files store the size of each sample as big-endian 32-bit integers, which is essential for correct recomposition.

### Metadata Preservation

The tool preserves:
- All DJI proprietary metadata (`djmd`, `dbgi`, `camd`)
- Embedded thumbnail image
- File type and compatibility information
- User data boxes

## Limitations

- Video/audio extraction uses FFmpeg (must be installed and in PATH)
- The recomposed file may have minor size differences due to box structure variations
- Does not support extended size boxes (files > 4GB per box)

## Dependencies

- **FFmpeg**: Required for extracting video and audio tracks as valid MP4 files

Install FFmpeg:
- macOS: `brew install ffmpeg`
- Ubuntu/Debian: `sudo apt install ffmpeg`
- Windows: Download from [ffmpeg.org](https://ffmpeg.org/download.html)

## License

MIT License

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.
