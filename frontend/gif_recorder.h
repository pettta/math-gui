#pragma once

#include <cstdint>
#include <cstdio>
#include <string>

// Records the app window to a GIF by piping raw BGRA frames into ffmpeg, then
// copies the finished GIF to the macOS clipboard. macOS-only (implemented in
// gif_recorder.mm). All ffmpeg/Cocoa specifics are hidden behind this interface
// so it can be used from the platform-agnostic engine loop.
class GifRecorder {
public:
    // Opens an ffmpeg pipe sized to width x height at the given fps. Returns
    // false (and sets status()) if ffmpeg is unavailable or the pipe fails.
    bool start(int width, int height, int fps);

    // Writes one raw BGRA frame (width*height*4 bytes). Frames whose dimensions
    // differ from the size locked in at start() are skipped.
    void addFrame(const uint8_t* bgra, int width, int height);

    // Closes the pipe, lets ffmpeg finalize the GIF, then copies it to the
    // clipboard. Returns true on success; sets status() either way.
    bool stop();

    bool isRecording() const { return pipe_ != nullptr; }
    const std::string& status() const { return status_; }

private:
    FILE* pipe_ = nullptr;
    int w_ = 0;
    int h_ = 0;
    int fps_ = 15;
    std::string outputPath_;
    std::string status_;
};
