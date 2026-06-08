#include "gif_recorder.h"

#include <cstdlib>
#include <cstdio>
#include <filesystem>

#ifdef _WIN32
    // On Windows the C runtime spells these with a leading underscore, and the
    // pipe MUST be opened in binary mode or raw frame bytes get mangled by
    // CRLF translation.
    #define GIF_POPEN  _popen
    #define GIF_PCLOSE _pclose
    #define GIF_PIPE_MODE "wb"
#else
    #define GIF_POPEN  popen
    #define GIF_PCLOSE pclose
    #define GIF_PIPE_MODE "w"
#endif

namespace {
// Returns true if an ffmpeg executable can be found on PATH.
bool ffmpegAvailable()
{
#ifdef _WIN32
    return std::system("where ffmpeg >nul 2>nul") == 0;
#else
    return std::system("command -v ffmpeg >/dev/null 2>&1") == 0;
#endif
}

// A writable temp path for the encoded GIF (e.g. /tmp or %TEMP%).
std::string defaultOutputPath()
{
    std::error_code ec;
    std::filesystem::path dir = std::filesystem::temp_directory_path(ec);
    if (ec)
    {
        // Fall back to the current directory if the temp dir is unknown.
        dir = std::filesystem::path(".");
    }
    return (dir / "mathgui_recording.gif").string();
}
}

bool GifRecorder::start(int width, int height, int fps)
{
    if (pipe_)
    {
        return true; // already recording
    }

    if (width <= 0 || height <= 0)
    {
        status_ = "GIF: invalid frame size";
        return false;
    }

    if (!ffmpegAvailable())
    {
#ifdef _WIN32
        status_ = "GIF: ffmpeg not found (add ffmpeg to PATH)";
#else
        status_ = "GIF: ffmpeg not found (brew install ffmpeg)";
#endif
        return false;
    }

    w_ = width;
    h_ = height;
    fps_ = fps > 0 ? fps : 15;
    outputPath_ = defaultOutputPath();

    // Read rawvideo BGRA from stdin; produce an optimized GIF using a single-pass
    // palettegen/paletteuse. scale=iw/2 halves the (retina) pixel dimensions.
    char cmd[1024];
    std::snprintf(cmd, sizeof(cmd),
        "ffmpeg -y -loglevel error -f rawvideo -pix_fmt bgra -s %dx%d -r %d -i - "
        "-vf \"fps=%d,scale=iw/2:-1:flags=lanczos,split[s0][s1];"
        "[s0]palettegen[p];[s1][p]paletteuse\" \"%s\"",
        w_, h_, fps_, fps_, outputPath_.c_str());

    pipe_ = GIF_POPEN(cmd, GIF_PIPE_MODE);
    if (!pipe_)
    {
        status_ = "GIF: failed to launch ffmpeg";
        return false;
    }

    status_ = "\xE2\x97\x8F recording...";
    return true;
}

void GifRecorder::addFrame(const uint8_t* bgra, int width, int height)
{
    if (!pipe_ || !bgra)
    {
        return;
    }
    if (width != w_ || height != h_)
    {
        return; // window resized mid-recording; skip off-size frames
    }
    const size_t bytes = static_cast<size_t>(w_) * static_cast<size_t>(h_) * 4u;
    std::fwrite(bgra, 1, bytes, pipe_);
}

bool GifRecorder::stop()
{
    if (!pipe_)
    {
        return false;
    }

    const int rc = GIF_PCLOSE(pipe_);
    pipe_ = nullptr;

    if (rc != 0)
    {
        status_ = "GIF: ffmpeg failed to encode";
        return false;
    }

    return gif_platform::copyToClipboard(outputPath_, status_);
}
