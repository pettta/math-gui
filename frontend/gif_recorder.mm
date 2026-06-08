#include "gif_recorder.h"

#import <Foundation/Foundation.h>
#import <AppKit/AppKit.h>

#include <cstdlib>

namespace {
// UTI for GIF data on the pasteboard. Avoids depending on UTType SDK availability.
NSString* const kGifPasteboardType = @"com.compuserve.gif";
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

    if (std::system("command -v ffmpeg >/dev/null 2>&1") != 0)
    {
        status_ = "GIF: ffmpeg not found (brew install ffmpeg)";
        return false;
    }

    w_ = width;
    h_ = height;
    fps_ = fps > 0 ? fps : 15;
    outputPath_ = "/tmp/mathgui_recording.gif";

    // Read rawvideo BGRA from stdin; produce an optimized GIF using a single-pass
    // palettegen/paletteuse. scale=iw/2 halves the (retina) pixel dimensions.
    char cmd[1024];
    std::snprintf(cmd, sizeof(cmd),
        "ffmpeg -y -loglevel error -f rawvideo -pix_fmt bgra -s %dx%d -r %d -i - "
        "-vf \"fps=%d,scale=iw/2:-1:flags=lanczos,split[s0][s1];"
        "[s0]palettegen[p];[s1][p]paletteuse\" \"%s\"",
        w_, h_, fps_, fps_, outputPath_.c_str());

    pipe_ = popen(cmd, "w");
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

    const int rc = pclose(pipe_);
    pipe_ = nullptr;

    if (rc != 0)
    {
        status_ = "GIF: ffmpeg failed to encode";
        return false;
    }

    NSString* path = [NSString stringWithUTF8String:outputPath_.c_str()];
    NSData* gifData = [NSData dataWithContentsOfFile:path];
    if (!gifData)
    {
        status_ = "GIF: encoded file missing";
        return false;
    }

    NSPasteboard* pb = [NSPasteboard generalPasteboard];
    [pb clearContents];

    // 1) File reference (pastes into Finder/Slack/Messages as an attachment).
    NSURL* fileURL = [NSURL fileURLWithPath:path];
    BOOL wroteURL = [pb writeObjects:@[ fileURL ]];

    // 2) Raw GIF image data (some apps paste it inline as an animated image).
    [pb addTypes:@[ kGifPasteboardType ] owner:nil];
    BOOL wroteData = [pb setData:gifData forType:kGifPasteboardType];

    if (wroteURL || wroteData)
    {
        status_ = std::string("GIF copied to clipboard (") + outputPath_ + ")";
        return true;
    }

    status_ = "GIF: clipboard copy failed";
    return false;
}
