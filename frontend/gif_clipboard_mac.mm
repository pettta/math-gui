#include "gif_recorder.h"

#import <Foundation/Foundation.h>
#import <AppKit/AppKit.h>

namespace {
// UTI for GIF data on the pasteboard. Avoids depending on UTType SDK availability.
NSString* const kGifPasteboardType = @"com.compuserve.gif";
}

bool gif_platform::copyToClipboard(const std::string& path, std::string& status)
{
    NSString* nsPath = [NSString stringWithUTF8String:path.c_str()];
    NSData* gifData = [NSData dataWithContentsOfFile:nsPath];
    if (!gifData)
    {
        status = "GIF: encoded file missing";
        return false;
    }

    NSPasteboard* pb = [NSPasteboard generalPasteboard];
    [pb clearContents];

    // 1) File reference (pastes into Finder/Slack/Messages as an attachment).
    NSURL* fileURL = [NSURL fileURLWithPath:nsPath];
    BOOL wroteURL = [pb writeObjects:@[ fileURL ]];

    // 2) Raw GIF image data (some apps paste it inline as an animated image).
    [pb addTypes:@[ kGifPasteboardType ] owner:nil];
    BOOL wroteData = [pb setData:gifData forType:kGifPasteboardType];

    if (wroteURL || wroteData)
    {
        status = std::string("GIF copied to clipboard (") + path + ")";
        return true;
    }

    status = "GIF: clipboard copy failed";
    return false;
}
