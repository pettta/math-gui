#include "gif_recorder.h"

#ifdef _WIN32

#include <windows.h>
#include <shlobj.h> // DROPFILES
#include <vector>

bool gif_platform::copyToClipboard(const std::string& path, std::string& status)
{
    // Convert the UTF-8 path to a wide string for the Win32 clipboard APIs.
    const int wlen = MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, nullptr, 0);
    if (wlen <= 0)
    {
        status = "GIF: clipboard copy failed (bad path)";
        return false;
    }
    std::vector<wchar_t> wpath(static_cast<size_t>(wlen));
    MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, wpath.data(), wlen);

    // Build a CF_HDROP payload: a DROPFILES header followed by a
    // double-null-terminated list of wide file paths. This pastes the GIF as a
    // file attachment into Explorer/Slack/Discord/etc.
    const size_t pathBytes = static_cast<size_t>(wlen) * sizeof(wchar_t);
    const size_t totalBytes = sizeof(DROPFILES) + pathBytes + sizeof(wchar_t); // extra null terminates the list

    HGLOBAL hGlobal = GlobalAlloc(GHND, totalBytes);
    if (!hGlobal)
    {
        status = "GIF: clipboard copy failed (alloc)";
        return false;
    }

    auto* dropFiles = static_cast<DROPFILES*>(GlobalLock(hGlobal));
    if (!dropFiles)
    {
        GlobalFree(hGlobal);
        status = "GIF: clipboard copy failed (lock)";
        return false;
    }
    dropFiles->pFiles = sizeof(DROPFILES);
    dropFiles->fWide = TRUE;
    auto* dest = reinterpret_cast<wchar_t*>(reinterpret_cast<char*>(dropFiles) + sizeof(DROPFILES));
    memcpy(dest, wpath.data(), pathBytes);
    dest[wlen] = L'\0'; // second terminator (GHND already zeroed the buffer)
    GlobalUnlock(hGlobal);

    if (!OpenClipboard(nullptr))
    {
        GlobalFree(hGlobal);
        status = "GIF: clipboard copy failed (open)";
        return false;
    }
    EmptyClipboard();
    if (SetClipboardData(CF_HDROP, hGlobal))
    {
        // Ownership of hGlobal transfers to the clipboard on success.
        CloseClipboard();
        status = std::string("GIF copied to clipboard (") + path + ")";
        return true;
    }

    CloseClipboard();
    GlobalFree(hGlobal);
    status = "GIF: clipboard copy failed";
    return false;
}

#else // Linux / other Unix

#include <cstdlib>

bool gif_platform::copyToClipboard(const std::string& path, std::string& status)
{
    // No portable clipboard API on Linux; best-effort copy via xclip if present,
    // otherwise just report where the GIF was saved.
    if (std::system("command -v xclip >/dev/null 2>&1") == 0)
    {
        const std::string cmd =
            "xclip -selection clipboard -t image/gif -i \"" + path + "\" >/dev/null 2>&1 &";
        if (std::system(cmd.c_str()) == 0)
        {
            status = std::string("GIF copied to clipboard (") + path + ")";
            return true;
        }
    }

    status = std::string("GIF saved to ") + path;
    return true;
}

#endif
