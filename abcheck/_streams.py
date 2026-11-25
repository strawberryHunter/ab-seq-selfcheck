"""Output-stream setup shared by the CLI and the scripts."""

from __future__ import annotations

import sys


def configure_streams() -> None:
    """Never crash on output encoding.

    On Windows, a redirected or piped stream uses the ANSI code page (cp1252, or cp949 on
    Korean Windows), which cannot encode Korean text or an em dash, so the run died with a
    traceback. Redirected output is written as UTF-8 (what JSON readers and editors expect);
    a terminal keeps its encoding and replaces anything it cannot show.
    """
    for stream in (sys.stdout, sys.stderr):
        if not hasattr(stream, "reconfigure"):
            continue
        try:
            if stream.isatty():
                stream.reconfigure(errors="replace")
            else:
                stream.reconfigure(encoding="utf-8")
        except (OSError, ValueError):
            pass
