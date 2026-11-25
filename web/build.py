"""Build the single-file browser app (web/abcheck.html) from the template and the engine.

The artifact viewer only runs inline scripts, so the JavaScript engine and the example
sequences are inlined into one HTML file. tests/test_web_app.py checks the committed file
is current.

    python3 web/build.py          # write web/abcheck.html
    python3 web/build.py --check  # exit 1 if web/abcheck.html is stale
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

WEB = Path(__file__).resolve().parent
ROOT = WEB.parent
TEMPLATE = WEB / "app.html"
ENGINE = WEB / "abcheck-engine.js"
APPLIB = WEB / "app-lib.js"
OUTPUT = WEB / "abcheck.html"
CORPUS = ROOT / "data" / "synthetic" / "sequences.json"

#: Example input shown on first load: synthetic sequences that cover the three triage
#: outcomes (no rule matched / needs review / cannot evaluate), both chains, and a spread
#: of pI for the property plot.
EXAMPLE_IDS = ["SYN-VH-001", "SYN-VH-003", "SYN-VH-007", "SYN-VH-005", "SYN-VH-013", "SYN-VH-902",
               "SYN-VL-001", "SYN-VL-004", "SYN-VL-005"]


def example_fasta() -> str:
    records = {r["id"]: r["sequence"] for r in json.loads(CORPUS.read_text(encoding="utf-8"))["sequences"]}
    return "".join(f">{i}\n{records[i]}\n" for i in EXAMPLE_IDS)


def build() -> str:
    template = TEMPLATE.read_text(encoding="utf-8")
    engine = ENGINE.read_text(encoding="utf-8")
    applib = APPLIB.read_text(encoding="utf-8")
    for marker in ("/*__ENGINE__*/", "/*__APPLIB__*/", '/*__EXAMPLE__*/""'):
        if template.count(marker) != 1:
            raise SystemExit(f"web/app.html must contain {marker} exactly once")
    for name, code in (("engine", engine), ("app-lib", applib)):
        if "</script" in code.lower():
            raise SystemExit(f"{name} contains '</script' and cannot be inlined")
    return (template
            .replace("/*__ENGINE__*/", engine)
            .replace("/*__APPLIB__*/", applib)
            .replace('/*__EXAMPLE__*/""', json.dumps(example_fasta(), ensure_ascii=False)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    html = build()
    if args.check:
        current = OUTPUT.read_text(encoding="utf-8") if OUTPUT.exists() else ""
        if current != html:
            print("web/abcheck.html is stale: run python3 web/build.py", file=sys.stderr)
            return 1
        print("web/abcheck.html is current")
        return 0
    OUTPUT.write_text(html, encoding="utf-8")
    print(f"wrote {OUTPUT.relative_to(ROOT)} ({len(html) // 1024} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
