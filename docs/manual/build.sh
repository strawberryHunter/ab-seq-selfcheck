#!/usr/bin/env bash
# Rebuild the user manual PDF.
#
# Requires: python3, node, Playwright (Chromium)
#   Playwright is used only to build the manual. abcheck itself uses only the standard library.
#
#   bash docs/manual/build.sh
#
# Steps
#   1) build_manual.py — generate manual.html from the current code's actual output
#   2) render.js       — screenshot the sample_output report + print the PDF
#   3) --stamp         — record the hash of the source HTML the PDF was built from (tests check for staleness)
set -euo pipefail
cd "$(dirname "$0")/../.."

mkdir -p docs/manual/img
python3 docs/manual/build_manual.py
node docs/manual/render.js
python3 docs/manual/build_manual.py --stamp
