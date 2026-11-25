#!/usr/bin/env bash
# 사용 설명서 PDF 를 다시 만든다.
#
# 필요: python3, node, Playwright(Chromium)
#   Playwright 는 설명서를 만들 때만 쓴다. abcheck 자체는 표준 라이브러리만 쓴다.
#
#   bash docs/manual/build.sh
#
# 순서
#   1) build_manual.py — 현재 코드의 실제 출력으로 manual.html 생성
#   2) render.js       — sample_output 리포트 스크린샷 + PDF 인쇄
#   3) --stamp         — PDF 를 만든 원본 HTML 의 해시 기록 (테스트가 낡았는지 확인)
set -euo pipefail
cd "$(dirname "$0")/../.."

mkdir -p docs/manual/img
python3 docs/manual/build_manual.py
node docs/manual/render.js
python3 docs/manual/build_manual.py --stamp
