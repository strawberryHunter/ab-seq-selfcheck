"""The browser engine (web/abcheck-engine.js) must give exactly the Python output.

Runs scripts/js_parity.py's checks through `node`. Skipped when node is not installed,
because Node is only a test-time tool here, not a dependency of abcheck.
"""

from __future__ import annotations

import importlib.util
import json
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

@pytest.fixture(scope="module")
def parity():
    """scripts/js_parity.py as a module; every test using it is skipped without node."""
    if shutil.which("node") is None:
        pytest.skip("node 가 PATH 에 없어 JS 엔진 비교를 건너뛴다 (Node.js 를 설치하면 실행된다)")
    spec = importlib.util.spec_from_file_location("js_parity", ROOT / "scripts" / "js_parity.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_analyze_matches_python_on_corpus_and_edge_cases(parity):
    n, diffs = parity.analyze_cases(parity.standard_cases())
    assert n > 400
    assert diffs == [], "\n".join(diffs[:20])


def test_edge_cases_reach_the_fallbacks_they_name(parity):
    """The rescued-motif cases really go through the J-tail fallback in Python."""
    from abcheck.report import build_many

    cases = dict(parity.edge_cases())
    reports = {r.seq_id: r for r in build_many(list(cases.items()))}
    for seq_id, chain in (("vh-j-tail", "H"), ("vl-j-tail-lambda", "L"), ("vl-j-tail-kappa", "L")):
        report = reports[seq_id]
        assert report.chain == chain
        assert any("J 유전자 꼬리" in w for w in report.annotation.warnings)
    assert not reports["stop"].ok
    assert "'*'" in reports["stop"].errors[0]
    assert "'\\ufeff'" in reports["bom"].errors[0]


def test_ph_variants_change_charge(parity):
    """pH 5.0 / 9.0 cases are not accidentally identical to pH 7.4."""
    vh = dict(parity.corpus_entries())["SYN-VH-001"]
    answers = parity.run_node([
        {"op": "analyze", "entries": [["x", vh]], "ph": ph} for ph in (5.0, 7.4, 9.0)
    ])
    charges = [json.loads(a["text"])["sequences"][0]["properties"]["도메인 전체"]["net_charge"]
               for a in answers]
    assert charges[0] > charges[1] > charges[2]
    assert '"pH": 5.0' in answers[0]["text"] and '"pH": 9.0' in answers[2]["text"]


def test_number_helpers_match_python(parity):
    n, diffs = parity.helper_diffs()
    assert n > 100000
    assert diffs == [], "\n".join(diffs[:20])


def _js_read(parity, text, filename=None):
    return parity.run_node([{"op": "readInput", "text": text, "filename": filename}])[0]


def test_read_input_fasta_string(parity):
    text = ">VH-1 heavy chain\nQVQLVQSG\nAEVKK\n\n>VL-1\r\nDIQMTQ\r\n>\nEVQL\n"
    py = parity.py_read_input(text, "후보.fasta")
    assert py == _js_read(parity, text, "후보.fasta")
    assert py["entries"] == [["VH-1", "QVQLVQSGAEVKK"], ["VL-1", "DIQMTQ"], ["seq3", "EVQL"]]
    assert py["source"] == "FASTA (후보.fasta)"
    # pasted, without a file name
    assert parity.py_read_input(text) == _js_read(parity, text)
    assert _js_read(parity, text)["source"] == "붙여넣은 FASTA"


def test_read_input_csv_with_korean_headers(parity):
    text = "﻿후보,아미노산 서열,비고\r\n  B1 ,QVQLVQ,x\r\n,DIQMTQ,\r\n\r\n"
    py = parity.py_read_input(text, "표.csv")
    assert py == _js_read(parity, text, "표.csv")
    assert py["entries"] == [["B1", "QVQLVQ"], ["seq2", "DIQMTQ"]]
    tsv = "이름\t서열\nT1\tQVQ\n"
    assert parity.py_read_input(tsv, "t.tsv") == _js_read(parity, tsv, "t.tsv")
    no_header = "a,b,c\n1,QVQ,3\n"
    py_err = parity.py_read_input(no_header, "x.csv")
    assert not py_err["ok"] and "서열 열을 찾지 못했다" in py_err["error"]
    assert py_err == _js_read(parity, no_header, "x.csv")


def test_read_input_json_string(parity):
    text = json.dumps({"sequences": [{"id": "J1", "sequence": "QVQ"}, {"sequence": "DIQ"},
                                     {"id": 7, "sequence": "EVQ"}]})
    py = parity.py_read_input(text, "c.json")
    assert py == _js_read(parity, text, "c.json")
    assert py["entries"] == [["J1", "QVQ"], ["seq2", "DIQ"], ["7", "EVQ"]]
    broken = '{"sequences": [\n  {"id": "a" "sequence": "Q"}\n]}'
    py_err = parity.py_read_input(broken, "b.json")
    assert py_err == {"ok": False, "error": (
        "JSON 형식이 깨져 있다 (2번째 줄 14번째 글자 근처). "
        "쉼표·따옴표·괄호가 빠지지 않았는지 확인하라.")}
    assert py_err == _js_read(parity, broken, "b.json")


def test_read_input_raw_sequence(parity):
    text = "  qvqlvqsg aevkk\n"
    py = parity.py_read_input(text)
    assert py == _js_read(parity, text)
    assert py == {"ok": True, "entries": [["입력 서열", text]], "source": "붙여넣은 서열"}


def test_read_input_fixed_and_fuzzed_cases(parity):
    cases = parity.read_input_cases() + parity.fuzz_read_input_cases()
    diffs = parity.read_input_diffs(cases)
    assert diffs == [], "\n".join(diffs[:20])


def test_engine_is_a_plain_script():
    """One file, no imports or exports, so a <script> tag can load it as-is."""
    source = (ROOT / "web" / "abcheck-engine.js").read_text(encoding="utf-8")
    for line in source.splitlines():
        stripped = line.strip()
        assert not stripped.startswith(("import ", "export ")), line
    assert "require(" not in source
