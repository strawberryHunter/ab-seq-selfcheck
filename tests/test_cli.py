"""CLI tests: actually run via subprocess and check exit codes and outputs."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from abcheck.cli import load_input, main, read_fasta

from .conftest import CORPUS_PATH, by_id

REPO_ROOT = CORPUS_PATH.resolve().parents[2]


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "abcheck.cli", *args],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )


# --------------------------------------------------------------------------
# Acceptance criterion: corpus -> HTML
# --------------------------------------------------------------------------

def test_corpus_to_html_creates_report(tmp_path):
    result = _run(str(CORPUS_PATH), "--html", "--out-dir", str(tmp_path))

    assert result.returncode == 0, result.stderr
    report = tmp_path / "report.html"
    assert report.exists()

    out = report.read_text(encoding="utf-8")
    assert out.startswith("<!doctype html>")
    assert "Arial" in out
    # Every sequence in the corpus appears in the report
    corpus = json.loads(CORPUS_PATH.read_text(encoding="utf-8"))
    for rec in corpus["sequences"]:
        assert rec["id"] in out
    assert "HTML 리포트를 썼다" in result.stderr


def test_html_mode_does_not_dump_table_to_stdout(tmp_path):
    result = _run(str(CORPUS_PATH), "--html", "--out-dir", str(tmp_path))
    assert result.stdout.strip() == ""


def test_not_confident_sequences_are_listed_on_stderr(tmp_path):
    result = _run(str(CORPUS_PATH), "--html", "--out-dir", str(tmp_path))
    assert "주석을 신뢰할 수 없는 서열 5건" in result.stderr
    for seq_id in ["SYN-VH-901", "SYN-VH-902", "SYN-VH-903", "SYN-VH-013", "SYN-VH-014"]:
        assert seq_id in result.stderr
    assert "ANARCI/IMGT 재확인 필요" in result.stderr


# --------------------------------------------------------------------------
# Pasting a single sequence
# --------------------------------------------------------------------------

def test_single_sequence_prints_table(records):
    seq = by_id(records, "SYN-VH-001")["sequence"]
    result = _run(seq)

    assert result.returncode == 0, result.stderr
    assert "CDR-H1" in result.stdout
    assert "CDR-H3" in result.stdout
    assert "중쇄 (VH)" in result.stdout
    assert "numbering 엔진이" in result.stdout


def test_single_sequence_json_mode(records):
    seq = by_id(records, "SYN-VL-001")["sequence"]
    result = _run(seq, "--json")

    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["count"] == 1
    entry = payload["sequences"][0]
    assert entry["annotation"]["chain"] == "L"
    assert set(entry["annotation"]["cdrs"]) == {"CDR-L1", "CDR-L2", "CDR-L3"}
    assert payload["limitations"]


def test_lowercase_input_is_accepted(records):
    seq = by_id(records, "SYN-VH-001")["sequence"].lower()
    result = _run(seq, "--json")
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["sequences"][0]["annotation"]["chain"] == "H"


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------

def test_fasta_file_is_read(tmp_path, records):
    recs = [by_id(records, i) for i in ["SYN-VH-001", "SYN-VL-001"]]
    path = tmp_path / "input.fasta"
    path.write_text(
        "".join(f">{r['id']} 설명\n{r['sequence']}\n" for r in recs), encoding="utf-8"
    )

    result = _run(str(path), "--json")
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert [e["id"] for e in payload["sequences"]] == ["SYN-VH-001", "SYN-VL-001"]


def test_read_fasta_handles_wrapped_lines():
    text = ">a desc\nACDE\nFGHI\n>b\nKLMN\n"
    assert read_fasta(text) == [("a", "ACDEFGHI"), ("b", "KLMN")]


def test_read_fasta_ignores_blank_lines():
    assert read_fasta("\n>a\n\nACDE\n\n") == [("a", "ACDE")]


def test_pasted_fasta_is_detected(records):
    rec = by_id(records, "SYN-VH-001")
    result = _run(f">{rec['id']}\n{rec['sequence']}", "--json")
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["sequences"][0]["id"] == "SYN-VH-001"


# --------------------------------------------------------------------------
# Error handling
# --------------------------------------------------------------------------

def test_invalid_sequence_exits_with_one():
    result = _run("ACDEFGHIK")
    assert result.returncode == 1
    assert "입력 검증에서 걸러진 서열" in result.stderr
    assert "너무 짧다" in result.stderr


def test_dna_input_is_rejected_with_korean_message():
    result = _run("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG" * 3)
    assert result.returncode == 1
    assert "DNA/RNA" in result.stderr


def test_missing_file_exits_with_two():
    result = _run("no_such_file.fasta")
    assert result.returncode == 2
    assert "파일을 찾을 수 없다" in result.stderr


def test_malformed_json_exits_with_two(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text('{"sequences": [{"no_sequence_key": 1}]}', encoding="utf-8")
    result = _run(str(path))
    assert result.returncode == 2
    assert "'sequence' 키가 없다" in result.stderr


def test_valid_input_with_unconfident_sequence_still_exits_zero(records):
    """Low annotation confidence still exits successfully if the input itself is valid."""
    seq = by_id(records, "SYN-VH-901")["sequence"]
    result = _run(seq, "--json")
    assert result.returncode == 0
    payload = json.loads(result.stdout)
    assert payload["sequences"][0]["annotation"]["confident"] is False


# --------------------------------------------------------------------------
# Input type detection
# --------------------------------------------------------------------------

def test_load_input_detects_corpus_json():
    entries, source = load_input(str(CORPUS_PATH))
    assert len(entries) == 24
    assert "코퍼스 JSON" in source


def test_load_input_treats_bare_sequence_as_sequence(records):
    seq = by_id(records, "SYN-VH-001")["sequence"]
    entries, source = load_input(seq)
    assert entries == [("입력 서열", seq)]
    assert source == "붙여넣은 서열"


def test_load_input_accepts_plain_list_json(tmp_path):
    path = tmp_path / "list.json"
    path.write_text(json.dumps([{"id": "x", "sequence": "ACDE"}]), encoding="utf-8")
    entries, _ = load_input(str(path))
    assert entries == [("x", "ACDE")]


def test_load_input_accepts_string_list_json(tmp_path):
    path = tmp_path / "strings.json"
    path.write_text(json.dumps(["ACDE", "FGHI"]), encoding="utf-8")
    entries, _ = load_input(str(path))
    assert entries == [("seq1", "ACDE"), ("seq2", "FGHI")]


# --------------------------------------------------------------------------
# Options
# --------------------------------------------------------------------------

def test_custom_ph_changes_net_charge(records):
    seq = by_id(records, "SYN-VH-001")["sequence"]

    low = json.loads(_run(seq, "--json", "--ph", "4.0").stdout)
    high = json.loads(_run(seq, "--json", "--ph", "10.0").stdout)

    low_charge = low["sequences"][0]["properties"]["도메인 전체"]["net_charge"]
    high_charge = high["sequences"][0]["properties"]["도메인 전체"]["net_charge"]
    assert low_charge > high_charge
    assert low["sequences"][0]["properties"]["도메인 전체"]["pH"] == 4.0


def test_custom_out_name(tmp_path):
    result = _run(str(CORPUS_PATH), "--html", "--out-dir", str(tmp_path),
                  "--out-name", "custom.html")
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "custom.html").exists()


def test_out_dir_is_created_when_missing(tmp_path, records):
    target = tmp_path / "deep" / "nested"
    seq = by_id(records, "SYN-VH-001")["sequence"]
    result = _run(seq, "--html", "--out-dir", str(target))
    assert result.returncode == 0, result.stderr
    assert (target / "report.html").exists()


def test_help_lists_three_examples():
    result = _run("--help")
    assert result.returncode == 0
    assert result.stdout.count("python -m abcheck.cli") >= 4   # prog line + 3 examples


def test_main_is_callable_in_process(tmp_path):
    code = main([str(CORPUS_PATH), "--html", "--out-dir", str(tmp_path)])
    assert code == 0
    assert (tmp_path / "report.html").exists()
