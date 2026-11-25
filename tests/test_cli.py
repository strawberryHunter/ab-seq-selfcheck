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
        encoding="utf-8",
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
    assert "입력 검증에서 걸러진 서열" in result.stdout
    assert "너무 짧다" in result.stdout


def test_all_rejected_input_shows_only_the_fix():
    """When nothing passed validation, show the error and how to proceed -- not a full report.

    A usability pass found a 3-line fix printed under ~45 lines of notices and empty tables.
    """
    result = _run("ACDEFGHIK")
    lines = [l for l in result.stdout.splitlines() if l.strip()]
    assert len(lines) <= 6, result.stdout
    assert "다시 실행" in result.stdout
    assert "[이 도구의 한계]" not in result.stdout


def test_dna_input_is_rejected_with_korean_message():
    result = _run("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG" * 3)
    assert result.returncode == 1
    assert "DNA/RNA" in result.stdout


def test_trailing_stop_codon_gets_a_specific_hint(records):
    result = _run(records[0]["sequence"] + "*")
    assert result.returncode == 1
    assert "정지 코돈" in result.stdout


def test_mixed_batch_still_lists_rechecks_and_summary(tmp_path):
    """One bad sequence must not hide the recheck list for the others."""
    corpus = json.loads(CORPUS_PATH.read_text(encoding="utf-8"))
    by_id = {r["id"]: r["sequence"] for r in corpus["sequences"]}
    fasta = tmp_path / "mixed.fasta"
    fasta.write_text(
        f">good\n{by_id['SYN-VH-001']}\n>shaky\n{by_id['SYN-VH-903']}\n>bad\nACDEFGHIK\n",
        encoding="utf-8",
    )
    result = _run(str(fasta))
    assert result.returncode == 1
    assert "입력 검증에서 걸러진 서열 1건" in result.stderr
    assert "ANARCI/IMGT 재확인 필요): shaky" in result.stderr
    assert "3건 중 2건 분석, 1건은 입력 오류로 제외" in result.stderr


def test_argument_errors_are_korean_with_an_example():
    result = _run()
    assert result.returncode == 2
    assert "검사할 서열이나 파일을 넣지 않았다" in result.stderr
    assert "사용 예" in result.stderr
    assert "required" not in result.stderr


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
    assert result.stdout.count("python3 -m abcheck.cli") >= 4   # prog line + 3 examples


def test_main_is_callable_in_process(tmp_path):
    code = main([str(CORPUS_PATH), "--html", "--out-dir", str(tmp_path)])
    assert code == 0
    assert (tmp_path / "report.html").exists()


# --------------------------------------------------------------------------
# Files as researchers' tools save them (Windows Notepad / Excel / Korean Windows)
# --------------------------------------------------------------------------

import os  # noqa: E402

import pytest  # noqa: E402


def _fasta_text(records, n=3):
    return "".join(f">{r['id']} 설명\n{r['sequence']}\n" for r in records[:n])


@pytest.mark.parametrize("encoding", ["utf-8", "utf-8-sig", "utf-16", "cp949"])
@pytest.mark.parametrize("newline", ["\n", "\r\n"])
def test_fasta_is_read_in_every_common_encoding(tmp_path, records, encoding, newline):
    """A UTF-8 BOM once made the first record vanish silently; UTF-16/cp949 crashed."""
    path = tmp_path / "후보 목록.fasta"
    path.write_bytes(_fasta_text(records).replace("\n", newline).encode(encoding))
    entries, _ = load_input(str(path))
    assert [e[0] for e in entries] == [r["id"] for r in records[:3]]
    assert [e[1] for e in entries] == [r["sequence"] for r in records[:3]]


@pytest.mark.parametrize("header", ["id,sequence", "이름,서열", "Name,Sequence", None])
def test_csv_saved_from_excel(tmp_path, records, header):
    lines = ([header] if header else []) + [f"{r['id']},{r['sequence']}" for r in records[:3]]
    path = tmp_path / "후보.csv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8-sig")   # Excel "CSV UTF-8"
    entries, source = load_input(str(path))
    assert "표 파일" in source
    assert entries == [(r["id"], r["sequence"]) for r in records[:3]]


def test_tsv_with_extra_columns(tmp_path, records):
    path = tmp_path / "t.tsv"
    path.write_text("번호\t메모\t서열\n" + f"7\t좋음\t{records[0]['sequence']}\n", encoding="utf-8")
    assert load_input(str(path))[0] == [("7", records[0]["sequence"])]


def test_xlsx_gets_save_as_csv_guidance(tmp_path):
    path = tmp_path / "후보.xlsx"
    path.write_bytes(b"PK\x03\x04")
    result = _run(str(path))
    assert result.returncode == 2
    assert "CSV UTF-8" in result.stderr


def test_folder_instead_of_file(tmp_path):
    result = _run(str(tmp_path))
    assert result.returncode == 2
    assert "폴더가 아니라 파일" in result.stderr


def test_broken_json_message_is_korean_with_location(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text('{"sequences": [ {"id": "a" ', encoding="utf-8")
    result = _run(str(path))
    assert result.returncode == 2
    assert "JSON 형식이 깨져 있다" in result.stderr
    assert "Traceback" not in result.stderr


def test_unwritable_out_dir_is_reported_not_crashed(tmp_path):
    blocker = tmp_path / "file"
    blocker.write_text("x", encoding="utf-8")
    result = _run(str(CORPUS_PATH), "--html", "--out-dir", str(blocker / "sub"))
    assert result.returncode == 2
    assert "저장하지 못했다" in result.stderr
    assert "Traceback" not in result.stderr


@pytest.mark.parametrize("mode", [[], ["--json"]])
def test_redirected_output_on_korean_windows_does_not_crash(mode):
    """cp949 streams cannot encode an em dash; `--json > result.json` used to die."""
    env = dict(os.environ, PYTHONIOENCODING="cp949")
    result = subprocess.run(
        [sys.executable, "-m", "abcheck.cli", str(CORPUS_PATH), *mode],
        capture_output=True, env=env, cwd=Path(__file__).resolve().parents[1],
    )
    assert result.returncode == 0, result.stderr.decode("utf-8", "replace")
    text = result.stdout.decode("utf-8")              # written as UTF-8 when redirected
    if mode:
        assert json.loads(text)["count"] == 24


def test_old_python_gets_a_clear_message():
    """The version guard must stay parseable by old interpreters and say what to do."""
    import ast

    source = (Path(__file__).resolve().parents[1] / "abcheck" / "__init__.py").read_text(encoding="utf-8")
    ast.parse(source, feature_version=(3, 4))
    assert "sys.version_info < (3, 9)" in source
    assert "Python 3.9 이상이 필요하다" in source


# --------------------------------------------------------------------------
# Sequences copied from IMGT/V-QUEST or an alignment (gap characters)
# --------------------------------------------------------------------------

@pytest.mark.parametrize("gapped", [
    lambda s: s[:30] + "...." + s[30:70] + "---" + s[70:],    # IMGT-style gaps
    lambda s: s[:50] + "." + s[50:],                          # a single gap
])
def test_gapped_sequence_is_analyzed_not_mistaken_for_a_file(records, gapped):
    """A pasted sequence containing '.' used to be reported as a missing file."""
    seq = records[0]["sequence"]
    result = _run(gapped(seq), "--json")
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["sequences"][0]["sequence"] == seq


@pytest.mark.parametrize("name", ["후보.fast", "candidates.fast", "nofile.fasta", "x.csv"])
def test_filename_typos_still_say_file_not_found(name):
    result = _run(name)
    assert result.returncode == 2
    assert "파일을 찾을 수 없다" in result.stderr


def test_mostly_punctuation_is_still_rejected():
    result = _run("....----..AC")
    assert result.returncode == 1


def test_broken_csv_is_a_korean_message_not_a_traceback(tmp_path):
    """csv.Error is not a ValueError; it used to escape as a traceback."""
    path = tmp_path / "bad.csv"
    path.write_bytes(b"id,sequence\na,QVQL\rEVQL\nb,QQ\n")
    result = _run(str(path))
    assert result.returncode == 2
    assert "Traceback" not in result.stderr
    assert "번째 줄 근처에서 형식이 깨져 있다" in result.stderr


def test_gravy_is_exact_on_every_python_version():
    """GRAVY used to depend on sum()'s float algorithm, which changed in Python 3.12."""
    from abcheck.properties import KYTE_DOOLITTLE, gravy

    seq = "".join(KYTE_DOOLITTLE) * 7 + "AILV"
    tenths = sum(round(KYTE_DOOLITTLE[r] * 10) for r in seq)
    assert gravy(seq) == tenths / (10 * len(seq))
