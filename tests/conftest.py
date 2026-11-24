"""테스트 공용 fixture — 합성 코퍼스를 읽어 온다."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

CORPUS_PATH = Path(__file__).resolve().parents[1] / "data" / "synthetic" / "sequences.json"


@pytest.fixture(scope="session")
def corpus() -> dict:
    assert CORPUS_PATH.exists(), (
        f"{CORPUS_PATH} 가 없다. 먼저 `python3 scripts/make_synthetic_sequences.py` 를 실행하라."
    )
    return json.loads(CORPUS_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="session")
def records(corpus: dict) -> list[dict]:
    return corpus["sequences"]


def by_id(records: list[dict], seq_id: str) -> dict:
    for rec in records:
        if rec["id"] == seq_id:
            return rec
    raise KeyError(f"코퍼스에 {seq_id} 가 없다")
