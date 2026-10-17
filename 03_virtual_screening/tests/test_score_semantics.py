"""Score semantics: engine columns keep their real units and provenance."""

from __future__ import annotations

import os

import pytest

from zvss.ingest.poses import SCORE_SEMANTICS, normalise_scores, parse_gnina_log, read_gnina_log

REAL_RUN = r"F:\MemPilot\runs\public20_engineering\entries\P000002_S01_C01_seed20261015\stdout.log"


def test_cnn_affinity_documented_as_pk() -> None:
    entry = SCORE_SEMANTICS["CNNaffinity"]
    text = " ".join(str(value) for value in entry.values()).lower()
    assert "pk" in text or "not kcal" in text


def test_normalise_scores_keeps_engine_column_names() -> None:
    canonical, sources = normalise_scores(
        {
            "minimizedAffinity": "-9.0",
            "CNNscore": "0.5",
            "CNNaffinity": "7.5",
            "CNN_VS": "0.01",
            "not_a_gnina_column": "1",
        }
    )
    assert canonical["minimizedAffinity"] == pytest.approx(-9.0)
    assert canonical["CNNscore"] == pytest.approx(0.5)
    assert canonical["CNNaffinity"] == pytest.approx(7.5)
    assert "not_a_gnina_column" not in canonical
    assert sources["minimizedAffinity"].lower() == "minimizedaffinity"


def test_missing_values_become_none_not_zero() -> None:
    canonical, _ = normalise_scores({"minimizedAffinity": "", "CNNscore": "  "})
    assert canonical["minimizedAffinity"] is None
    assert canonical["CNNscore"] is None


@pytest.mark.skipif(not os.path.exists(REAL_RUN), reason="real GNINA output not available")
def test_real_gnina_log_parses_modes() -> None:
    log = read_gnina_log(REAL_RUN)
    assert log.modes, "expected at least one scored mode in a real GNINA log"
    first = log.modes[0]
    payload = first.to_dict() if hasattr(first, "to_dict") else vars(first)
    assert payload, "mode must serialise"
    assert log.argv and log.argv[0].lower().endswith("gnina")


def test_log_without_scoring_table_reports_no_modes() -> None:
    log = parse_gnina_log("# Commandline: gnina -r r.pdb\nno table here\n")
    assert log.modes == []
    assert log.warnings or log.argv
