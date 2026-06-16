"""QD smoke: 4D evolving geodesic metrics must drive MAP-Elites placement."""

from __future__ import annotations

import json

import pytest

from grteclyn_wrapper.core.evaluation import Evaluation
from grteclyn_wrapper.search import qd_search as qd_module
from grteclyn_wrapper.search.ftl_retention import FTL_CHAMPIONS_FILE
from grteclyn_wrapper.search.optimize import SearchDimension
from grteclyn_wrapper.search.qd_search import QDArchive, run_qd_search


def _base_components(**extra: float) -> dict[str, float]:
    return {
        "constraint_growth": 0.0,
        "operational_ftl": 0.2,
        "ftl_persistence": 0.0,
        "ftl_geo_evolving": 0.0,
        "operational_ftl_geodesic": 0.0,
        "ftl_geo_timeavg": 0.0,
        **extra,
    }


def _frozen_artifact_metrics() -> dict:
    return {
        "ftl_timeseries": {
            "n_frames": 10,
            "f_geo_peak": 0.05,
            "ftl_lifetime_fraction": 0.9,
        },
        "evolving_geodesic": {
            "f_geo": 0.0,
            "h_quality_ok": True,
            "n_rays": 5,
            "n_reached": 5,
        },
    }


def _four_d_shortcut_metrics(f_geo: float = 0.12) -> dict:
    return {
        "ftl_timeseries": {
            "n_frames": 10,
            "f_geo_peak": 0.001,
            "ftl_lifetime_fraction": 0.1,
        },
        "evolving_geodesic": {
            "f_geo": f_geo,
            "h_quality_ok": True,
            "n_rays": 3,
            "n_reached": 3,
        },
    }


@pytest.fixture
def qd_4d_smoke_dir(tmp_path, monkeypatch):
    """Run a 2-eval QD batch with mocked GPU results contrasting frozen vs 4D geo."""

    def fake_evaluate_overrides(overrides, *, out_dir, name, **_kwargs):
        idx = int(name.rsplit("_", 1)[1])
        episode = out_dir / name
        episode.mkdir(parents=True, exist_ok=True)
        if idx == 1:
            return Evaluation(
                score=-10.0,
                components=_base_components(),
                notes=[],
                episode_path=str(episode),
                exit_code=0,
                preflight_rejected=False,
                reason=None,
                metrics=_frozen_artifact_metrics(),
            )
        return Evaluation(
            score=50.0,
            components=_base_components(ftl_geo_evolving=0.12),
            notes=[],
            episode_path=str(episode),
            exit_code=0,
            preflight_rejected=False,
            reason=None,
            metrics=_four_d_shortcut_metrics(),
        )

    monkeypatch.setattr(qd_module, "evaluate_overrides", fake_evaluate_overrides)

    qd_dir = tmp_path / "qd_4d_smoke"
    run_qd_search(
        runs_dir=tmp_path,
        name="qd_4d_smoke",
        iterations=0,
        init_random=2,
        batch_size=2,
        bins=8,
        seed=7,
        descriptor_mode="ftl_lifetime",
        search_space=[SearchDimension("toy_param", 0.0, 1.0)],
        consume_plotfiles=False,
        check_params=False,
        ftl_retention_enabled=True,
    )
    return qd_dir


def test_qd_map_elites_bins_on_4d_not_frozen_peak(qd_4d_smoke_dir) -> None:
    rows = [
        json.loads(line)
        for line in (qd_4d_smoke_dir / "trajectory.jsonl").read_text(encoding="utf-8").splitlines()
    ]
    assert len(rows) == 2

    artifact = next(row for row in rows if row["eval"] == 1)
    genuine = next(row for row in rows if row["eval"] == 2)

    art_details = artifact["descriptor_details"]
    gen_details = genuine["descriptor_details"]

    assert art_details["x"] == 0.0
    assert art_details["y"] == 0.0
    assert art_details["f_geo_evol"] == 0.0
    assert art_details["ftl_geo_timeavg"] == 0.0

    assert gen_details["x"] > 0.0
    assert gen_details["y"] == 1.0
    assert gen_details["f_geo_evol"] == pytest.approx(0.12)
    assert gen_details["ftl_geo_evolving"] == pytest.approx(0.12)

    assert artifact["cell"] != genuine["cell"]
    assert genuine["cell"][1] > artifact["cell"][1]


def test_qd_archive_elites_carry_4d_descriptor_details(qd_4d_smoke_dir) -> None:
    archive = QDArchive.from_dict(
        json.loads((qd_4d_smoke_dir / "archive.json").read_text(encoding="utf-8"))
    )
    assert len(archive.cells) == 2

    by_score = sorted(archive.cells.values(), key=lambda elite: elite.score)
    low, high = by_score

    assert low.descriptor_details["x"] == 0.0
    assert low.descriptor_details["f_geo_evol"] == 0.0
    assert high.descriptor_details["x"] > 0.0
    assert high.descriptor_details["y"] == 1.0
    assert high.descriptor_details["f_geo_evol"] == pytest.approx(0.12)


def test_qd_ftl_retention_crowns_4d_geo_not_frozen_peak(qd_4d_smoke_dir) -> None:
    champs = json.loads(
        (qd_4d_smoke_dir / FTL_CHAMPIONS_FILE).read_text(encoding="utf-8")
    )
    assert "f_geo_evol" in champs
    assert champs["f_geo_evol"]["eval"] == 2
    assert champs["f_geo_evol"]["value"] == pytest.approx(0.12)
    assert "f_geo_peak" not in champs
    assert "ftl_geo_timeavg" not in champs
