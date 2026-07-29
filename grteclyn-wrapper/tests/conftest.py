from __future__ import annotations

import numpy as np
import pytest

from grteclyn_wrapper.core import scratch


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(42)


@pytest.fixture(autouse=True)
def isolated_plotfile_scratch(tmp_path_factory, monkeypatch):
    """Keep the suite out of the real ``/tmp/grteclyn_scratch``.

    Anything that writes params or builds a consumer command now resolves a
    node-local scratch directory, and without this every such test would leave
    a directory in the scratch root a live campaign is using.  Scratch stays
    *enabled* -- pointed somewhere disposable -- so tests still exercise the
    real path; a test that cares about a specific root overrides the variable
    itself.
    """
    root = tmp_path_factory.mktemp("grteclyn_scratch")
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(root))
    scratch._root_cache.clear()
    yield root
    scratch._root_cache.clear()
