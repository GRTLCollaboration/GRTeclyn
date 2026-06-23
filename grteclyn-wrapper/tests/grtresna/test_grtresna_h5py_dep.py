"""GRTresna Chombo→gridinit bridge requires h5py in the wrapper venv."""

from __future__ import annotations


def test_h5py_import_available() -> None:
    import h5py  # noqa: F401 — required by grtresna/io.py at solve time


def test_convert_chombo_imports_h5py_lazily() -> None:
    from grteclyn_wrapper.grtresna.io import convert_chombo_to_gridinit

    assert callable(convert_chombo_to_gridinit)
