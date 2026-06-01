"""GRTresna elliptic solver bridge: Chombo HDF5 to GRTeclyn gridinit."""

__all__ = ["GRTresnaConfig", "convert_chombo_to_gridinit", "solve"]


def __getattr__(name: str):
    if name == "convert_chombo_to_gridinit":
        from .io import convert_chombo_to_gridinit

        return convert_chombo_to_gridinit
    if name in {"GRTresnaConfig", "solve"}:
        from . import solver

        return getattr(solver, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
