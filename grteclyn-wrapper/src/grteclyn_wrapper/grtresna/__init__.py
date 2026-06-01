"""GRTresna elliptic solver bridge: Chombo HDF5 to GRTeclyn gridinit."""

__all__ = [
    "GRTresnaConfig",
    "convert_chombo_to_gridinit",
    "read_chombo_domain",
    "solve",
]


def __getattr__(name: str):
    if name in {"convert_chombo_to_gridinit", "read_chombo_domain"}:
        from .io import convert_chombo_to_gridinit, read_chombo_domain

        _map = {
            "convert_chombo_to_gridinit": convert_chombo_to_gridinit,
            "read_chombo_domain": read_chombo_domain,
        }
        return _map[name]
    if name in {"GRTresnaConfig", "solve"}:
        from . import solver

        return getattr(solver, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
