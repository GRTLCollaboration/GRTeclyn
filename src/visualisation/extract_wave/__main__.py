"""
Wave-extraction entry point.

This repository now uses a Psi4-only workflow (extract and analyze Weyl scalar Ψ4)
to deduce GW frequency content without reconstructing strain h(t).

`python -m src.visualisation.extract_wave ...` is kept as a convenience alias for:
`python -m src.visualisation.extract_wave.plot_psi4 ...`
"""

from .plot_psi4 import main


if __name__ == "__main__":
    main()