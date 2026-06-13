"""CLI entry point for ``python -m grteclyn_wrapper`` and the ``grteclyn-wrapper`` console script."""

from __future__ import annotations

from .cli import build_parser, main

__all__ = ["build_parser", "main"]

if __name__ == "__main__":
    raise SystemExit(main())
