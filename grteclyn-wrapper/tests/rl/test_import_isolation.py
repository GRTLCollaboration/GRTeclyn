from __future__ import annotations

from pathlib import Path


def test_core_modules_do_not_import_rl() -> None:
    root = Path("src/grteclyn_wrapper")
    forbidden_roots = ("core", "search", "metrics")
    for path in root.rglob("*.py"):
        rel = path.relative_to(root)
        if rel.parts[0] not in forbidden_roots:
            continue
        text = path.read_text(encoding="utf-8")
        assert "grteclyn_wrapper.rl" not in text, f"{path} imports rl package"
