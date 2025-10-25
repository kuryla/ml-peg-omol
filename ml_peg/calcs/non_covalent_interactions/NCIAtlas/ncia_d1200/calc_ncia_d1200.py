from __future__ import annotations

# Wrapper to keep legacy location discoverable under new category
from ml_peg.calcs.NCIAtlas.ncia_d1200.calc_ncia_d1200 import build_project as _build, test_ncia_d1200 as _test


def build_project(repro: bool = False) -> None:  # noqa: D401
    return _build(repro)


def test_ncia_d1200():  # noqa: D401
    return _test()

