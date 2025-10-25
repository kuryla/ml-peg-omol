from ml_peg.calcs.NCIAtlas.ncia_r739x5.calc_ncia_r739x5 import build_project as _build, test_ncia_r739x5 as _test


def build_project(repro: bool = False) -> None:
    return _build(repro)


def test_ncia_r739x5():
    return _test()

