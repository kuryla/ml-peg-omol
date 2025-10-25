from ml_peg.calcs.NCIAtlas.ncia_d442x10.calc_ncia_d442x10 import build_project as _build, test_ncia_d442x10 as _test


def build_project(repro: bool = False) -> None:
    return _build(repro)


def test_ncia_d442x10():
    return _test()

