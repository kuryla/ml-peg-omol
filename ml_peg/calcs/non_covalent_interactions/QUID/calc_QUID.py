"""Run calculations for the QUID non‑covalent interactions benchmark."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Iterable

from ase import Atoms
import h5py
import mlipx
from mlipx.abc import NodeWithCalculator
import numpy as np
import pytest
from tqdm import tqdm
import zntrack

from ml_peg.calcs.utils.utils import chdir
from ml_peg.models.get_models import load_models
from ml_peg.models.models import current_models


EV_TO_KCAL = 1.0 / 0.04336414
MODELS = load_models(current_models)


def get_paths() -> tuple[Path, Path, Path]:
    base = Path(__file__).parent
    data = base / "data" / "QUID.h5"
    cache = base / "cache_QUID"
    cache.mkdir(parents=True, exist_ok=True)
    if not data.exists():
        raise FileNotFoundError("QUID.h5 not found under data/")
    return base, cache, data


def _atoms_hash(numbers: np.ndarray, positions: np.ndarray) -> str:
    vec = np.concatenate([numbers.astype(np.int64), positions.reshape(-1)])
    return hashlib.md5(vec.tobytes()).hexdigest()


def energy_cached(numbers: np.ndarray, positions: np.ndarray, calc, cache_dir: Path, cache_tag: str) -> float:
    reset = getattr(calc, "reset", None)
    if callable(reset):
        reset()
    key = _atoms_hash(numbers, positions)
    fpath = cache_dir / f"{cache_tag}_{key}.json"
    if fpath.exists():
        try:
            return json.loads(fpath.read_text())["E"]
        except Exception:
            pass

    atoms = Atoms(numbers=numbers, positions=positions)
    atoms.calc = calc
    e = float(atoms.get_potential_energy())  # eV
    fpath.write_text(json.dumps({"E": e}))
    return e


def iter_dimer_ids(h5: h5py.File) -> Iterable[str]:
    # Top-level keys are equilibrium dimer labels
    return h5.keys()


class QUIDBenchmark(zntrack.Node):
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    def run(self):
        base, cache_dir, h5_path = get_paths()
        calc = self.model.get_calculator()

        out_dir = base / "outputs" / self.model_name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_csv = out_dir / "quid_results.csv"

        with h5py.File(h5_path, "r") as f, open(out_csv, "w", newline="") as fo:
            w = csv.writer(fo)
            w.writerow([
                "id",
                "ref_eint_eV",
                "calc_eint_eV",
                "ref_eint_kcal",
                "calc_eint_kcal",
                "error_kcal",
                "abs_error_kcal",
            ])
            for name in tqdm(list(iter_dimer_ids(f)), desc=f"QUID: {self.model_name}"):
                g = f[name]
                # Reference interaction energy: prefer CCSDT -> CCSD -> PBE0+MBD, else first available
                ref_grp = g["Eint"]
                ref_val = None
                for pref in ["CCSDT", "CCSD", "PBE0+MBD", "PBE0+XDM", "PBE+MBD"]:
                    if pref in ref_grp:
                        ref_val = float(ref_grp[pref][()])
                        break
                if ref_val is None:
                    # Fall back to SAPT if nothing else available
                    for k in ref_grp.keys():
                        ref_val = float(ref_grp[k][()])
                        break

                # Build geometries
                nums_d = g["atoms"]["dimer"][()]
                pos_d = g["positions"]["dimer"][()]
                nums_b = g["atoms"]["big_monomer"][()]
                pos_b = g["positions"]["big_monomer"][()]
                nums_s = g["atoms"]["small_monomer"][()]
                pos_s = g["positions"]["small_monomer"][()]

                e_d = energy_cached(nums_d, pos_d, calc, cache_dir, self.model_name)
                e_b = energy_cached(nums_b, pos_b, calc, cache_dir, self.model_name)
                e_s = energy_cached(nums_s, pos_s, calc, cache_dir, self.model_name)
                calc_eint_eV = e_d - e_b - e_s

                w.writerow([
                    name,
                    ref_val,
                    calc_eint_eV,
                    ref_val * EV_TO_KCAL,
                    calc_eint_eV * EV_TO_KCAL,
                    (calc_eint_eV - ref_val) * EV_TO_KCAL,
                    abs((calc_eint_eV - ref_val) * EV_TO_KCAL),
                ])


def build_project(repro: bool = False) -> None:
    project = mlipx.Project()
    nodes: dict[str, QUIDBenchmark] = {}
    for model_name, model in MODELS.items():
        with project.group(model_name):
            nodes[model_name] = QUIDBenchmark(model=model, model_name=model_name)
    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()
        for n in nodes.values():
            n.run()


def test_quid():
    build_project(repro=False)
