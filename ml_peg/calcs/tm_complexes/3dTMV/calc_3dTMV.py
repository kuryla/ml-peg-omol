"""Run calculations for the 3dTMV benchmark (vertical ionization potentials)."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from ase.io import read
import mlipx
from mlipx.abc import NodeWithCalculator
import numpy as np
import pytest
from tqdm import tqdm
import zntrack

from ml_peg.calcs.utils.utils import chdir
from ml_peg.models.get_models import load_models
from ml_peg.models.models import current_models


KCAL_TO_EV = 0.04336414
EV_TO_KCAL = 1.0 / KCAL_TO_EV

# Models selected via CLI/pytest flag
MODELS = load_models(current_models)


def get_dirs() -> tuple[Path, Path]:
    """Locate dataset and cache directories for 3dTMV.

    Returns
    -------
    tuple[Path, Path]
        (structures_dir, cache_dir)
    """
    # Prefer packaged data
    candidates = [Path(__file__).parent / "data" / "3dtmv_structures"]
    env = os.environ.get("THREEDTMV_DIR") or os.environ.get("3DTMV_DIR")
    if env:
        candidates.append(Path(env))
    candidates += [
        Path("/Users/Lilyes/Documents/Work/notebooks_mp/3dtmv/3dtmv_structures"),
        Path("/Users/Lilyes/Documents/GitHub/benchmarks-mp/3dtmv/3dtmv_structures"),
    ]

    structures_dir = None
    for cand in candidates:
        if (cand).exists():
            structures_dir = cand
            break
    if structures_dir is None:
        raise FileNotFoundError(
            "3dTMV structures not found. Set 3DTMV_DIR to path containing 1/struc.xyz ... 28/struc.xyz."
        )

    cache_env = os.environ.get("3DTMV_CACHE_DIR")
    cache_dir = Path(cache_env) if cache_env else Path(__file__).parent / "cache_3dTMV"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return structures_dir, cache_dir


# Reference IPs (kcal/mol) from the notebook
REFERENCE_IPS = {
    1: 188.4, 2: 158.3, 3: 119.6, 4: 152.3, 5: 142.2, 6: 315.9,
    7: 191.1, 8: 259.6, 9: 276.2, 10: 284.1, 11: 198.5, 12: 230.3,
    13: 120.9, 14: 148.1, 15: 140.4, 16: 164.1, 17: 130.9, 18: 136.3,
    19: 300.7, 20: 186.4, 21: 125.3, 22: 161.2,
    23: 198.9, 24: 166.0, 25: 215.8, 26: 192.9, 27: 68.6, 28: 43.6,
}

# Molecular charge/multiplicity for neutral (in) and oxidized (ox) states
MOLECULAR_DATA = {
    1: {"charge_ox": 1, "charge_in": 0, "mult_ox": 2, "mult_in": 1, "subset": "SR"},
    2: {"charge_ox": 1, "charge_in": 0, "mult_ox": 1, "mult_in": 2, "subset": "SR"},
    3: {"charge_ox": 1, "charge_in": 0, "mult_ox": 4, "mult_in": 3, "subset": "SR"},
    4: {"charge_ox": 1, "charge_in": 0, "mult_ox": 2, "mult_in": 1, "subset": "SR"},
    5: {"charge_ox": 1, "charge_in": 0, "mult_ox": 2, "mult_in": 1, "subset": "SR"},
    6: {"charge_ox": 2, "charge_in": 1, "mult_ox": 2, "mult_in": 1, "subset": "SR"},
    7: {"charge_ox": 1, "charge_in": 0, "mult_ox": 2, "mult_in": 1, "subset": "SR"},
    8: {"charge_ox": 2, "charge_in": 1, "mult_ox": 2, "mult_in": 1, "subset": "SR"},
    9: {"charge_ox": 2, "charge_in": 1, "mult_ox": 2, "mult_in": 1, "subset": "SR"},
    10: {"charge_ox": 2, "charge_in": 1, "mult_ox": 2, "mult_in": 1, "subset": "SR"},
    11: {"charge_ox": 2, "charge_in": 1, "mult_ox": 1, "mult_in": 2, "subset": "SR"},
    12: {"charge_ox": 2, "charge_in": 1, "mult_ox": 1, "mult_in": 2, "subset": "SR"},
    13: {"charge_ox": 1, "charge_in": 0, "mult_ox": 1, "mult_in": 2, "subset": "SR/MR"},
    14: {"charge_ox": 1, "charge_in": 0, "mult_ox": 2, "mult_in": 3, "subset": "SR/MR"},
    15: {"charge_ox": 1, "charge_in": 0, "mult_ox": 2, "mult_in": 1, "subset": "SR/MR"},
    16: {"charge_ox": 1, "charge_in": 0, "mult_ox": 2, "mult_in": 1, "subset": "SR/MR"},
    17: {"charge_ox": 1, "charge_in": 0, "mult_ox": 2, "mult_in": 1, "subset": "SR/MR"},
    18: {"charge_ox": 1, "charge_in": 0, "mult_ox": 2, "mult_in": 1, "subset": "SR/MR"},
    19: {"charge_ox": 2, "charge_in": 1, "mult_ox": 2, "mult_in": 1, "subset": "SR/MR"},
    20: {"charge_ox": 1, "charge_in": 0, "mult_ox": 3, "mult_in": 2, "subset": "SR/MR"},
    21: {"charge_ox": 1, "charge_in": 0, "mult_ox": 3, "mult_in": 2, "subset": "SR/MR"},
    22: {"charge_ox": 1, "charge_in": 0, "mult_ox": 2, "mult_in": 1, "subset": "SR/MR"},
    23: {"charge_ox": 1, "charge_in": 0, "mult_ox": 2, "mult_in": 1, "subset": "MR"},
    24: {"charge_ox": 1, "charge_in": 0, "mult_ox": 3, "mult_in": 4, "subset": "MR"},
    25: {"charge_ox": 1, "charge_in": 0, "mult_ox": 3, "mult_in": 6, "subset": "MR"},
    26: {"charge_ox": 1, "charge_in": 0, "mult_ox": 2, "mult_in": 1, "subset": "MR"},
    27: {"charge_ox": 0, "charge_in": -1, "mult_ox": 2, "mult_in": 3, "subset": "MR"},
    28: {"charge_ox": 0, "charge_in": -1, "mult_ox": 1, "mult_in": 2, "subset": "MR"},
}


def _atoms_hash(atoms) -> str:
    vec = np.concatenate(
        [
            atoms.get_positions().ravel(),
            atoms.get_cell().ravel(),
            atoms.get_atomic_numbers(),
            np.array([atoms.info.get("charge", 0)]),
            np.array([atoms.info.get("spin", 1)]),
        ]
    )
    return hashlib.md5(vec.tobytes()).hexdigest()


def energy_cached(atoms, calc, cache_dir: Path, cache_tag: str) -> float:
    """Compute potential energy with per-structure caching (eV)."""
    reset = getattr(calc, "reset", None)
    if callable(reset):
        reset()

    key = _atoms_hash(atoms)
    fpath = cache_dir / f"{cache_tag}_{key}.json"
    if fpath.exists():
        try:
            return json.loads(fpath.read_text())["E"]
        except Exception:
            pass

    cpy = atoms.copy()
    cpy.calc = calc

    # Pass charge/multiplicity to calculators that support it
    charge = int(cpy.info.get("charge", 0))
    mult = int(cpy.info.get("spin", 1))
    if hasattr(calc, "set_charge"):
        try:
            calc.set_charge(charge)
        except Exception:
            pass
    if hasattr(calc, "set_multiplicity"):
        try:
            calc.set_multiplicity(mult)
        except Exception:
            pass

    e = float(cpy.get_potential_energy())
    fpath.write_text(json.dumps({"E": e}))
    return e


def load_structure(structures_dir: Path, complex_id: int) -> Optional["ase.Atoms"]:
    structure_path = structures_dir / str(complex_id) / "struc.xyz"
    if not structure_path.exists():
        return None
    return read(str(structure_path))


def prepare_atoms_for_calculation(atoms, charge: int, multiplicity: int):
    cpy = atoms.copy()
    cpy.info["charge"] = int(charge)
    cpy.info["spin"] = int(multiplicity)
    return cpy


@dataclass
class TMVResult:
    complex_id: int
    subset: str
    IP_ref: float
    IP_model: float
    error: float
    abs_error: float
    atoms_count: int
    charge_in: int
    charge_ox: int
    mult_in: int
    mult_ox: int


class TMVBenchmark(zntrack.Node):
    """3dTMV benchmark node."""

    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    def run(self):
        structures_dir, cache_dir = get_dirs()
        calc = self.model.get_calculator()

        out_dir = Path(__file__).parent / "outputs" / self.model_name
        out_dir.mkdir(parents=True, exist_ok=True)

        results: list[TMVResult] = []
        for complex_id in tqdm(range(1, 29), desc=f"3dTMV: {self.model_name}"):
            atoms = load_structure(structures_dir, complex_id)
            if atoms is None:
                continue
            data = MOLECULAR_DATA[complex_id]
            subset = data["subset"]

            atoms_in = prepare_atoms_for_calculation(
                atoms, data["charge_in"], data["mult_in"]
            )
            atoms_ox = prepare_atoms_for_calculation(
                atoms, data["charge_ox"], data["mult_ox"]
            )

            e_in = energy_cached(atoms_in, calc, cache_dir, self.model_name)
            e_ox = energy_cached(atoms_ox, calc, cache_dir, self.model_name)
            ip_kcal = (e_ox - e_in) * EV_TO_KCAL
            ref = REFERENCE_IPS[complex_id]
            err = ip_kcal - ref

            results.append(
                TMVResult(
                    complex_id=complex_id,
                    subset=subset,
                    IP_ref=ref,
                    IP_model=ip_kcal,
                    error=err,
                    abs_error=abs(err),
                    atoms_count=len(atoms),
                    charge_in=data["charge_in"],
                    charge_ox=data["charge_ox"],
                    mult_in=data["mult_in"],
                    mult_ox=data["mult_ox"],
                )
            )

        # Write detailed CSV
        with open(out_dir / "tmv_results.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(
                [
                    "complex_id",
                    "subset",
                    "IP_ref",
                    "IP_model",
                    "error",
                    "abs_error",
                    "atoms_count",
                    "charge_in",
                    "charge_ox",
                    "mult_in",
                    "mult_ox",
                ]
            )
            for r in results:
                w.writerow(
                    [
                        r.complex_id,
                        r.subset,
                        r.IP_ref,
                        r.IP_model,
                        r.error,
                        r.abs_error,
                        r.atoms_count,
                        r.charge_in,
                        r.charge_ox,
                        r.mult_in,
                        r.mult_ox,
                    ]
                )


def build_project(repro: bool = False) -> None:
    project = mlipx.Project()
    nodes: dict[str, TMVBenchmark] = {}
    for model_name, model in MODELS.items():
        with project.group(model_name):
            nodes[model_name] = TMVBenchmark(model=model, model_name=model_name)
    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()
        for node in nodes.values():
            node.run()


def test_3dtmv():
    build_project(repro=False)

