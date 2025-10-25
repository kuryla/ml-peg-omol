"""Run calculations for the MOR41 benchmark (closed‑shell reactions)."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

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


EV_TO_KCAL = 1.0 / 0.04336414

MODELS = load_models(current_models)


def get_dirs() -> tuple[Path, Path]:
    structures_dir = Path(__file__).parent / "data" / "geometries"
    cache_dir = Path(__file__).parent / "cache_MOR41"
    cache_dir.mkdir(parents=True, exist_ok=True)
    if not structures_dir.exists():
        raise FileNotFoundError("MOR41 geometries not found under data/geometries")
    return structures_dir, cache_dir


# Reference values (kcal/mol) from notebook (DLPNO‑CCSD(T)/CBS)
REFERENCE_VALUES = {
    "MOR41_01": -43.1,
    "MOR41_02": -46.6,
    "MOR41_03": -27.6,
    "MOR41_04": -62.5,
    "MOR41_05": 3.7,
    "MOR41_06": -23.2,
    "MOR41_07": -16.2,
    "MOR41_08": -17.2,
    "MOR41_09": -18.7,
    "MOR41_10": -22.6,
    "MOR41_11": 27.0,
    "MOR41_12": -29.8,
    "MOR41_13": -43.2,
    "MOR41_14": -52.0,
    "MOR41_15": -4.1,
    "MOR41_16": -39.8,
    "MOR41_17": -16.1,
    "MOR41_18": -34.2,
    "MOR41_19": -40.1,
    "MOR41_20": -30.2,
    "MOR41_21": -15.1,
    "MOR41_22": -35.9,
    "MOR41_23": -55.0,
    "MOR41_24": -41.6,
    "MOR41_25": -45.9,
    "MOR41_26": -36.4,
    "MOR41_27": -21.8,
    "MOR41_28": -36.3,
    "MOR41_29": -28.3,
    "MOR41_30": -14.0,
    "MOR41_31": -29.9,
    "MOR41_32": -1.8,
    "MOR41_33": -10.7,
    "MOR41_34": -25.6,
    "MOR41_35": -30.9,
    "MOR41_36": -39.8,
    "MOR41_37": -14.0,
    "MOR41_38": -64.4,
    "MOR41_39": -63.9,
    "MOR41_40": -65.8,
    "MOR41_41": -3.2,
}


REACTION_TYPES = {
    "σ-donative complexation": [1, 2, 3, 4, 5, 6, 7, 8, 9, 16, 24, 25],
    "π-acceptor complexation": [12, 13, 14],
    "oxidative addition": [10, 11, 17, 18, 19, 20, 21, 22, 23],
    "ligand exchange": [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38],
    "insertion": [15],
    "dimerization": [39],
    "metathesis": [41],
    "special": [40],
}


def create_reaction_mapping() -> dict[str, dict]:
    # Extracted from the MOR41 notebook
    return {
        "MOR41_01": {"reactants": [("ED01", 1), ("CO", 1)], "products": [("PR01", 1)], "description": "σ-donative complexation"},
        "MOR41_02": {"reactants": [("ED02", 1), ("CO", 1)], "products": [("PR02", 1)], "description": "σ-donative complexation"},
        "MOR41_03": {"reactants": [("ED03", 1), ("CO", 1)], "products": [("PR03", 1)], "description": "σ-donative complexation"},
        "MOR41_04": {"reactants": [("ED04", 1), ("CO", 1)], "products": [("PR04", 1)], "description": "σ-donative complexation"},
        "MOR41_05": {"reactants": [("ED05", 1), ("CO", 1)], "products": [("PR05", 1)], "description": "σ-donative complexation"},
        "MOR41_06": {"reactants": [("ED01", 1), ("H2", 1)], "products": [("PR06", 1)], "description": "σ-donative complexation"},
        "MOR41_07": {"reactants": [("ED07", 1), ("H2", 1)], "products": [("PR07", 1)], "description": "σ-donative complexation"},
        "MOR41_08": {"reactants": [("ED08", 1), ("H2", 1)], "products": [("PR08", 1)], "description": "σ-donative complexation"},
        "MOR41_09": {"reactants": [("ED09", 1), ("H2", 1)], "products": [("PR09", 1)], "description": "σ-donative complexation"},
        "MOR41_10": {"reactants": [("ED10", 1), ("C3H8", 1)], "products": [("PR10", 1)], "description": "oxidative addition"},
        "MOR41_11": {"reactants": [("ED11", 1), ("C2H6", 1)], "products": [("PR11", 1)], "description": "oxidative addition"},
        "MOR41_12": {"reactants": [("ED01", 1), ("C2H4", 1)], "products": [("PR12", 1)], "description": "π-acceptor complexation"},
        "MOR41_13": {"reactants": [("ED13", 1), ("C2H4", 1)], "products": [("PR13", 1)], "description": "π-acceptor complexation"},
        "MOR41_14": {"reactants": [("ED14", 1), ("COD", 1)], "products": [("PR14", 1)], "description": "π-acceptor complexation"},
        "MOR41_15": {"reactants": [("ED15", 1), ("CO2", 1)], "products": [("PR15", 1)], "description": "insertion reaction"},
        "MOR41_16": {"reactants": [("ED16a", 1), ("ED16b", 1)], "products": [("PR16", 1)], "description": "σ-donative complexation"},
        "MOR41_17": {"reactants": [("ED17", 1), ("MeOH", 1)], "products": [("PR17", 1)], "description": "oxidative addition"},
        "MOR41_18": {"reactants": [("ED18", 1), ("MeI", 1)], "products": [("PR18", 1)], "description": "oxidative addition"},
        "MOR41_19": {"reactants": [("ED18", 1), ("AcI", 1)], "products": [("PR19", 1)], "description": "oxidative addition"},
        "MOR41_20": {"reactants": [("ED18", 1), ("AcCl", 1)], "products": [("PR20", 1)], "description": "oxidative addition"},
        "MOR41_21": {"reactants": [("ED21", 1), ("MeI", 1)], "products": [("PR21", 1)], "description": "oxidative addition"},
        "MOR41_22": {"reactants": [("ED22", 1), ("I2", 1)], "products": [("PR22", 1)], "description": "oxidative addition"},
        "MOR41_23": {"reactants": [("ED18", 1), ("I2", 1)], "products": [("PR23", 1)], "description": "oxidative addition"},
        "MOR41_24": {"reactants": [("ED24", 1), ("PCy3", 1)], "products": [("PR24", 1)], "description": "σ-donative complexation"},
        "MOR41_25": {"reactants": [("ED25", 1), ("PCy3", 1)], "products": [("PR25", 1)], "description": "σ-donative complexation"},
        "MOR41_26": {"reactants": [("ED26", 1), ("PMe3", 2)], "products": [("PR26", 1), ("COD", 1)], "description": "ligand exchange"},
        "MOR41_27": {"reactants": [("ED26", 1), ("ED27", 1)], "products": [("PR27", 1), ("COD", 1)], "description": "ligand exchange"},
        "MOR41_28": {"reactants": [("ED26", 1), ("ED28", 1)], "products": [("PR28", 1), ("COD", 1)], "description": "ligand exchange"},
        "MOR41_29": {"reactants": [("ED29", 1), ("ED28", 1)], "products": [("PR29", 1), ("COD", 1)], "description": "ligand exchange"},
        "MOR41_30": {"reactants": [("ED29", 1), ("ED30", 1)], "products": [("PR30", 1), ("COD", 1)], "description": "ligand exchange"},
        "MOR41_31": {"reactants": [("ED29", 1), ("ED31", 1)], "products": [("PR31", 1), ("COD", 1)], "description": "ligand exchange"},
        "MOR41_32": {"reactants": [("ED32", 1), ("COD", 1)], "products": [("PR32", 1), ("C2H4", 2)], "description": "ligand exchange"},
        "MOR41_33": {"reactants": [("ED33", 1), ("PhOH", 1)], "products": [("PR33", 1), ("C2H4", 1)], "description": "ligand exchange"},
        "MOR41_34": {"reactants": [("ED33", 1), ("PhSH", 1)], "products": [("PR34", 1), ("C2H4", 1)], "description": "ligand exchange"},
        "MOR41_35": {"reactants": [("ED33", 1), ("PhSeH", 1)], "products": [("PR35", 1), ("C2H4", 1)], "description": "ligand exchange"},
        "MOR41_36": {"reactants": [("ED36", 1), ("I2", 1)], "products": [("PR36", 1), ("C2H4", 1)], "description": "ligand exchange"},
        "MOR41_37": {"reactants": [("ED37", 1), ("MeCN", 3)], "products": [("PR37", 1), ("Bz", 1)], "description": "ligand exchange"},
        "MOR41_38": {"reactants": [("ED37", 1), ("PMe3", 3)], "products": [("PR38", 1), ("Bz", 1)], "description": "ligand exchange"},
        "MOR41_39": {"reactants": [("ED39", 2)], "products": [("PR39", 1)], "description": "dimerization"},
        "MOR41_40": {"reactants": [("ED40a", 1), ("ED40b", 1)], "products": [("PR40", 1)], "description": "Tebbe reagent formation"},
        "MOR41_41": {"reactants": [("ED41", 1)], "products": [("PR41", 1)], "description": "ring opening metathesis"},
    }


def _atoms_hash(atoms) -> str:
    vec = np.concatenate([atoms.get_positions().ravel(), atoms.get_cell().ravel(), atoms.get_atomic_numbers()])
    return hashlib.md5(vec.tobytes()).hexdigest()


def energy_cached(atoms, calc, cache_dir: Path, cache_tag: str) -> float:
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
    e = float(cpy.get_potential_energy())
    fpath.write_text(json.dumps({"E": e}))
    return e


def load_molecule(structures_dir: Path, name: str):
    # Handle known naming discrepancy: ED41 is stored as EDPR41 in dataset
    alias = {"ED41": "EDPR41"}
    mol_dir = structures_dir / (alias.get(name, name))
    xyz = mol_dir / "mol.xyz"
    coord = mol_dir / "coord"
    if xyz.exists():
        return read(str(xyz))
    if coord.exists():
        return read(str(coord))
    raise FileNotFoundError(f"No coordinates for {name}")


class MOR41Benchmark(zntrack.Node):
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    def run(self):
        structures_dir, cache_dir = get_dirs()
        mapping = create_reaction_mapping()
        calc = self.model.get_calculator()

        out_dir = Path(__file__).parent / "outputs" / self.model_name
        out_dir.mkdir(parents=True, exist_ok=True)

        with open(out_dir / "mor41_results.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow([
                "reaction_id",
                "reaction_number",
                "reaction_type",
                "description",
                "ref_energy",
                "calc_energy",
                "error",
                "abs_error",
                "max_atoms",
            ])
            for rid, data in tqdm(mapping.items(), desc=f"MOR41: {self.model_name}"):
                coeffs: list[tuple[float, str]] = []
                for mol, sto in data["reactants"]:
                    coeffs.append((-float(sto), mol))
                for mol, sto in data["products"]:
                    coeffs.append((float(sto), mol))

                e_comb = 0.0
                max_atoms = 0
                for coeff, mol in coeffs:
                    atoms = load_molecule(structures_dir, mol)
                    e = energy_cached(atoms, calc, cache_dir, self.model_name)
                    e_comb += coeff * e
                    max_atoms = max(max_atoms, len(atoms))

                calc_kcal = e_comb * EV_TO_KCAL
                ref = REFERENCE_VALUES[rid]
                err = calc_kcal - ref
                num = int(rid.split("_")[-1])
                rtype = next((k for k, lst in REACTION_TYPES.items() if num in lst), "")
                w.writerow([rid, num, rtype, data["description"], ref, calc_kcal, err, abs(err), max_atoms])


def build_project(repro: bool = False) -> None:
    project = mlipx.Project()
    nodes: dict[str, MOR41Benchmark] = {}
    for model_name, model in MODELS.items():
        with project.group(model_name):
            nodes[model_name] = MOR41Benchmark(model=model, model_name=model_name)
    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()
        for n in nodes.values():
            n.run()


def test_mor41():
    build_project(repro=False)
