"""Run calculations for the ROST61 benchmark (open-shell TM complexes)."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

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


# Models selected via CLI/pytest flag (see conftest.py)
MODELS = load_models(current_models)

# Energy conversion constants (as used in the notebook)
KCAL_TO_EV = 0.04336414
EV_TO_KCAL = 1.0 / KCAL_TO_EV


def get_dirs() -> tuple[Path, Path]:
    """Locate dataset and cache directories for ROST61.

    Returns
    -------
    tuple[Path, Path]
        (geometries_dir, cache_dir)
    """
    # Dataset dir precedence: local packaged data, env var, known user paths
    candidates = [Path(__file__).parent / "data" / "ROST61"]
    env = os.environ.get("ROST61_DIR")
    if env:
        candidates.append(Path(env))
    candidates.append(Path("/Users/Lilyes/Documents/Work/notebooks_mp/rost61/ROST61"))
    candidates.append(Path("/Users/Lilyes/Documents/GitHub/benchmarks-mp/rost61/ROST61"))

    geometries_dir = None
    for cand in candidates:
        if (cand).exists():
            geometries_dir = cand
            break
    if geometries_dir is None:
        raise FileNotFoundError(
            "ROST61 geometry directory not found. Set ROST61_DIR to the path containing m1, m2, ... and mol.xyz files."
        )

    # Cache dir: prefer local cache folder under benchmark
    cache_env = os.environ.get("ROST61_CACHE_DIR")
    if cache_env:
        cache_dir = Path(cache_env)
    else:
        cache_dir = Path(__file__).parent / "cache_ROST61"
    cache_dir.mkdir(parents=True, exist_ok=True)

    return geometries_dir, cache_dir


REFERENCE_VALUES: dict[str, float] = {
    "R1": -40.05,
    "R2": -37.71,
    "R3": -12.59,
    "R4": -17.04,
    "R5": -10.66,
    "R6": -12.66,
    "R7": -14.73,
    "R8": -10.51,
    "R9": -12.76,
    "R10": -21.20,
    "R11": -8.44,
    "R12": -2.79,
    "R13": -178.96,
    "R14": -42.33,
    "R15": -7.17,
    "R16": -193.43,
    "R17": -200.10,
    "R18": -203.09,
    "R19": -4.94,
    "R20": -27.02,
    "R21": -25.09,
    "R22": -46.70,
    "R23": -115.35,
    "R24": -46.37,
    "R25": -19.39,
    "R26": -47.73,
    "R27": -25.86,
    "R28": -28.77,
    "R29": -11.13,
    "R30": -141.32,
    "R31": -63.73,
    "R32": -49.75,
    "R33": -10.64,
    "R34": -0.65,
    "R35": -5.46,
    "R36": -49.04,
    "R37": -40.62,
    "R38": -27.70,
    "R39": -44.99,
    "R40": -24.51,
    "R41": -5.96,
    "R42": -26.96,
    "R43": -80.58,
    "R44": -66.15,
    "R45": -42.83,
    "R46": -33.27,
    "R47": -58.47,
    "R48": -4.24,
    "R49": -35.26,
    "R50": -2.74,
    "R51": -29.77,
    "R52": -42.72,
    "R53": -30.78,
    "R54": -12.33,
    "R55": -4.61,
    "R56": -13.65,
    "R57": -31.89,
    "R58": -39.30,
    "R59": -29.89,
    "R60": -66.00,
    "R61": -69.48,
}


def create_reaction_mapping() -> dict[str, dict]:
    """Notebook metadata (metal, description) for each reaction.

    Stoichiometry is sourced from the `.res` file; this mapping is used only
    for metadata (metal and description) and reaction-type classification.
    """
    reaction_mapping = {
        "R1": {
            "reactants": [("m1", 1), ("m2", 1)],
            "products": [("m3", 1)],
            "metal": "Ti",
            "description": "1-electron oxidative addition",
        },
        "R2": {
            "reactants": [("m3", 1)],
            "products": [("m4", 1)],
            "metal": "Ti",
            "description": "1-electron reductive elimination",
        },
        "R3": {
            "reactants": [("m1", 1), ("m5", 1)],
            "products": [("m4", 1), ("m6", 1)],
            "metal": "Ti",
            "description": "σ-bond-metathesis",
        },
        "R4": {
            "reactants": [("m7", 1), ("m8", 1)],
            "products": [("m9", 1)],
            "metal": "Cr",
            "description": "Ligand coordination",
        },
        "R5": {
            "reactants": [("m9", 1), ("m8", 1)],
            "products": [("m10", 1)],
            "metal": "Cr",
            "description": "Ligand coordination",
        },
        "R6": {
            "reactants": [("m10", 1)],
            "products": [("m11", 1)],
            "metal": "Cr",
            "description": "Cycloaddition",
        },
        "R7": {
            "reactants": [("m11", 1), ("m8", 1)],
            "products": [("m12", 1)],
            "metal": "Cr",
            "description": "Ligand coordination",
        },
        "R8": {
            "reactants": [("m12", 1)],
            "products": [("m13", 1)],
            "metal": "Cr",
            "description": "Migratory insertion",
        },
        "R9": {
            "reactants": [("m14", 1)],
            "products": [("m13", 1)],
            "metal": "Cr",
            "description": "Migratory insertion",
        },
        "R10": {
            "reactants": [("m14", 1)],
            "products": [("m15", 1)],
            "metal": "Cr",
            "description": "Reductive elimination",
        },
        "R11": {
            "reactants": [("m13", 1)],
            "products": [("m15", 1)],
            "metal": "Cr",
            "description": "Intramolecular H atom abstraction",
        },
        "R12": {
            "reactants": [("m9", 1), ("m16", 1)],
            "products": [("m15", 1), ("m8", 1)],
            "metal": "Cr",
            "description": "Ligand exchange",
        },
        "R13": {
            "reactants": [("m19", 1), ("m20", 1), ("m21", 1)],
            "products": [("m22", 1), ("m23", 1)],
            "metal": "Co",
            "description": "Reductive ligand exchange",
        },
        "R14": {
            "reactants": [("m24", 1), ("m25", 1)],
            "products": [("m26", 1)],
            "metal": "Cr",
            "description": "Migratory insertion",
        },
        "R15": {
            "reactants": [("m26", 1), ("m27", 1)],
            "products": [("m28", 1), ("m29", 1)],
            "metal": "Cr",
            "description": "σ-bond-metathesis",
        },
        "R16": {
            "reactants": [("m30", 1), ("m31", 1)],
            "products": [("m32", 1)],
            "metal": "V",
            "description": "Ligand coordination",
        },
        "R17": {
            "reactants": [("m33", 1), ("m31", 1)],
            "products": [("m34", 1)],
            "metal": "Ta",
            "description": "Ligand coordination",
        },
        "R18": {
            "reactants": [("m35", 1), ("m31", 1)],
            "products": [("m36", 1)],
            "metal": "Nb",
            "description": "Ligand coordination",
        },
        "R19": {
            "reactants": [("m37", 1), ("m38", 0.5)],
            "products": [("m39", 1)],
            "metal": "Cr",
            "description": "1-electron oxidative addition",
        },
        "R20": {
            "reactants": [("m40", 1), ("m41", 1)],
            "products": [("m42", 1)],
            "metal": "Mo",
            "description": "Ligand coordination",
        },
        "R21": {
            "reactants": [("m43", 1), ("m44", 1)],
            "products": [("m45", 1), ("m46", 1)],
            "metal": "Zr",
            "description": "1-electron oxidative addition",
        },
        "R22": {
            "reactants": [("m43", 1), ("m46", 1)],
            "products": [("m47", 1)],
            "metal": "Zr",
            "description": "1-electron oxidative addition",
        },
        "R23": {
            "reactants": [("m49", 1), ("m50", 1)],
            "products": [("m51", 1)],
            "metal": "Cu",
            "description": "Ligand exchange",
        },
        "R24": {
            "reactants": [("m51", 1), ("m52", 1)],
            "products": [("m53", 1)],
            "metal": "Cu",
            "description": "Cycloaddition",
        },
        "R25": {
            "reactants": [("m54", 1), ("m55", 2)],
            "products": [("m56", 1), ("m17", 2)],
            "metal": "Cr",
            "description": "Ligand exchange",
        },
        "R26": {
            "reactants": [("m56", 1), ("m57", 1)],
            "products": [("m58", 1), ("m18", 1), ("m59", 1)],
            "metal": "Cr",
            "description": "Ligand exchange",
        },
        "R27": {
            "reactants": [("m60", 1), ("m48", 1), ("m61", 1), ("m62", 1)],
            "products": [("m25", 1), ("m48", 1), ("m62", 1)],
            "metal": "Ru",
            "description": "Ligand exchange",
        },
        "R28": {
            "reactants": [("m63", 1), ("m48", 1), ("m64", 1), ("m62", 1)],
            "products": [("m25", 1), ("m48", 1), ("m62", 1)],
            "metal": "Os",
            "description": "Ligand exchange",
        },
        "R29": {
            "reactants": [("m65", 1)],
            "products": [("m66", 1)],
            "metal": "Ti",
            "description": "Ligand coordination",
        },
        "R30": {
            "reactants": [("m67", 1), ("m41", 1), ("m65", 1)],
            "products": [("m29", 1)],
            "metal": "Ti",
            "description": "Ligand coordination",
        },
        "R31": {
            "reactants": [("m68", 1), ("m31", 1), ("m69", 1)],
            "products": [("m32", 1)],
            "metal": "Ti",
            "description": "Ligand coordination",
        },
        "R32": {
            "reactants": [("m70", 1), ("m20", 1), ("m71", 1), ("m41", 1)],
            "products": [("m32", 1), ("m41", 1)],
            "metal": "W",
            "description": "Ligand exchange",
        },
        "R33": {
            "reactants": [("m72", 1), ("m73", 1)],
            "products": [("m74", 1)],
            "metal": "Ti",
            "description": "1-electron reductive elimination",
        },
        "R34": {
            "reactants": [("m75", 1), ("m76", 1)],
            "products": [("m77", 1), ("m17", 1)],
            "metal": "Ti",
            "description": "Ligand exchange",
        },
        "R35": {
            "reactants": [("m78", 1)],
            "products": [("m77", 1)],
            "metal": "Ti",
            "description": "Radical translocation",
        },
        "R36": {
            "reactants": [("m79", 1), ("m80", 4)],
            "products": [("m81", 1), ("m21", 4)],
            "metal": "Cu",
            "description": "Ligand exchange",
        },
        "R37": {
            "reactants": [("m82", 1), ("m83", 1)],
            "products": [("m84", 1)],
            "metal": "Ru",
            "description": "Ligand coordination",
        },
        "R38": {
            "reactants": [("m85", 1)],
            "products": [("m86", 1), ("m87", 1)],
            "metal": "Ni",
            "description": "1-electron oxidative addition",
        },
        "R39": {
            "reactants": [("m88", 1), ("m89", 1)],
            "products": [("m90", 1), ("m91", 1)],
            "metal": "V",
            "description": "σ-bond-metathesis",
        },
        "R40": {
            "reactants": [("m92", 1), ("m93", 1)],
            "products": [("m94", 1), ("m91", 1)],
            "metal": "V",
            "description": "σ-bond-metathesis",
        },
        "R41": {
            "reactants": [("m94", 1), ("m95", 1)],
            "products": [("m96", 1), ("m93", 1)],
            "metal": "V",
            "description": "σ-bond-metathesis",
        },
        "R42": {
            "reactants": [("m96", 1), ("m97", 1)],
            "products": [("m98", 1), ("m99", 0.5)],
            "metal": "V",
            "description": "1-electron oxidative addition",
        },
        "R43": {
            "reactants": [("m100", 1), ("m101", 1), ("m102", 1), ("m17", 2)],
            "products": [("m103", 1), ("m104", 1), ("m105", 1)],
            "metal": "Mo",
            "description": "Ligand exchange",
        },
        "R44": {
            "reactants": [("m106", 1), ("m102", 1), ("m107", 1), ("m101", 1)],
            "products": [("m108", 1), ("m104", 1), ("m105", 1)],
            "metal": "W",
            "description": "Ligand exchange",
        },
        "R45": {
            "reactants": [("m109", 1), ("m83", 1)],
            "products": [("m110", 1)],
            "metal": "Os",
            "description": "Other",
        },
        "R46": {
            "reactants": [("m111", 1), ("m112", 1)],
            "products": [("m113", 1)],
            "metal": "Zr",
            "description": "Migratory insertion",
        },
        "R47": {
            "reactants": [("m114", 1), ("m115", 1)],
            "products": [("m116", 1), ("m41", 2), ("m107", 1)],
            "metal": "Nb",
            "description": "Ligand exchange",
        },
        "R48": {
            "reactants": [("m117", 1), ("m91", 1)],
            "products": [("m118", 1)],
            "metal": "Ti",
            "description": "σ-bond-metathesis",
        },
        "R49": {
            "reactants": [("m119", 1), ("m120", 1)],
            "products": [("m121", 1), ("m41", 1)],
            "metal": "Ti",
            "description": "Ligand exchange",
        },
        "R50": {
            "reactants": [("m121", 1)],
            "products": [("m122", 1), ("m91", 1)],
            "metal": "Ti",
            "description": "σ-bond-metathesis",
        },
        "R51": {
            "reactants": [("m122", 1), ("m123", 1)],
            "products": [("m124", 1)],
            "metal": "Ti",
            "description": "σ-bond-metathesis",
        },
        "R52": {
            "reactants": [("m125", 1)],
            "products": [("m126", 1)],
            "metal": "Rh",
            "description": "Reductive elimination",
        },
        "R53": {
            "reactants": [("m127", 1), ("m120", 1)],
            "products": [("m128", 1), ("m41", 1)],
            "metal": "Zr",
            "description": "Ligand exchange",
        },
        "R54": {
            "reactants": [("m129", 1), ("m130", 1)],
            "products": [("m131", 1)],
            "metal": "Pd",
            "description": "Oxidative addition",
        },
        "R55": {
            "reactants": [("m132", 1)],
            "products": [("m133", 1)],
            "metal": "Tc",
            "description": "Other",
        },
        "R56": {
            "reactants": [("m134", 1), ("m135", 1)],
            "products": [("m136", 1), ("m137", 1)],
            "metal": "Ag",
            "description": "Ligand exchange",
        },
        "R57": {
            "reactants": [("m138", 1)],
            "products": [("m139", 1)],
            "metal": "Rh",
            "description": "Ligand coordination",
        },
        "R58": {
            "reactants": [("m140", 1)],
            "products": [("m141", 1)],
            "metal": "Ir",
            "description": "Ligand coordination",
        },
        "R59": {
            "reactants": [("m142", 1), ("m143", 1)],
            "products": [("m144", 1)],
            "metal": "Ir",
            "description": "Ligand coordination",
        },
        "R60": {
            "reactants": [("m145", 1), ("m146", 1)],
            "products": [("m147", 1)],
            "metal": "Y",
            "description": "Ligand coordination",
        },
        "R61": {
            "reactants": [("m148", 1), ("m149", 1)],
            "products": [("m150", 1)],
            "metal": "Sc",
            "description": "Ligand coordination",
        },
    }
    return reaction_mapping


def parse_res_file(res_path: Path) -> list[tuple[str, list[tuple[float, str]]]]:
    """Parse ROST61 `.res` file into per-reaction coefficient lists.

    Each line encodes a reaction's linear combination of species energies:
    `$tmer mX/$f ... x c1 c2 ... $w ref idx`

    Returns
    -------
    list[tuple[str, list[tuple[float, str]]]]
        List of (reaction_id, [(coeff, molecule), ...]). Negative coeffs are
        reactants, positive coeffs are products.
    """
    results: list[tuple[str, list[tuple[float, str]]]] = []
    with open(res_path, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0] != "$tmer":
                continue
            try:
                x_idx = parts.index("x")
                w_idx = parts.index("$w")
            except ValueError:
                continue
            mol_specs = parts[1:x_idx]
            coeff_specs = parts[x_idx + 1 : w_idx]
            # ref = float(parts[w_idx + 1])  # not used (we have REFERENCE_VALUES)
            idx = int(float(parts[w_idx + 2]))

            mol_labels = [spec.split("/")[0] for spec in mol_specs]
            coeffs = [float(c) for c in coeff_specs]
            if len(mol_labels) != len(coeffs):
                raise ValueError(
                    f"Malformed .res line (mismatch species/coeffs): {line}"
                )
            results.append((f"R{idx}", list(zip(coeffs, mol_labels, strict=True))))
    return results


METAL_CLASSIFICATION: dict[str, str] = {
    # 3d
    "Ti": "3d",
    "Cr": "3d",
    "Mn": "3d",
    "Fe": "3d",
    "Co": "3d",
    "Ni": "3d",
    "Cu": "3d",
    # 4d
    "Zr": "4d",
    "Nb": "4d",
    "Mo": "4d",
    "Tc": "4d",
    "Ru": "4d",
    "Rh": "4d",
    "Pd": "4d",
    "Ag": "4d",
    # 5d
    "Ta": "5d",
    "W": "5d",
    "Os": "5d",
    "Ir": "5d",
    # Group 3 elements
    "Y": "4d",
    "Sc": "3d",
}


REACTION_TYPES: dict[str, list[str]] = {
    "ligand_coordination": [
        "R4",
        "R5",
        "R16",
        "R17",
        "R18",
        "R20",
        "R29",
        "R30",
        "R31",
        "R37",
        "R57",
        "R58",
        "R59",
        "R60",
        "R61",
    ],
    "oxidative_addition": ["R1", "R19", "R21", "R22", "R38", "R42", "R54"],
    "sigma_bond_metathesis": ["R3", "R15", "R39", "R40", "R41", "R48", "R50", "R51"],
    "ligand_exchange": [
        "R12",
        "R13",
        "R23",
        "R25",
        "R26",
        "R27",
        "R28",
        "R32",
        "R34",
        "R36",
        "R43",
        "R44",
        "R47",
        "R49",
        "R53",
        "R56",
    ],
    "reductive_elimination": ["R2", "R10", "R33", "R52"],
    "migratory_insertion": ["R8", "R9", "R14", "R46"],
    "cycloaddition": ["R6", "R24"],
    "intramolecular_rearrangement": ["R11", "R35"],
    "other": ["R7", "R45", "R55"],
}


# Molecule properties (charge, multiplicity, atom count) from the notebook
MOLECULE_PROPERTIES: dict[str, dict] = {
    "m1": {"charge": 0, "multiplicity": 2, "atoms": 22},
    "m2": {"charge": 0, "multiplicity": 1, "atoms": 17},
    "m3": {"charge": 0, "multiplicity": 2, "atoms": 39},
    "m4": {"charge": 0, "multiplicity": 2, "atoms": 39},
    "m5": {"charge": 0, "multiplicity": 1, "atoms": 32},
    "m6": {"charge": 0, "multiplicity": 1, "atoms": 15},
    "m7": {"charge": 1, "multiplicity": 6, "atoms": 25},
    "m8": {"charge": 0, "multiplicity": 1, "atoms": 6},
    "m9": {"charge": 1, "multiplicity": 6, "atoms": 31},
    "m10": {"charge": 1, "multiplicity": 6, "atoms": 37},
    "m11": {"charge": 1, "multiplicity": 4, "atoms": 37},
    "m12": {"charge": 1, "multiplicity": 4, "atoms": 43},
    "m13": {"charge": 1, "multiplicity": 4, "atoms": 43},
    "m14": {"charge": 1, "multiplicity": 4, "atoms": 43},
    "m15": {"charge": 1, "multiplicity": 6, "atoms": 43},
    "m16": {"charge": 0, "multiplicity": 1, "atoms": 18},
    "m17": {"charge": 0, "multiplicity": 1, "atoms": 13},
    "m18": {"charge": 0, "multiplicity": 1, "atoms": 17},
    "m19": {"charge": 1, "multiplicity": 2, "atoms": 38},
    "m20": {"charge": -1, "multiplicity": 1, "atoms": 2},
    "m21": {"charge": 0, "multiplicity": 1, "atoms": 3},
    "m22": {"charge": 0, "multiplicity": 2, "atoms": 37},
    "m23": {"charge": 0, "multiplicity": 1, "atoms": 6},
    "m24": {"charge": 0, "multiplicity": 4, "atoms": 50},
    "m25": {"charge": 0, "multiplicity": 1, "atoms": 14},
    "m26": {"charge": 0, "multiplicity": 4, "atoms": 64},
    "m27": {"charge": 0, "multiplicity": 1, "atoms": 14},
    "m28": {"charge": 0, "multiplicity": 4, "atoms": 40},
    "m29": {"charge": 0, "multiplicity": 1, "atoms": 38},
    "m30": {"charge": 1, "multiplicity": 2, "atoms": 13},
    "m31": {"charge": -1, "multiplicity": 1, "atoms": 10},
    "m32": {"charge": 0, "multiplicity": 2, "atoms": 23},
    "m33": {"charge": 1, "multiplicity": 2, "atoms": 13},
    "m34": {"charge": 0, "multiplicity": 2, "atoms": 23},
    "m35": {"charge": 1, "multiplicity": 2, "atoms": 13},
    "m36": {"charge": 0, "multiplicity": 2, "atoms": 23},
    "m37": {"charge": 0, "multiplicity": 4, "atoms": 62},
    "m38": {"charge": 0, "multiplicity": 1, "atoms": 2},
    "m39": {"charge": 0, "multiplicity": 3, "atoms": 63},
    "m40": {"charge": 0, "multiplicity": 2, "atoms": 20},
    "m41": {"charge": -1, "multiplicity": 1, "atoms": 1},
    "m42": {"charge": -1, "multiplicity": 2, "atoms": 21},
    "m43": {"charge": 0, "multiplicity": 3, "atoms": 48},
    "m44": {"charge": 0, "multiplicity": 1, "atoms": 14},
    "m45": {"charge": 0, "multiplicity": 2, "atoms": 49},
    "m46": {"charge": 0, "multiplicity": 2, "atoms": 13},
    "m47": {"charge": 0, "multiplicity": 2, "atoms": 61},
    "m48": {"charge": 0, "multiplicity": 1, "atoms": 13},
    "m49": {"charge": 2, "multiplicity": 2, "atoms": 34},
    "m50": {"charge": 0, "multiplicity": 1, "atoms": 17},
    "m51": {"charge": 2, "multiplicity": 2, "atoms": 51},
    "m52": {"charge": 0, "multiplicity": 1, "atoms": 16},
    "m53": {"charge": 2, "multiplicity": 2, "atoms": 67},
    "m54": {"charge": 0, "multiplicity": 4, "atoms": 77},
    "m55": {"charge": -1, "multiplicity": 1, "atoms": 16},
    "m56": {"charge": 0, "multiplicity": 4, "atoms": 81},
    "m57": {"charge": 1, "multiplicity": 1, "atoms": 31},
    "m58": {"charge": 1, "multiplicity": 4, "atoms": 80},
    "m59": {"charge": 0, "multiplicity": 1, "atoms": 15},
    "m60": {"charge": 0, "multiplicity": 2, "atoms": 61},
    "m61": {"charge": 0, "multiplicity": 2, "atoms": 72},
    "m62": {"charge": 0, "multiplicity": 1, "atoms": 2},
    "m63": {"charge": 0, "multiplicity": 2, "atoms": 61},
    "m64": {"charge": 0, "multiplicity": 2, "atoms": 72},
    "m65": {"charge": 0, "multiplicity": 2, "atoms": 51},
    "m66": {"charge": 0, "multiplicity": 2, "atoms": 51},
    "m67": {"charge": 1, "multiplicity": 2, "atoms": 50},
    "m68": {"charge": 0, "multiplicity": 2, "atoms": 13},
    "m69": {"charge": -1, "multiplicity": 2, "atoms": 23},
    "m70": {"charge": 0, "multiplicity": 2, "atoms": 31},
    "m71": {"charge": 0, "multiplicity": 2, "atoms": 32},
    "m72": {"charge": 0, "multiplicity": 2, "atoms": 60},
    "m73": {"charge": 0, "multiplicity": 2, "atoms": 22},
    "m74": {"charge": 0, "multiplicity": 1, "atoms": 38},
    "m75": {"charge": 0, "multiplicity": 2, "atoms": 35},
    "m76": {"charge": 0, "multiplicity": 1, "atoms": 16},
    "m77": {"charge": 0, "multiplicity": 2, "atoms": 38},
    "m78": {"charge": 0, "multiplicity": 2, "atoms": 38},
    "m79": {"charge": 2, "multiplicity": 2, "atoms": 19},
    "m80": {"charge": 0, "multiplicity": 1, "atoms": 4},
    "m81": {"charge": 2, "multiplicity": 2, "atoms": 23},
    "m82": {"charge": 0, "multiplicity": 2, "atoms": 13},
    "m83": {"charge": 0, "multiplicity": 1, "atoms": 34},
    "m84": {"charge": 0, "multiplicity": 2, "atoms": 47},
    "m85": {"charge": 1, "multiplicity": 2, "atoms": 51},
    "m86": {"charge": 0, "multiplicity": 2, "atoms": 4},
    "m87": {"charge": 1, "multiplicity": 3, "atoms": 55},
    "m88": {"charge": 0, "multiplicity": 3, "atoms": 25},
    "m89": {"charge": 0, "multiplicity": 1, "atoms": 26},
    "m90": {"charge": 0, "multiplicity": 3, "atoms": 46},
    "m91": {"charge": 0, "multiplicity": 1, "atoms": 5},
    "m92": {"charge": 0, "multiplicity": 4, "atoms": 39},
    "m93": {"charge": 0, "multiplicity": 1, "atoms": 25},
    "m94": {"charge": 0, "multiplicity": 4, "atoms": 59},
    "m95": {"charge": 0, "multiplicity": 1, "atoms": 25},
    "m96": {"charge": 0, "multiplicity": 4, "atoms": 59},
    "m97": {"charge": 0, "multiplicity": 1, "atoms": 35},
    "m98": {"charge": 0, "multiplicity": 3, "atoms": 60},
    "m99": {"charge": 0, "multiplicity": 1, "atoms": 68},
    "m100": {"charge": -1, "multiplicity": 3, "atoms": 45},
    "m101": {"charge": 0, "multiplicity": 1, "atoms": 2},
    "m102": {"charge": 1, "multiplicity": 1, "atoms": 23},
    "m103": {"charge": 0, "multiplicity": 3, "atoms": 31},
    "m104": {"charge": 0, "multiplicity": 1, "atoms": 43},
    "m105": {"charge": 0, "multiplicity": 1, "atoms": 22},
    "m106": {"charge": -1, "multiplicity": 3, "atoms": 45},
    "m107": {"charge": 0, "multiplicity": 1, "atoms": 16},
    "m108": {"charge": 0, "multiplicity": 3, "atoms": 21},
    "m109": {"charge": 1, "multiplicity": 2, "atoms": 52},
    "m110": {"charge": 1, "multiplicity": 2, "atoms": 86},
    "m111": {"charge": 0, "multiplicity": 2, "atoms": 22},
    "m112": {"charge": 0, "multiplicity": 1, "atoms": 12},
    "m113": {"charge": 0, "multiplicity": 2, "atoms": 34},
    "m114": {"charge": 0, "multiplicity": 3, "atoms": 20},
    "m115": {"charge": -2, "multiplicity": 1, "atoms": 74},
    "m116": {"charge": 0, "multiplicity": 3, "atoms": 76},
    "m117": {"charge": 0, "multiplicity": 2, "atoms": 50},
    "m118": {"charge": 0, "multiplicity": 2, "atoms": 55},
    "m119": {"charge": 0, "multiplicity": 2, "atoms": 68},
    "m120": {"charge": -1, "multiplicity": 1, "atoms": 4},
    "m121": {"charge": 0, "multiplicity": 2, "atoms": 71},
    "m122": {"charge": 0, "multiplicity": 2, "atoms": 66},
    "m123": {"charge": 0, "multiplicity": 1, "atoms": 27},
    "m124": {"charge": 0, "multiplicity": 2, "atoms": 93},
    "m125": {"charge": 0, "multiplicity": 2, "atoms": 55},
    "m126": {"charge": 0, "multiplicity": 2, "atoms": 55},
    "m127": {"charge": 0, "multiplicity": 2, "atoms": 77},
    "m128": {"charge": 0, "multiplicity": 2, "atoms": 80},
    "m129": {"charge": 1, "multiplicity": 2, "atoms": 33},
    "m130": {"charge": 0, "multiplicity": 1, "atoms": 5},
    "m131": {"charge": 1, "multiplicity": 2, "atoms": 38},
    "m132": {"charge": 0, "multiplicity": 2, "atoms": 59},
    "m133": {"charge": 0, "multiplicity": 2, "atoms": 59},
    "m134": {"charge": 0, "multiplicity": 2, "atoms": 33},
    "m135": {"charge": 0, "multiplicity": 1, "atoms": 29},
    "m136": {"charge": 0, "multiplicity": 2, "atoms": 45},
    "m137": {"charge": 0, "multiplicity": 1, "atoms": 17},
    "m138": {"charge": 2, "multiplicity": 2, "atoms": 61},
    "m139": {"charge": 2, "multiplicity": 2, "atoms": 61},
    "m140": {"charge": 2, "multiplicity": 2, "atoms": 61},
    "m141": {"charge": 2, "multiplicity": 2, "atoms": 61},
    "m142": {"charge": 0, "multiplicity": 2, "atoms": 36},
    "m143": {"charge": 0, "multiplicity": 1, "atoms": 31},
    "m144": {"charge": 0, "multiplicity": 2, "atoms": 67},
    "m145": {"charge": 0, "multiplicity": 2, "atoms": 45},
    "m146": {"charge": -1, "multiplicity": 1, "atoms": 22},
    "m147": {"charge": -1, "multiplicity": 2, "atoms": 67},
    "m148": {"charge": 0, "multiplicity": 2, "atoms": 55},
    "m149": {"charge": -1, "multiplicity": 1, "atoms": 27},
    "m150": {"charge": -1, "multiplicity": 2, "atoms": 82},
}


def _atoms_hash(atoms) -> str:
    vec = np.concatenate(
        [atoms.get_positions().ravel(), atoms.get_cell().ravel(), atoms.get_atomic_numbers()]
    )
    return hashlib.md5(vec.tobytes()).hexdigest()


def energy_cached(atoms, calc, cache_dir: Path, cache_tag: str) -> float:
    """Compute potential energy with per-structure caching."""
    # Reset calculator state (mirrors notebook behaviour)
    reset = getattr(calc, "reset", None)
    if callable(reset):
        reset()

    key = _atoms_hash(atoms)
    fpath = cache_dir / f"{cache_tag}_{key}.json"
    if fpath.exists():
        try:
            return json.loads(fpath.read_text())["E"]
        except Exception:  # pragma: no cover - robust fallback
            pass

    atoms_copy = atoms.copy()
    atoms_copy.calc = calc
    e = float(atoms_copy.get_potential_energy())
    fpath.write_text(json.dumps({"E": e}))
    return e


def load_structure(geometries_dir: Path, molecule_id: str):
    """Load structure with charge and multiplicity metadata.

    Follows notebook behavior:
    - prefer `mol.xyz`, fallback to `coord` (TURBOMOLE)
    - charge/multiplicity from MOLECULE_PROPERTIES; fallback to .CHRG and .UHF
    - annotate atoms.info with molecule_id, charge, spin(multiplicity)
    """
    mol_dir = geometries_dir / molecule_id
    mol_xyz = mol_dir / "mol.xyz"
    coord = mol_dir / "coord"

    if mol_xyz.exists():
        atoms = read(str(mol_xyz))
    elif coord.exists():
        atoms = read(str(coord))
    else:
        raise FileNotFoundError(f"No coordinates found for {molecule_id} in {mol_dir}")

    # Defaults
    charge = 0
    multiplicity = 1
    props = MOLECULE_PROPERTIES.get(molecule_id)
    if props:
        charge = int(props.get("charge", 0))
        multiplicity = int(props.get("multiplicity", 1))
    else:
        chrg = mol_dir / ".CHRG"
        if chrg.exists():
            try:
                charge = int(chrg.read_text().strip())
            except Exception:  # pragma: no cover
                pass
        uhf = mol_dir / ".UHF"
        if uhf.exists():
            try:
                unpaired = int(uhf.read_text().strip())
                multiplicity = unpaired + 1
            except Exception:  # pragma: no cover
                pass

    atoms.info["molecule_id"] = molecule_id
    atoms.info["charge"] = charge
    atoms.info["spin"] = multiplicity
    # As in notebook: set per-atom formal charges to zero
    atoms.set_initial_charges([0] * len(atoms))
    return atoms


@dataclass
class Reaction:
    reaction_id: str
    coeffs: list[tuple[float, str]]  # (coeff, molecule)
    metal: str | None = None
    description: str | None = None


def build_reactions_from_res(geometries_dir: Path) -> list[Reaction]:
    res_path = geometries_dir / ".res"
    coeffs_list = parse_res_file(res_path)
    meta = create_reaction_mapping()
    reactions: list[Reaction] = []
    for rid, coeffs in coeffs_list:
        m = meta.get(rid, {})
        reactions.append(
            Reaction(
                reaction_id=rid,
                coeffs=coeffs,
                metal=m.get("metal"),
                description=m.get("description"),
            )
        )
    return reactions


class ROST61Benchmark(zntrack.Node):
    """ROST61 benchmark using notebook-styled definitions and metadata."""

    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    def run(self):
        geometries_dir, cache_dir = get_dirs()
        reactions = build_reactions_from_res(geometries_dir)
        reaction_types = {rid: key for key, rids in REACTION_TYPES.items() for rid in rids}

        calc = self.model.get_calculator()

        out_dir = Path(__file__).parent / "outputs" / self.model_name
        out_dir.mkdir(parents=True, exist_ok=True)

        # Write per-molecule energies for reference/debug
        mol_names = sorted({mol for rxn in reactions for _, mol in rxn.coeffs})
        with open(out_dir / "molecule_energies.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["molecule", "energy_eV"])
            for mol in mol_names:
                e = energy_cached(load_structure(geometries_dir, mol), calc, cache_dir, self.model_name)
                w.writerow([mol, e])

        # Detailed CSV mirroring notebook fields
        with open(out_dir / "rost61_results.csv", "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(
                [
                    "reaction_id",
                    "reaction_type",
                    "metal",
                    "metal_block",
                    "description",
                    "ref_energy",
                    "calc_energy",
                    "error",
                    "abs_error",
                    "max_atoms",
                ]
            )

            for rxn in tqdm(reactions, desc=f"ROST61: {self.model_name}"):
                # Linear combination directly from .res coefficients
                e_comb_eV = 0.0
                max_n = 0
                for coeff, mol in rxn.coeffs:
                    atoms = load_structure(geometries_dir, mol)
                    e_mol = energy_cached(atoms, calc, cache_dir, self.model_name)
                    e_comb_eV += coeff * e_mol
                    props = MOLECULE_PROPERTIES.get(mol, {})
                    max_n = max(max_n, int(props.get("atoms", len(atoms))))

                delta_e_kcal = e_comb_eV * EV_TO_KCAL
                ref = REFERENCE_VALUES[rxn.reaction_id]
                err = delta_e_kcal - ref

                # Determine max atoms among participating structures
                rtype = reaction_types.get(rxn.reaction_id, "unknown")
                mblock = METAL_CLASSIFICATION.get(rxn.metal or "", "unknown")
                
                writer.writerow(
                    [
                        rxn.reaction_id,
                        rtype,
                        rxn.metal or "",
                        mblock,
                        rxn.description or "",
                        ref,
                        delta_e_kcal,
                        err,
                        abs(err),
                        max_n,
                    ]
                )


def build_project(repro: bool = False) -> None:
    project = mlipx.Project()
    nodes: dict[str, ROST61Benchmark] = {}
    for model_name, model in MODELS.items():
        with project.group(model_name):
            nodes[model_name] = ROST61Benchmark(model=model, model_name=model_name)
    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        # Build metadata, then run directly without DVC
        project.build()
        for node in nodes.values():
            node.run()


def test_rost61():
    """Run ROST61 benchmark via pytest."""
    build_project(repro=False)
