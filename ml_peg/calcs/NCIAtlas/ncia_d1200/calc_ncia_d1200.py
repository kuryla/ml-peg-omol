import numpy as np
from ase.io import read, write
from ase import units
from ase.calculators.calculator import Calculator
from pathlib import Path
import os
from glob import glob
from mace.calculators import MACECalculator
from mace.calculators import mace_off, mace_omol
from tqdm import tqdm
from fnmatch import fnmatch
import json
from collections import defaultdict
import pandas as pd
import logging
from scipy.stats import kendalltau, pearsonr, spearmanr
from typing import Dict, List, Tuple, Optional
from matplotlib import pyplot as plt
import mlipx
from mlipx.abc import NodeWithCalculator
import zntrack
from copy import deepcopy
from ml_peg.models.get_models import load_models
from ml_peg.models.models import current_models
from ml_peg.calcs.utils.utils import chdir, get_benchmark_data

MODELS = load_models(current_models)

KCAL_TO_EV = units.kcal / units.mol
EV_TO_KCAL = 1 / KCAL_TO_EV

OUT_PATH = Path(__file__).parent / "outputs"
EXCLUDE_NOBLE_GASES = True


class NCIA_D1200_Benchmark(zntrack.Node):
    """
    Benchmarking NCIA_D1200 ionic hydrogen bonds benchmark dataset
    """
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    def get_ref_energies(self, data_path):
        self.ref_energies = {}
        with open(data_path / 'NCIA_D1200_benchmark.txt') as lines:
            for i, line in enumerate(lines):
                if i == 0:
                    continue
                items = line.strip().split()
                label = items[0]
                ref_energy = float(items[1]) * KCAL_TO_EV
                self.ref_energies[label] = ref_energy

    @staticmethod
    def get_monomers(atoms):
        if isinstance(atoms.info['selection_a'], str):
            a_ids = [int(id) for id in atoms.info['selection_a'].split('-')]
            a_ids[0] -= 1
        else:
            a_ids = [int(atoms.info['selection_a'])-1, int(atoms.info['selection_a'])]

        if isinstance(atoms.info['selection_b'], str):
            b_ids = [int(id) for id in atoms.info['selection_b'].split('-')]
            b_ids[0] -= 1
        else:
            b_ids = [int(atoms.info['selection_b'])-1, int(atoms.info['selection_b'])]

        atoms_a = atoms[a_ids[0] : a_ids[1]]
        atoms_b = atoms[b_ids[0] : b_ids[1]]
        assert len(atoms_a) + len(atoms_b) == len(atoms)

        atoms_a.info['charge'] = int(atoms.info['charge_a'])
        atoms_a.info['spin'] = 1

        atoms_b.info['charge'] = int(atoms.info['charge_b'])
        atoms_b.info['spin'] = 1
        return (atoms_a, atoms_b)


    def run(self):
        """Run new benchmark."""

        # Read in data and attach calculator
        data_path = get_benchmark_data("NCIA_D1200.zip") / "NCIA_D1200"
        self.get_ref_energies(data_path)
        
        calc = self.model.get_calculator()

        for label, ref_energy in self.get_ref_energies.items():
            xyz_fname = f"{label}_100.xyz"
            atoms = read(data_path / 'geometries' / xyz_fname)
            # exclude noble gases
            if np.isin(atoms.numbers, [2, 10, 18, 36, 54, 86]).any() and EXCLUDE_NOBLE_GASES:
                continue
            atoms_a, atoms_b = self.get_monomers(atoms)
            atoms.calc = calc
            atoms_a.calc = calc
            atoms_b.calc = calc

            atoms.info['model_int_energy'] = atoms.get_potential_energy() - atoms_a.get_potential_energy() - atoms_b.get_potential_energy()
            atoms.info['ref_int_energy'] = ref_energy
            
            write_dir = OUT_PATH / self.model_name
            write_dir.mkdir(parents=True, exist_ok=True)
            write(write_dir / f"{label}.xyz", atoms)


def build_project(repro: bool = False) -> None:
    """
    Build mlipx project.

    Parameters
    ----------
    repro
        Whether to call dvc repro -f after building.
    """
    project = mlipx.Project()
    benchmark_node_dict = {}

    for model_name, model in MODELS.items():
        with project.group(model_name):
            benchmark = NCIA_D1200_Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_ncia_d1200():
    """Run NCIA_D1200 energies benchmark via pytest."""
    build_project(repro=True)
