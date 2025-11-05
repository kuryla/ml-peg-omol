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

DATA_PATH = Path(__file__).parent / "data"
OUT_PATH = Path(__file__).parent / "outputs"


class MPCONF196_Benchmark(zntrack.Node):
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    @staticmethod
    def get_atoms(atoms_path):
        atoms = read(atoms_path)
        atoms.info['charge'] = 0
        atoms.info['spin'] = 1
        return atoms

    def get_ref_energies(self):
        df = pd.read_excel(DATA_PATH / 'Energies_CCSD(T).xlsx', sheet_name='Old vs New CCSD(T) Reference Va', header=1)
        self.ref_energies = {}
        for row in df.iterrows():
            label = row[1][0]
            if label[-1].isnumeric():
                label = label.replace('_', '')
            E_ref = float(row[1][2]) * KCAL_TO_EV
            self.ref_energies[label] = E_ref

    def run(self):
        """Run new benchmark."""
        MOLECULES = ['FGG', 'GFA', 'GGF', 'WG', 'WGG', 'CAMVES', 'CHPSAR', 'COHVAW', 'GS464992', 'GS557577', 'POXTRD', 'SANGLI', 'YIVNOG']
        self.get_ref_energies()
        # Read in data and attach calculator        
        calc = self.model.get_calculator()

        for molecule in MOLECULES:
            model_abs_energies = []
            ref_abs_energies = []
            current_molecule_labels = []
            for label, E_ref in self.ref_energies.items():
                group_label = label.split('_')[0]
                if label[-1].isnumeric():
                    label = label.replace('_', '')
                if molecule != group_label:
                    continue
                atoms = self.get_atoms(DATA_PATH / f'{label}.xyz')
                atoms.translate(-atoms.get_center_of_mass())
                atoms.calc = calc
                model_abs_energies.append(atoms.get_potential_energy())
                ref_abs_energies.append(E_ref)
                current_molecule_labels.append(label)

            for label, E_model in zip(current_molecule_labels, model_abs_energies):
                atoms = self.get_atoms((DATA_PATH / f'{label}.xyz'))
                atoms.translate(-atoms.get_center_of_mass())
                atoms.info['ref_rel_energy'] = self.ref_energies[label] - np.mean(ref_abs_energies)
                atoms.info['model_rel_energy'] = E_model - np.mean(model_abs_energies)
                
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
            benchmark = MPCONF196_Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_mpconf196():
    """Run MPCONF196 benchmark via pytest."""
    build_project(repro=True)
