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


class ACONFL_Benchmark(zntrack.Node):
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    @staticmethod
    def get_atoms(atoms_path):
        atoms = read(atoms_path)
        atoms.info['charge'] = 0
        atoms.info['spin'] = 1
        return atoms
    """
    def get_ref_energies(self):
        df = pd.read_excel(DATA_PATH / 'jp2c02439_si_001.xlsx', sheet_name='ACONFL', skiprows=range(1, 15))
        self.ref_energies = {}

        for idx, row in df.iterrows():
            if idx < 12:
                label = f'12{int(row[1])}'
            elif idx < 29:
                label = f'16{int(row[1])}'
            elif idx < 50:
                label = f'20{int(row[1])}'

            E_ref = float(row[2]) * KCAL_TO_EV
            self.ref_energies[label] = E_ref
    """
    def get_ref_energies(self):
        self.ref_energies = {}
        with open(DATA_PATH / 'ACONFL-main' / '.res') as res_file:
            for line in res_file:
                if '$tmer' in line:
                    items = line.strip().split()
                    label = items[2].replace('/$f', '').replace('A', '')
                    E_ref = float(items[7]) * KCAL_TO_EV
                    self.ref_energies[label] = E_ref

    def run(self):
        """Run new benchmark."""
        self.get_ref_energies()
        # Read in data and attach calculator        
        calc = self.model.get_calculator()

        # Get zero conformer energies
        zero_conf_energies = {12:0, 16:0, 20:0}
        atoms = self.get_atoms(DATA_PATH / 'ACONFL-main' / 'A120/struc.xyz')
        atoms.calc = calc
        zero_conf_energies[12] = atoms.get_potential_energy()

        atoms = self.get_atoms(DATA_PATH / 'ACONFL-main' / 'A160/struc.xyz')
        atoms.calc = calc
        zero_conf_energies[16] = atoms.get_potential_energy()

        atoms = self.get_atoms(DATA_PATH / 'ACONFL-main' / 'A2000/struc.xyz')
        atoms.calc = calc
        zero_conf_energies[20] = atoms.get_potential_energy()

        for label, E_ref in self.ref_energies.items():
            num_carbons = int(label[:2])
            atoms = self.get_atoms(DATA_PATH / 'ACONFL-main' / f'A{label}/struc.xyz')
            atoms.calc = calc
            atoms.info['ref_rel_energy'] = E_ref
            atoms.info['model_rel_energy'] = atoms.get_potential_energy() - zero_conf_energies[num_carbons]

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
            benchmark = ACONFL_Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_aconfl():
    """Run ACONFL benchmark via pytest."""
    build_project(repro=True)
