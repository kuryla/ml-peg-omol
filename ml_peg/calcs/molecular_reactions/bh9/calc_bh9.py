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


def process_atoms(path):
    with open(path) as lines:
        for i, line in enumerate(lines):
            if i == 1:
                items = line.strip().split()
                charge = int(items[0])
                spin = int(items[1])

    atoms = read(path)
    del atoms.info
    atoms.info['charge'] = charge
    atoms.info['spin'] = spin
    return atoms


def parse_cc_energy(fname):
    with open(fname) as lines:
        for line in lines:
            if 'ref' in line:
                items = line.strip().split()
                return float(items[1]) * KCAL_TO_EV


class BH9_Benchmark(zntrack.Node):
    """
    Benchmarking BH9 reaction benchmark
    """
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    def get_ref_energies(self):
        self.ref_energies = {}
        labels = [path.stem.replace('TS', '') for path in sorted((DATA_PATH / 'BH9_SI' / 'XYZ_FILES').glob('*TS.xyz'))]
        rxn_count = 0
        for label in labels:
            for direction in ['forward', 'reverse']:
                rxn_count += 1
                ref_fname = DATA_PATH / 'BH9_SI' / 'DB_FILES' / 'BH' / f'BH9-BH_{rxn_count}_{direction}.db'
                self.ref_energies[label][direction] = parse_cc_energy(ref_fname)

    def run(self):
        """Run new benchmark."""

        # Read in data and attach calculator
        self.get_ref_energies()
        calc = self.model.get_calculator()

        for fname in sorted((DATA_PATH / 'BH9_SI' / 'XYZ_FILES').glob('*TS.xyz')):
            atoms = process_atoms(fname)
            atoms.calc = calc
            atoms.info['model_energy'] = atoms.get_potential_energy()
            
            # Write both forward and reverse barriers, only forward will be used in analysis here
            label = fname.stem
            if 'TS' in label:
                atoms.info['ref_forward_barrier'] = self.ref_energies[label]['forward']
                atoms.info['ref_reverse_barrier'] = self.ref_energies[label]['reverse']

            
            write_dir = OUT_PATH / self.model_name
            write_dir.mkdir(parents=True, exist_ok=True)
            write(write_dir / f"{fname.stem}.xyz", atoms)


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
            benchmark = BH9_Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_bh9_barrier_heights():
    """Run BH9 barriers benchmark via pytest."""
    build_project(repro=True)
