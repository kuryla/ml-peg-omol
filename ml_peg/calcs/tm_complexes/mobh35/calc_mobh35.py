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


class MOBH35_Benchmark(zntrack.Node):
    """
    Benchmarking MOBH35 reaction benchmark
    """
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    def get_ref_energies(self):
        self.ref_barriers = {}
        with open(DATA_PATH / 'reference.txt') as lines:
            for i, line in enumerate(lines):
                if i < 2:
                    continue
                items = line.strip().split()
                label = int(items[0])
            
                self.ref_barriers[label]['forward'] = float(items[4]) * KCAL_TO_EV
                self.ref_barriers[label]['reverse'] = float(items[5]) * KCAL_TO_EV

    def run(self):
        """Run new benchmark."""

        # Read in data and attach calculator
        self.get_ref_energies()
        
        calc = self.model.get_calculator()

        for fname in (DATA_PATH / 'structures').glob('*'):
            atoms = read(fname)
            if '+' in fname.stem:
                atoms.info['charge'] = 1
            else:
                atoms.info['charge'] = 0

            atoms.calc = calc
            atoms.info['model_energy'] = atoms.get_potential_energy()

            # Write ref barriers on the TS atoms
            if 'ts' in fname.stem:
                label = int(fname.stem.replace('ts', '').replace('+', ''))
                atoms.info['ref_forward_barrier'] = self.ref_barriers[label]['forward']
                atoms.info['ref_reverse_barrier'] = self.ref_barriers[label]['reverse']
            
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
            benchmark = MOBH35_Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_mobh35_barrier_heights():
    """Run MOBH35 barriers benchmark via pytest."""
    build_project(repro=True)
