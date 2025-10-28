import numpy as np
from ase.io import read, write
from ase import Atoms, Atom
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


class RDB7Benchmark(zntrack.Node):
    """
    Benchmarking RDB7 reaction barriers from 10.1038/s41597-022-01529-6
    """
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    @staticmethod
    def get_cc_energy(fname):
        with open(fname, 'r') as lines:
            for line in lines:
                if 'CCSD(T)-F12/cc-pVDZ-F12 energy' in line:
                    energy = float(line.strip().split()[-1]) * units.Hartree
                    return energy

    @staticmethod
    def get_atoms_from_molpro(fname):
        atoms = Atoms(None)
        with open(fname, 'r') as lines:
            read_started = False
            for i, line in enumerate(lines):
                if 'ATOMIC COORDINATES' in line:
                    read_started = True
                    xyz_start = i + 4
                if read_started:
                    if i >= xyz_start:
                        items = line.strip().split()
                        if len(items) == 0:
                            break
                        position = np.array([float(items[3]), float(items[4]), float(items[5])]) * units.Bohr
                        atoms += Atom(symbol=items[1], position=position)
        atoms.info['charge'] = 0
        atoms.info['spin'] = 1
        return atoms

    def run(self):
        """Run new benchmark."""


        # Read in data and attach calculator
        
        calc = self.model.get_calculator()

        for i in range(0, 11961):
            BH_forward_ref = 0
            BH_forward_model = 0
            label = str(i).zfill(6)
            for qm_path in (DATA_PATH / 'ccsd/qm_logs').glob('r*'):
                BH_forward_ref -= self.get_cc_energy(qm_path)
                atoms = self.get_atoms_from_molpro(qm_path)
                atoms.calc = calc
                BH_forward_model -= atoms.get_potential_energy()
            for qm_path in (DATA_PATH / 'ccsd/qm_logs').glob('ts*'):
                BH_forward_model += self.get_cc_energy(qm_path)
                atoms = self.get_atoms_from_molpro(qm_path)
                atoms.calc = calc
                BH_forward_model += atoms.get_potential_energy()

                atoms.info['model_forward_barrier'] = BH_forward_model
                atoms.info['ref_forward_barrier'] = BH_forward_ref

                write_dir = OUT_PATH / self.model_name
                write_dir.mkdir(parents=True, exist_ok=True)
                write(write_dir / f"{label}_ts.xyz", atoms)
            

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
            benchmark = RDB7Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_rdb7_barrier_heights():
    """Run RDB7 benchmark via pytest."""
    build_project(repro=True)
