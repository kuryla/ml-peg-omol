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
KJ_TO_EV = units.kJ / units.mol
EV_TO_KJ = 1 / KJ_TO_EV

DATA_PATH = Path(__file__).parent / "data"
OUT_PATH = Path(__file__).parent / "outputs"


class Criegee22_Benchmark(zntrack.Node):
    """
    
    """
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()


    def run(self):
        """Run new benchmark."""

        # Read in data and attach calculator
        
        calc = self.model.get_calculator()

        with open(DATA_PATH / 'reference.txt', 'r') as lines:
            for i, line in enumerate(lines):
                if i == 0:
                    continue
                items = line.strip().split()
                label = items[0]
                BH_ref = float(items[8]) * KJ_TO_EV
                atoms_reac = read(DATA_PATH / 'structures' / f'{label}-reac.xyz')
                atoms_reac.calc = calc
                atoms_reac.info['model_energy'] = atoms_reac.get_potential_energy()
                atoms_reac.info['ref_energy'] = 0

                atoms_ts = read(DATA_PATH / 'structures' / f'{label}-TS.xyz')
                atoms_ts.calc = calc
                atoms_ts.info['model_energy'] = atoms_ts.get_potential_energy()
                atoms_ts.info['ref_energy'] = BH_ref
                
                write_dir = OUT_PATH / self.model_name
                write_dir.mkdir(parents=True, exist_ok=True)
                write(write_dir / f"{label}_rct.xyz", atoms_reac)
                write(write_dir / f"{label}_ts.xyz", atoms_ts)


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
            benchmark = Criegee22_Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_criegee22_barrier_heights():
    """Run Criegee22 benchmark via pytest."""
    build_project(repro=True)
