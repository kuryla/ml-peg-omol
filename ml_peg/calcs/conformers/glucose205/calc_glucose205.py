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


class Glucose205_Benchmark(zntrack.Node):
    """

    """
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()


    def get_labels(self):
        self.labels = []
        for system_path in sorted((DATA_PATH / 'Glucose_structures').glob('*.xyz')):
            self.labels.append(system_path.stem)


    def get_ref_energies(self):
        df = pd.read_csv(DATA_PATH / 'glucose.csv')
        self.get_labels()
        self.ref_energies = {}
        for i, label in enumerate(self.labels):
            self.ref_energies[label] = df[' dlpno/cbs(3-4)'][i] * KCAL_TO_EV


    def run(self):
        """Run new benchmark."""
        self.get_ref_energies()
        # Read in data and attach calculator        
        calc = self.model.get_calculator()

        lowest_conf_label = 'alpha_002'

        conf_lowest = read(DATA_PATH / 'Glucose_structures' / f'{lowest_conf_label}.xyz')
        conf_lowest.calc = calc
        E_conf_lowest_model = conf_lowest.get_potential_energy()

        for label, E_ref in self.ref_energies.items():
            # Skip the reference conformer for which the error is automatically zero 
            if label == lowest_conf_label:
                continue
            
            atoms = read(DATA_PATH / 'Glucose_structures' / f'{label}.xyz')
            atoms.calc = calc
            atoms.info['model_rel_energy'] = atoms.get_potential_energy() - E_conf_lowest_model
            atoms.info['ref_energy'] = E_ref
            
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
            benchmark = Glucose205_Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_glucose205():
    """Run Glucose205 benchmark via pytest."""
    build_project(repro=True)
