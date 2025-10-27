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


class TMCONF_Benchmark(zntrack.Node):
    """
    Benchmarking TMCONF40 CC subset ionic hydrogen bonds benchmark dataset
    """
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    def get_ref_energies(self, data_path):
        self.ref_energies = {}
        with open(data_path / 'TMCONF5_Econf.csv', 'r') as lines:
            conf_id = ''
            molecule = ''

            for line in lines:
                if ',,' in line:
                    conf_id = ''
                    molecule = ''

                elif 'Econf' in line:
                    molecule = line.strip().split(',')[0].lower()
                    zero_conf = read(data_path / f'{molecule}_1-cc.xyz')

                elif 'CC' in line:
                    items = line.strip().split(',')
                    conf_id = items[0].split('-')[1]
                    # Get CC reference
                    E_rel_ref = float(items[1]) * KCAL_TO_EV
                    label = f"{molecule}_{conf_id}"
                    self.ref_energies[label] = E_rel_ref

    def run(self):
        """Run new benchmark."""

        # Read in data and attach calculator
        self.get_ref_energies(DATA_PATH / 'TMCONF5')
        
        calc = self.model.get_calculator()

        for label, ref_energy in self.get_ref_energies.items():
            xyz_path = DATA_PATH / f"{label}_cc.xyz"
            atoms = read(xyz_path)
            atoms.calc = calc

            atoms.info['model_energy'] = atoms.get_potential_energy()
            atoms.info['ref_energy'] = ref_energy
            
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
            benchmark = TMCONF_Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_tmconf():
    """Run TMCONF conformation energies benchmark via pytest."""
    build_project(repro=True)
