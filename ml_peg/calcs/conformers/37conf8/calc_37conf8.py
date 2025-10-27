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


class Benchmark_37CONF8(zntrack.Node):
    """
    Benchmarking 37CONF8 dataset
    """
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    def run(self):
        """Run new benchmark."""

        df = pd.read_excel(DATA_PATH / '37Conf8_data.xlsx', sheet_name='Rel_Energy_SP', header=2)
        calc = self.model.get_calculator()

        write_dir = OUT_PATH / self.model_name
        write_dir.mkdir(parents=True, exist_ok=True)

        for i in range(0, len(df)-3):
            molecule_name = df.iloc[i][0].strip()
            conf_id = int(df.iloc[i][1])
            label = f"{molecule_name}_{conf_id}"
            if conf_id == 1:
                zero_conf = read(DATA_PATH / 'PBEPBE-D3' / f'{label}_PBEPBE-D3.xyz')
                zero_conf.info['charge'] = 0
                zero_conf.info['spin'] = 1
                zero_conf.calc = calc
                E_model_zero_conf = zero_conf.get_potential_energy()
            else:      
                atoms = read(DATA_PATH / 'PBEPBE-D3' / f'{label}_PBEPBE-D3.xyz')
                atoms.info['charge'] = 0
                atoms.info['spin'] = 1
                atoms.calc = calc
                atoms.info['model_rel_energy'] = atoms.get_potential_energy() - E_model_zero_conf
                atoms.info['ref_energy'] = float(df.iloc[i][2]) * KCAL_TO_EV
                write(write_dir/ f"{label}.xyz", atoms)
            

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
            benchmark = Benchmark_37CONF8(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_37conf8_conformer_energies():
    """Run 37CONF8 benchmark via pytest."""
    build_project(repro=True)
