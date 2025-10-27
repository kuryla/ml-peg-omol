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


class DIPCONFS_Benchmark(zntrack.Node):
    """
    Benchmarking DipCONFS conformers from 10.1021/acs.jctc.3c00176
    """
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()


    def run(self):
        """Run new benchmark."""

        # Read in data and attach calculator        
        calc = self.model.get_calculator()

        df = pd.read_excel(DATA_PATH / 'ct4c00801_si_004.xlsx', sheet_name='Conformational Energies in kcal')

        for zero_conf_label, label, E_rel_ref in zip(
                df['Reference Conformer'].tolist(),
                df['Conformer'].tolist(),
                df["PNO-LCCSD(T)-F12b/AVQZ’"].tolist()
        ):
            # Get reference energy
            E_rel_ref = float(E_rel_ref) * KCAL_TO_EV
            zero_conf_label = zero_conf_label.replace('/', '-')
            label = label.replace('/', '-')

            # Get zero ref conformer model energy
            zero_conf = read(DATA_PATH / zero_conf_label / 'struc.xyz')
            zero_conf.calc = calc
            E_model_zero_conf = zero_conf.get_potential_energy()
            
            # Get current conformer model energy
            atoms = read(DATA_PATH / label / 'struc.xyz')
            atoms.calc = calc
            atoms.info['model_rel_energy'] = atoms.get_potential_energy() - E_model_zero_conf
            atoms.info['ref_energy'] = E_rel_ref
            
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
            benchmark = DIPCONFS_Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_dipconfs():
    """Run DipCONFS benchmark via pytest."""
    build_project(repro=True)
