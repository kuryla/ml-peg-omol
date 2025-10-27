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


class Cyclo70_Benchmark(zntrack.Node):
    """
    
    """
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()


    def run(self):
        """Run new benchmark."""

        # Read in data and attach calculator
        
        calc = self.model.get_calculator()

        with open(DATA_PATH / 'dlpno-ccsdt-34.dat', 'r') as lines:
            for i, line in enumerate(lines):
                if i == 0:
                    continue
                items = line.strip().split()
                if len(items) == 0:
                    break
                rxn = items[0]
                
                BH_forward_ref = float(items[1]) * KCAL_TO_EV
                BH_reverse_ref = float(items[2]) * KCAL_TO_EV

                r_labels = [path.stem for path in (DATA_PATH / 'XYZ_CYCLO70' / rxn).glob('r*')]
                ts_labels = [path.stem for path in (DATA_PATH / 'XYZ_CYCLO70' / rxn).glob('TS*')]
                p_labels = [path.stem for path in (DATA_PATH / 'XYZ_CYCLO70' / rxn).glob('p*')]

                BH_forward_model = 0
                BH_reverse_model = 0

                write_dir = OUT_PATH / self.model_name
                write_dir.mkdir(parents=True, exist_ok=True)
    
                for atoms_label in r_labels:
                    atoms = read(DATA_PATH / 'XYZ_CYCLO70' / f"{atoms_label}.xyz")
                    atoms.calc = calc
                    BH_forward_model -= atoms.get_potential_energy()
                    write(write_dir / f"{atoms_label}.xyz", atoms)

                for atoms_label in p_labels:
                    atoms = read(DATA_PATH / 'XYZ_CYCLO70' / f"{atoms_label}.xyz")
                    atoms.calc = calc
                    BH_reverse_model -= atoms.get_potential_energy()
                    write(write_dir / f"{atoms_label}.xyz", atoms)

                for atoms_label in ts_labels:
                    atoms = read(DATA_PATH / 'XYZ_CYCLO70' / f"{atoms_label}.xyz")
                    atoms.calc = calc
                    BH_forward_model += atoms.get_potential_energy()
                    BH_reverse_model += atoms.get_potential_energy()
                    
                    atoms.info['ref_forward_bh'] = BH_forward_ref
                    atoms.info['ref_reverse_bh'] = BH_reverse_ref
                    atoms.info['model_forward_bh'] = BH_forward_model
                    atoms.info['model_reverse_bh'] = BH_reverse_model
                    write(write_dir / f"{atoms_label}.xyz", atoms)


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
            benchmark = Cyclo70_Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_cyclo70_barrier_heights():
    """Run Cyclo70 benchmark via pytest."""
    build_project(repro=True)
