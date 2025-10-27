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


class BH2O_36Benchmark(zntrack.Node):
    """
    Benchmarking hydrolysis reaction barriers from 10.1021/acs.jctc.3c00176
    """
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    @staticmethod
    def get_systems(info_path, xyz_dir):
        """ Get names and atoms objects. """
        
        replacements = {
            "amide11": "amide_1_1",
            "amide12": "amide_1_2",
            "amide21": "amide_2_1",
            "amide22": "amide_2_2",
            "basicepoxide1": "basic_epoxide_1",
            "basicepoxide2": "basic_epoxide_2" 
        }
        systems = {}
        
        with open(info_path, 'r') as f:
            data = json.load(f)
            for key in data.keys():
                if 'vacuum' not in key:
                    continue
                items = key.strip().split('_')
                if items[0] not in systems.keys():
                    systems[items[0]] = {}
                systems[items[0]][items[1]] = {}
                systems[items[0]][items[1]]['energy'] = data[key]

                xyz_prefix = items[0]
                if items[0] in replacements:
                    xyz_prefix = replacements[items[0]]
                xyz_prefix += f'_{items[1]}_'
                xyz_path = list(xyz_dir.glob(xyz_prefix+'*'))[0]
                systems[items[0]][items[1]]['xyz_path'] = xyz_path
                systems[items[0]][items[1]]['charge'] = int(str(xyz_path).replace('.xyz', '').split('_')[-1])

        return systems


    def run(self):
        """Run new benchmark."""

        # Read in data and attach calculator
        systems = self.get_systems(DATA_PATH / 'mp2_super.json', DATA_PATH / 'molecules/for_sp')
        
        calc = self.model.get_calculator()

        for identifier, system in systems.items():
            atoms_rct = read(system['rct']['xyz_path'])
            atoms_rct.info['charge'] = int(system['rct']['charge'])
            atoms_rct.info['spin'] = 1
            atoms_rct.info['ref_energy'] = system['rct']['energy']
            atoms_rct.calc = calc
            atoms_rct.info['pred_energy'] = atoms_rct.get_potential_energy()

            atoms_pro = read(system['pro']['xyz_path'])
            atoms_pro.info['charge'] = int(system['pro']['charge'])
            atoms_pro.info['spin'] = 1
            atoms_pro.info['ref_energy'] = system['pro']['energy']
            atoms_pro.calc = calc
            atoms_pro.info['pred_energy'] = atoms_pro.get_potential_energy()

            atoms_ts = read(system['ts']['xyz_path'])
            atoms_ts.info['charge'] = int(system['ts']['charge'])
            atoms_ts.info['spin'] = 1
            atoms_ts.info['ref_energy'] = system['ts']['energy']
            atoms_ts.calc = calc
            atoms_ts.info['pred_energy'] = atoms_ts.get_potential_energy()
            
            write_dir = OUT_PATH / self.model_name
            write_dir.mkdir(parents=True, exist_ok=True)
            write(write_dir / f"{identifier}_rct.xyz", atoms_rct)
            write(write_dir / f"{identifier}_pro.xyz", atoms_pro)
            write(write_dir / f"{identifier}_ts.xyz", atoms_ts)


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
            benchmark = BH2O_36Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_bh2o_36_barrier_heights():
    """Run BH2O-36 hydrolysis benchmark via pytest."""
    build_project(repro=True)
