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
from rdkit import Chem
from ase import Atoms


MODELS = load_models(current_models)

KCAL_TO_EV = units.kcal / units.mol
EV_TO_KCAL = 1 / KCAL_TO_EV

DATA_PATH = Path(__file__).parent / "data"
OUT_PATH = Path(__file__).parent / "outputs"


class OpenFF_Tors_Benchmark(zntrack.Node):
    """

    """
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    def run(self):
        # Read in data and attach calculator 
        calc = self.model.get_calculator()
        with open(DATA_PATH / 'MP2_heavy-aug-cc-pVTZ_torsiondrive_data.json', 'r') as file:
            data = json.load(file)

        for molecule_id, conf in data.items():
            charge = conf['metadata']['mol_charge']
            spin = conf['metadata']['mol_multiplicity']
            smiles = conf['metadata']['mapped_smiles']
            mol = Chem.MolFromSmiles(smiles)
            symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
            atom_map = {atom.GetIntProp("molAtomMapNumber"): idx for idx, atom in enumerate(mol.GetAtoms()) if atom.HasProp("molAtomMapNumber")}
            remapped_symbols = [symbols[atom_map[i]] for i in range(len(symbols))]

            for i, (ref_energy, positions) in enumerate(zip(conf['metadata']['final_energies'], conf['metadata']['final_geometries'])):
                label = f"{molecule_id}_{i}"
                atoms = Atoms(symbols=remapped_symbols, positions=np.array(positions)*units.Bohr)
                atoms.info['charge'] = charge
                atoms.info['spin'] = spin
                atoms.calc = calc

                if i == 0:
                    E_ref_zero_conf = ref_energy * units.Hartree
                    E_model_zero_conf = atoms.get_potential_energy()
                else:
                    atoms.info['ref_rel_energy'] = ref_energy * units.Hartree - E_ref_zero_conf
                    atoms.info['model_rel_energy'] = atoms.get_potential_energy() - E_model_zero_conf
                    
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
            benchmark = OpenFF_Tors_Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_openff_tors():
    """Run OpenFF-Tors benchmark via pytest."""
    build_project(repro=True)
