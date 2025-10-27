import numpy as np
from ase.io import read, write
from ase import units
from ase.calculators.calculator import Calculator
from pathlib import Path
import os
import ase
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


class PLA15_Benchmark(zntrack.Node):
    """
    Benchmarking PLA15 dataset
    """
    
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    # ------------------------------------------------------------
    # PDB processing functions
    # ------------------------------------------------------------
    def extract_charge_and_selections(pdb_path: Path) -> Tuple[float, float, float, str, str]:
        """Extract charge and selection information from PDB REMARK lines"""
        total_charge = qa = qb = 0.0
        selection_a = selection_b = ""
        
        with open(pdb_path, 'r') as f:
            for line in f:
                if not line.startswith('REMARK'):
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        break
                    continue
                
                parts = line.split()
                if len(parts) < 3:
                    continue
                    
                tag = parts[1].lower()
                
                if tag == 'charge':
                    total_charge = float(parts[2])
                elif tag == 'charge_a':
                    qa = float(parts[2])
                elif tag == 'charge_b': 
                    qb = float(parts[2])
                elif tag == 'selection_a':
                    selection_a = ' '.join(parts[2:])
                elif tag == 'selection_b':
                    selection_b = ' '.join(parts[2:])
        
        return total_charge, qa, qb, selection_a, selection_b

    def separate_protein_ligand_simple(pdb_path: Path):
        """Simple separation based on residue names"""
        import MDAnalysis as mda
        
        # Load with MDAnalysis
        u = mda.Universe(str(pdb_path))
        
        # Simple separation: ligand = UNK residues, protein = everything else
        protein_atoms = []
        ligand_atoms = []
        
        for atom in u.atoms:
            if atom.resname.strip().upper() in ['UNK', 'LIG', 'MOL']:
                ligand_atoms.append(atom)
            else:
                protein_atoms.append(atom)
        
        return u.atoms, protein_atoms, ligand_atoms

    def mda_atoms_to_ase(atom_list, charge: float, identifier: str) -> ase.Atoms:
        """Convert MDAnalysis atoms to ASE Atoms object"""
        from ase import Atoms
        
        if not atom_list:
            atoms = Atoms()
            atoms.info.update({'charge': charge, 'identifier': identifier})
            return atoms
        
        symbols = []
        positions = []
        
        for atom in atom_list:
            # Get element symbol
            try:
                elem = (atom.element or "").strip().title()
            except:
                elem = ""
            
            if not elem:
                # Fallback: first letter of atom name
                elem = "".join([c for c in atom.name if c.isalpha()])[:1].title() or "C"
            
            symbols.append(elem)
            positions.append(atom.position)
        
        atoms = Atoms(symbols=symbols, positions=np.array(positions))
        atoms.info.update({'charge': charge, 'identifier': identifier})
        return atoms

    def process_pdb_file(self, pdb_path: Path) -> Dict[str, ase.Atoms]:
        """Process one PDB file and return complex + separated fragments"""
        
        total_charge, charge_a, charge_b, _, _ = self.extract_charge_and_selections(pdb_path)
        
        try:
            all_atoms, protein_atoms, ligand_atoms = self.separate_protein_ligand_simple(pdb_path)
            
            if len(ligand_atoms) == 0:
                logging.warning(f"No ligand atoms found in {pdb_path.name}")
                return {}
            
            if len(protein_atoms) == 0:
                logging.warning(f"No protein atoms found in {pdb_path.name}")
                return {}
            
            base_id = pdb_path.stem
            
            # Create ASE objects
            complex_atoms = self.mda_atoms_to_ase(list(all_atoms), total_charge, base_id)
            protein_frag = self.mda_atoms_to_ase(protein_atoms, charge_a, base_id) 
            ligand = self.mda_atoms_to_ase(ligand_atoms, charge_b, base_id)
            
            return {
                'complex': complex_atoms,
                'protein': protein_frag,
                'ligand': ligand
            }
            
        except Exception as e:
            logging.warning(f"Error processing {pdb_path}: {e}")
            return {}

    # ------------------------------------------------------------
    # Reference energy parsing
    # ------------------------------------------------------------
    
    def parse_pla15_references(path: Path) -> Dict[str, float]:
        """Parse PLA15 reference total energies (kcal/mol -> eV)"""
        ref: Dict[str, float] = {}
        
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line or line.lower().startswith("no.") or line.startswith("-"):
                continue
                
            parts = line.split()
            if len(parts) < 3:
                continue
                
            try:
                energy_kcal = float(parts[-1])
            except ValueError:
                continue
                
            # Extract full identifier with residue type
            full_identifier = parts[1].replace(".pdb", "")
            
            # Extract base identifier by removing residue type suffix
            # Format: "1ABC_15_lys" -> "1ABC_15"
            identifier_parts = full_identifier.split('_')
            if len(identifier_parts) >= 3:
                # Assume last part is residue type (lys, arg, asp, etc.)
                base_identifier = '_'.join(identifier_parts[:-1])
            else:
                # Fallback: use full identifier if format is unexpected
                base_identifier = full_identifier
                
            energy_eV = energy_kcal * KCAL_TO_EV  # Convert to eV
            ref[base_identifier] = energy_eV
            
        return ref


    def run(self):
        """Run new benchmark."""

        # Read in data and attach calculator
        ref_energies = self.parse_plf547_references(DATA_PATH / 'reference_energies.txt')
        
        calc = self.model.get_calculator()

        for label, ref_energy in ref_energies.items():
            pdb_fname = DATA_PATH / f'{label}.pdb'
            
            fragments = self.process_pdb_file(pdb_fname)
            if not fragments:
                continue

            complex_atoms = fragments['complex']
            complex_atoms.calc = calc
            complex_energy = complex_atoms.get_potential_energy()

            protein_atoms = fragments['protein']
            protein_atoms.calc = calc
            protein_energy = protein_atoms.get_potential_energy()

            ligand_atoms = fragments['ligand']
            ligand_atoms.calc = calc
            ligand_energy = ligand_atoms.get_potential_energy()

            complex_atoms.info['model_int_energy'] = complex_energy - protein_energy - ligand_energy
            complex_atoms.info['ref_int_energy'] = ref_energy
            
            write_dir = OUT_PATH / self.model_name
            write_dir.mkdir(parents=True, exist_ok=True)
            write(write_dir / f"{label}_complex.xyz", complex_atoms)


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
            benchmark = PLA15_Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_pla15():
    """Run PLA15 benchmark via pytest."""
    build_project(repro=True)
