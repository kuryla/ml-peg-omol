"""
Compute the L14 complexes binding energy benchmark.

J. Chem. Phys. 161, 234103 (2024).
"""

from __future__ import annotations

from pathlib import Path

from ase import units
from ase.io import read, write
import mlipx
from mlipx.abc import NodeWithCalculator
from tqdm import tqdm
import zntrack

from ml_peg.calcs.utils.utils import chdir
from ml_peg.models.get_models import load_models
from ml_peg.models.models import current_models

MODELS = load_models(current_models)

KCAL_TO_EV = units.kcal / units.mol
EV_TO_KCAL = 1 / KCAL_TO_EV

OUT_PATH = Path(__file__).parent / "outputs"


class L14Benchmark(zntrack.Node):
    """Benchmark the L14 dataset."""

    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    """
    Starting indices of monomers.
    Some complexes contain more than two monomers.
    Starting atom index of each monomer is given below.
    """
    MONOMER_START_ATOMS = {
        "2a": [0, 72],
        "2b": [0, 72],
        "3b": [0, 98],
        "C2C2PD": [0, 36],
        "C3A": [0, 15],
        "C3GC": [0, 29],
        "CBH": [0, 56],
        "Da1": [0, 43],
        "Da2": [0, 44],
        "GCGC": [0, 29],
        "GGG": [0, 16, 32],
        "PHE": [0, 29, 58],
        "S8-1": [0, 68],
        "S8-2": [0, 68],
    }

    def get_ref_energies(self, data_path):
        """
        Get reference energies.

        Parameters
        ----------
        data_path
            Path to data.
        """
        self.ref_energies = {}
        for atoms_path in data_path.glob("*.xyz"):
            label = atoms_path.stem
            atoms = read(atoms_path)
            ref_energy = atoms.info["CCSD(T)/CBS"] * KCAL_TO_EV
            self.ref_energies[label] = ref_energy

    def get_monomers(self, atoms, label):
        """
        Get ASE atoms objects of the monomers.

        Parameters
        ----------
        atoms
            ASE atoms object of the structure.
        label
            String identifier of the system.

        Returns
        -------
        list[ASE.Atoms]
            A list containing ASE atoms of the monomers.
        """
        monomers_list = []
        monomer_indices = self.MONOMER_START_ATOMS[label] + [len(atoms)]
        for i in range(len(monomer_indices) - 1):
            start_id = monomer_indices[i]
            end_id = monomer_indices[i + 1]
            monomer = atoms[start_id:end_id]
            monomer.info["charge"] = 0
            monomer.info["spin"] = 1
            monomers_list.append(monomer)
        return monomers_list

    def run(self):
        """Run new benchmark."""
        # Read in data and attach calculator
        data_path = Path("/home/dk584/.cache/ml_peg/L14")
        """
        data_path = (
            download_s3_data(
                filename="L14.zip",
                key="inputs/supramolecular/L14.zip",
            )
            / "L14"
        )
        """
        self.get_ref_energies(data_path)
        calc = self.model.get_calculator()

        # Add D3 calculator for this test
        calc = self.model.add_d3_calculator(calc)

        for label, ref_energy in tqdm(self.ref_energies.items()):
            xyz_fname = f"{label}.xyz"
            atoms = read(data_path / xyz_fname)
            atoms.info["charge"] = 0
            atoms.info["spin"] = 1
            atoms.calc = calc
            atoms.translate(-atoms.get_center_of_mass())
            model_int_energy = atoms.get_potential_energy()

            monomers_list = self.get_monomers(atoms, label)
            for monomer in monomers_list:
                monomer.calc = calc
                monomer.translate(-monomer.get_center_of_mass())
                model_int_energy -= monomer.get_potential_energy()

            atoms.info["model_int_energy"] = model_int_energy
            atoms.info["ref_int_energy"] = ref_energy
            atoms.calc = None

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
            benchmark = L14Benchmark(
                model=model,
                model_name=model_name,
            )
            benchmark_node_dict[model_name] = benchmark

    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_l14():
    """Run L14 benchmark via pytest."""
    build_project(repro=True)
