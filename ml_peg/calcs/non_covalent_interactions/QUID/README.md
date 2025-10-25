QUID Benchmark (Non‑covalent Interactions)

Molecular dimers set QUID contains 42 equilibrium and 128 non-equilibrium structures consisting of a large and a small monomer. They were optimised at PBE0+MBD level and their interaction energies were calculated on different ab initio, DFT, semiempirical, and empirical levels, with further properties also provided at PBE0+MBD level, all stored in the QUID.h5 file. For each dimer label, e.g. F1B1, the chemical dimer formula is available under the key formula, while atomic numbers and their positions for dimer, big_monomer, and small_monomer are stored under the labels atoms and positions. Three further types of information are stored - interaction energies, SAPT components, and other physicochemical properties with the labels Eint, SAPT, and properties. The properties were generated at PBE0 or PBE0 plus MBD level, and van der Waals atomic force components at PBE0 plus MBD, plus D4, and plus XDM levels. The information for the non-equilibrium dimers can be found under the main key dissociation with the number to the dimer name indicating its location along the dissociation curve. The interaction energies and forces are provided in eV and eV/Angstrom, with further details available in Table S3 for the interaction energies and Table S4 for the properties of the Supplimentary Information of the paper Extending quantum-mechanical benchmark accuracy to biological ligand-pocket interactions, pre-print available at ChemRxiv 10.26434/chemrxiv-2025-f6615.

Notes
- Data file: `data/QUID.h5` (included for self‑contained runs).
- Outputs: `outputs/<model>/quid_results.csv` with reference and predicted interaction energies (eV and kcal/mol) and absolute error.
- Reference selection in calc: prioritizes CCSD(T) → CCSD → PBE0+MBD when multiple reference levels are available.

Running (example)
- Calculations (e.g. MACE‑OMOL in `torch_env`):
  - `conda run -n torch_env pytest -q ml_peg/calcs/non_covalent_interactions/QUID/calc_QUID.py --models mace-omol -s`
- Analysis (optional):
  - `conda run -n torch_env pytest -q ml_peg/analysis/non_covalent_interactions/QUID/analyse_QUID.py --models mace-omol -s`
