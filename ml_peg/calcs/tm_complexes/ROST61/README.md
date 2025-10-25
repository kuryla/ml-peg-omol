ROST61 Benchmark (Open‑shell TM Reactions)

- Purpose: Benchmark MLIPs on 61 open‑shell organometallic reactions (adsorption, insertion, ligand exchange, σ‑bond metathesis, etc.).
- Reference level: DLPNO‑CCSD(T1)/CBS adsorption/reaction energies (kcal/mol).
- Data: `data/ROST61/` with `m1..m150/mol.xyz` and `.res` mapping file (self‑contained).
- Outputs: `outputs/<model>/rost61_results.csv` with per‑reaction energies and errors.

Citation
- “Assessing Density Functional Theory for Chemically Relevant Open‑shell Transition Metal Reactions” (Maurer et al., 2021).

Running
- Env: use `torch_env` for MACE OMOL.
- Calculations:
  - `conda run -n torch_env pytest -q ml_peg/calcs/tm_complexes/ROST61/calc_ROST61.py --models mace-omol -s`
- Analysis (optional; ensure Plotly/Kaleido are compatible):
  - `conda run -n torch_env pytest -q ml_peg/analysis/tm_complexes/ROST61/analyse_ROST61.py --models mace-omol -s`

Data Layout
- `data/ROST61/.res` (stoichiometry and IDs), `data/ROST61/m*/mol.xyz`.
- Local cache: `cache_ROST61/` (JSON per‑structure energies) — safe to delete for fresh runs.

Notes
- Stoichiometry is parsed directly from `.res`; charge and multiplicity are taken from the notebook table (.CHRG/.UHF fallback) and passed to calculators where supported.
