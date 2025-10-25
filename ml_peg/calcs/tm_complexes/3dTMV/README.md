3dTMV Benchmark (Vertical Ionization Potentials)

- Purpose: Benchmark MLIPs on vertical ionization potentials (IPs) of 28 3d transition‑metal complexes relevant to electrocatalysis.
- Reference level: ph‑AFQMC IPs (kcal/mol).
- Data: `data/3dtmv_structures/1..28/struc.xyz` copied locally for self‑contained runs.
- Outputs: `outputs/<model>/tmv_results.csv` with per‑complex IPs, subset labels (SR, SR/MR, MR), and errors.

Citation
- Neugebauer et al., J. Chem. Theory Comput. 2023, DOI: 10.1021/acs.jctc.3c00617.

Running
- Env: use `torch_env` for MACE OMOL.
- Calculations:
  - `conda run -n torch_env pytest -q ml_peg/calcs/tm_complexes/3dTMV/calc_3dTMV.py --models mace-omol -s`
- Analysis (optional; ensure Plotly/Kaleido versions are compatible):
  - `conda run -n torch_env pytest -q ml_peg/analysis/tm_complexes/3dTMV/analyse_3dTMV.py --models mace-omol -s`

Data Layout
- `data/3dtmv_structures/<complex_id>/struc.xyz`
- Local cache: `cache_3dTMV/` (JSON per‑structure energies) — safe to delete for fresh runs.

Notes
- Charge/multiplicity are taken from the 3dTMV molecular table (per complex) and passed to calculators that support them.
- IP is computed as E(oxidized) − E(neutral) in eV, converted to kcal/mol using the same constants as the notebook.
