MME55 Metalloenzyme Benchmark

- Purpose: Benchmark MLIPs on metalloenzyme reaction steps (reaction energies and barrier heights) across 55 representative cases.
- Reference level: as defined in the accompanying notebook/results; energies in kcal/mol.
- Data: expected under `data/MME55_structures/` (enzyme subfolders with structures).

Citation
- Please insert the primary literature reference for the MME55 set used here (the local notebook does not include a formal citation). If you have a DOI or BibTeX entry, add it below.

Status
- Data copied locally for self‑contained runs.
- Calculation/analysis modules can be added following the ROST61/3dTMV templates; charge/multiplicity handling is supported similarly.

Data Layout
- `data/MME55_structures/<enzyme>/<structure>.xyz` (from your dataset).

Running
- To integrate fully, implement `calc_MME55.py` and `analyse_MME55.py` analogous to ROST61/3dTMV, then run via `ml_peg calc`/`ml_peg analyse`.
