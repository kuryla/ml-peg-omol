Biochemical Proton‑Transfer Benchmark

- Purpose: Benchmark MLIPs on biochemical proton‑transfer reactions (isolated and microsolvated cases).
- Reference level: as provided in the dataset; barriers/reaction energies in kcal/mol.
- Data: expected under `data/geometries/` with subtrees for isolated/microsolvated reactions.

Citation
- Arantes & Řezáč (2025), “Benchmark of approximate quantum‑chemical and machine‑learning potentials for biochemical proton transfer” (preprint). See the accompanying PDF in the source dataset.

Status
- Data can be copied to `data/geometries/` for self‑contained runs.
- Calculation/analysis modules can be added following the ROST61/3dTMV patterns (read mapping table, set charges/spins per case if applicable, compute barriers and reaction energies).

Data Layout
- `data/geometries/isolated_reactions/<case>/conf*.xyz`
- `data/geometries/microsolvated_reactions/<case>/qcmm/conf*.xyz`
