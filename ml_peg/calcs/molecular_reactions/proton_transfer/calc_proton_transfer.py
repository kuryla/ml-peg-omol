"""Run calculations for the biochemical proton‑transfer benchmark."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from ase.io import read
import mlipx
from mlipx.abc import NodeWithCalculator
import numpy as np
import pandas as pd
import pytest
from tqdm import tqdm
import zntrack

from ml_peg.calcs.utils.utils import chdir
from ml_peg.models.get_models import load_models
from ml_peg.models.models import current_models


EV_TO_KCAL = 1.0 / 0.04336414
MODELS = load_models(current_models)


def get_dirs() -> tuple[Path, Path, Optional[Path]]:
    base = Path(__file__).parent
    geo = base / "data" / "geometries"
    cache = base / "cache_proton"
    cache.mkdir(parents=True, exist_ok=True)
    ods = (base / "data" / "results.ods")
    if not geo.exists():
        raise FileNotFoundError("Proton transfer geometries not found under data/geometries")
    return geo, cache, (ods if ods.exists() else None)


def _atoms_hash(atoms) -> str:
    vec = np.concatenate([atoms.get_positions().ravel(), atoms.get_cell().ravel(), atoms.get_atomic_numbers()])
    return hashlib.md5(vec.tobytes()).hexdigest()


def energy_cached(atoms, calc, cache_dir: Path, cache_tag: str) -> float:
    reset = getattr(calc, "reset", None)
    if callable(reset):
        reset()
    key = _atoms_hash(atoms)
    fpath = cache_dir / f"{cache_tag}_{key}.json"
    if fpath.exists():
        try:
            return json.loads(fpath.read_text())["E"]
        except Exception:
            pass
    cpy = atoms.copy()
    cpy.calc = calc
    e = float(cpy.get_potential_energy())
    fpath.write_text(json.dumps({"E": e}))
    return e


def _load_case_confs(geo_root: Path, subset: str, case: str) -> list[Path]:
    if subset.lower().startswith("iso"):
        p = geo_root / "isolated_reactions" / case
        return sorted(p.glob("conf*.xyz"))
    # microsolvated default
    p = geo_root / "microsolvated_reactions" / case / "qcmm"
    return sorted(p.glob("conf*.xyz"))


def _pick_conf(paths: list[Path], label: Optional[str]) -> Optional[Path]:
    if not paths:
        return None
    if label:
        # accept labels like conf1, 1
        name = label if label.startswith("conf") else f"conf{label}"
        for p in paths:
            if p.stem == name:
                return p
    # fallback: first/middle/last heuristics based on keyword
    if label is None:
        return None
    lab = label.lower()
    if "react" in lab:
        return paths[0]
    if "prod" in lab:
        return paths[-1]
    if "ts" in lab:
        return paths[len(paths)//2]
    return None


def _read_ods(ods_path: Path) -> pd.DataFrame:
    try:
        df = pd.read_excel(ods_path, engine="odf")
        return df
    except Exception:
        return pd.DataFrame()


@dataclass
class PTRow:
    case_id: str
    subset: str
    ref_rxn: Optional[float]
    ref_barrier: Optional[float]
    reactant_label: Optional[str]
    ts_label: Optional[str]
    product_label: Optional[str]


def _rows_from_ods(ods: Optional[Path]) -> list[PTRow]:
    rows: list[PTRow] = []
    if ods is None:
        return rows
    df = _read_ods(ods)
    if df.empty:
        return rows
    # normalize columns
    cols = {c.lower(): c for c in df.columns}
    for _, rec in df.iterrows():
        cid = str(rec.get(cols.get("id") or cols.get("case")))
        subset = str(rec.get(cols.get("subset"), "isolated"))
        ref_rxn = rec.get(cols.get("e_ref_rxn_kcal")) if cols.get("e_ref_rxn_kcal") in rec else None
        ref_bar = rec.get(cols.get("e_ref_barrier_kcal")) if cols.get("e_ref_barrier_kcal") in rec else None
        react = rec.get(cols.get("reactant_conf"))
        ts = rec.get(cols.get("ts_conf"))
        prod = rec.get(cols.get("product_conf"))
        rows.append(PTRow(cid, subset, ref_rxn if pd.notna(ref_rxn) else None, ref_bar if pd.notna(ref_bar) else None,
                          str(react) if pd.notna(react) else None,
                          str(ts) if pd.notna(ts) else None,
                          str(prod) if pd.notna(prod) else None))
    return rows


class ProtonTransferBenchmark(zntrack.Node):
    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    def run(self):
        geo_root, cache_dir, ods = get_dirs()
        calc = self.model.get_calculator()
        rows = _rows_from_ods(ods)

        out_dir = Path(__file__).parent / "outputs" / self.model_name
        out_dir.mkdir(parents=True, exist_ok=True)

        with open(out_dir / "pt_results.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["case_id", "subset", "ref_rxn", "calc_rxn", "ref_barrier", "calc_barrier", "abs_err_rxn", "abs_err_barrier"])

            # If ODS missing, fall back to scanning cases under isolated_reactions
            if not rows:
                # construct minimal rows from folder names
                for subset in ("isolated_reactions", "microsolvated_reactions"):
                    base = geo_root / subset
                    if not base.exists():
                        continue
                    for case_dir in sorted([p for p in base.iterdir() if p.is_dir()]):
                        cid = case_dir.name
                        rows.append(PTRow(cid, "isolated" if subset.startswith("iso") else "microsolvated", None, None, "reactant", "ts", "product"))

            for row in tqdm(rows, desc=f"PT: {self.model_name}"):
                paths = _load_case_confs(geo_root, row.subset, row.case_id)
                if not paths:
                    continue
                p_react = _pick_conf(paths, row.reactant_label) or (paths[0])
                p_prod = _pick_conf(paths, row.product_label) or (paths[-1])
                p_ts = _pick_conf(paths, row.ts_label) or (paths[len(paths)//2])

                e_r = energy_cached(read(str(p_react)), calc, cache_dir, self.model_name)
                e_p = energy_cached(read(str(p_prod)), calc, cache_dir, self.model_name)
                e_ts = energy_cached(read(str(p_ts)), calc, cache_dir, self.model_name)

                rxn = (e_p - e_r) * EV_TO_KCAL
                barrier = (e_ts - e_r) * EV_TO_KCAL

                abs_rxn = abs(rxn - row.ref_rxn) if row.ref_rxn is not None else float("nan")
                abs_bar = abs(barrier - row.ref_barrier) if row.ref_barrier is not None else float("nan")
                w.writerow([row.case_id, row.subset, row.ref_rxn, rxn, row.ref_barrier, barrier, abs_rxn, abs_bar])


def build_project(repro: bool = False) -> None:
    project = mlipx.Project()
    nodes: dict[str, ProtonTransferBenchmark] = {}
    for model_name, model in MODELS.items():
        with project.group(model_name):
            nodes[model_name] = ProtonTransferBenchmark(model=model, model_name=model_name)
    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()
        for n in nodes.values():
            n.run()


def test_proton_transfer():
    build_project(repro=False)

