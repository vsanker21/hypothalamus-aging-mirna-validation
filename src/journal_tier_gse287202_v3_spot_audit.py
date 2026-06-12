"""
GSE287202 Visium aging hypothalamus — third-ventricle spot audit.

Downloads GEO supplementary metadata (spatial coordinates + tissue positions) when absent,
audits whether hypothalamus sections contain spots near the ventricular system (3V) despite
the authors excluding VS from formal DE.

Outputs:
  exploratory_gse287202_v3_spot_audit_summary.json
  exploratory_gse287202_v3_spot_audit_per_sample.csv
"""
from __future__ import annotations

import gzip
import io
import json
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests

GEO_SERIES = "GSE287202"
GEO_RAW_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287202/suppl/GSE287202_RAW.tar"
)


def _download_geo_raw(cache_dir: Path, log) -> Path | None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / "GSE287202_RAW.tar"
    if dest.is_file() and dest.stat().st_size > 1_000_000:
        return dest
    log(f"GSE287202 audit: downloading {GEO_RAW_URL} (~463 MB)...")
    try:
        with requests.get(GEO_RAW_URL, stream=True, timeout=600) as r:
            r.raise_for_status()
            with open(dest, "wb") as fh:
                for chunk in r.iter_content(1 << 20):
                    if chunk:
                        fh.write(chunk)
        return dest
    except Exception as exc:
        log(f"GSE287202 download failed: {exc}")
        return None


def _load_visium_positions_from_tar(tar_path: Path) -> dict[str, pd.DataFrame]:
    """Extract tissue_positions CSV per sample from GEO RAW tar (10x Visium)."""
    out: dict[str, pd.DataFrame] = {}
    with tarfile.open(tar_path, "r:*") as tf:
        for member in tf.getmembers():
            name = member.name
            lname = name.lower()
            if "tissue_positions" not in lname:
                continue
            if not (lname.endswith(".csv") or lname.endswith(".csv.gz")):
                continue
            f = tf.extractfile(member)
            if f is None:
                continue
            try:
                raw = f.read()
                if lname.endswith(".gz"):
                    raw = gzip.decompress(raw)
                df = pd.read_csv(io.BytesIO(raw))
            except Exception:
                continue
            base = Path(name).name.replace(".csv.gz", "").replace(".csv", "")
            sample_key = base.replace("_tissue_positions", "")
            out[sample_key] = df
    return out


def _audit_v3_proximity(pos: pd.DataFrame, sample_id: str) -> dict:
    """
    Heuristic: Visium spots with low pixel_y (dorsal) and central pixel_x may lie near 3V
    in coronal hypothalamus sections. Also flag spots in lowest 10% y (ventricular lining proxy).
    """
    df = pos.copy()
    # Standard 10x tissue_positions_list: barcode, in_tissue, array_row, array_col, px_row, px_col
    if df.shape[1] >= 6 and df.columns[0] != "barcode":
        df.columns = ["barcode", "in_tissue", "array_row", "array_col", "px_row", "px_col"][: df.shape[1]]
    px_y = pd.to_numeric(
        df.get(
            "px_row",
            df.get(
                "pxl_row_in_fullres",
                df.get("y", df.get("row", pd.Series(dtype=float))),
            ),
        ),
        errors="coerce",
    )
    px_x = pd.to_numeric(
        df.get(
            "px_col",
            df.get(
                "pxl_col_in_fullres",
                df.get("x", df.get("col", pd.Series(dtype=float))),
            ),
        ),
        errors="coerce",
    )
    in_tissue = df.get("in_tissue", df.get("in_tissue", pd.Series(1, index=df.index)))
    tissue_mask = in_tissue.astype(str).isin(["1", "True", "true", "1.0"]) | (in_tissue == 1)
    valid = tissue_mask & px_y.notna() & px_x.notna()
    n_tissue = int(valid.sum())
    if n_tissue < 10:
        return {"sample_id": sample_id, "n_tissue_spots": n_tissue, "audit_note": "insufficient_spots"}

    y = px_y[valid].astype(float)
    x = px_x[valid].astype(float)
    y_thr = float(np.quantile(y, 0.15))
    x_med = float(np.median(x))
    x_iqr = float(np.subtract(*np.quantile(x, [0.75, 0.25])))
    near_midline = (x >= x_med - x_iqr) & (x <= x_med + x_iqr)
    ventral_band = y <= y_thr
    v3_proxy = near_midline & ventral_band
    return {
        "sample_id": sample_id,
        "n_tissue_spots": n_tissue,
        "n_v3_proxy_spots": int(v3_proxy.sum()),
        "frac_v3_proxy_spots": float(v3_proxy.mean()),
        "px_y_15th_percentile": y_thr,
        "px_x_median": x_med,
        "audit_note": "heuristic_midline_ventral_band_not_anatomical_atlas",
    }


def run_gse287202_v3_spot_audit(out_dir: Path, log) -> None:
    log("\n=== GSE287202 Visium V3-proximity spot audit ===")
    cache = Path(__file__).resolve().parents[1] / "data" / "references" / "geo" / GEO_SERIES
    tar = _download_geo_raw(cache, log)
    if tar is None:
        summ = {
            "geo_series": GEO_SERIES,
            "status": "download_failed",
            "caveat": "Could not retrieve RAW tar; manual audit required.",
        }
        (out_dir / "exploratory_gse287202_v3_spot_audit_summary.json").write_text(
            json.dumps(summ, indent=2), encoding="utf-8"
        )
        log("GSE287202 audit: download failed; wrote failure summary.")
        return

    positions = _load_visium_positions_from_tar(tar)
    if not positions:
        summ = {
            "geo_series": GEO_SERIES,
            "status": "no_tissue_positions_in_tar",
            "tar_members_sample": [],
        }
        with tarfile.open(tar, "r:*") as tf:
            summ["tar_members_sample"] = [m.name for m in tf.getmembers()[:30]]
        (out_dir / "exploratory_gse287202_v3_spot_audit_summary.json").write_text(
            json.dumps(summ, indent=2), encoding="utf-8"
        )
        log("GSE287202 audit: no tissue_positions files found in tar.")
        return

    audit_rows = [_audit_v3_proximity(df, sid) for sid, df in positions.items()]
    audit_df = pd.DataFrame(audit_rows)
    csv_path = out_dir / "exploratory_gse287202_v3_spot_audit_per_sample.csv"
    audit_df.to_csv(csv_path, index=False)

    total_proxy = int(audit_df.get("n_v3_proxy_spots", pd.Series(dtype=int)).sum())
    summ = {
        "geo_series": GEO_SERIES,
        "reference": "Conacher et al. 2025 spatial aging mouse brain (Visium); PRJNA1211604",
        "n_samples_with_positions": len(audit_rows),
        "total_v3_proxy_spots_heuristic": total_proxy,
        "per_sample": audit_rows,
        "paper_caveat": "Authors excluded ventricular system from downstream DE due to insufficient spots.",
        "audit_conclusion": (
            "Heuristic midline-ventral spot proxy suggests some hypothalamus sections may retain "
            "3V-adjacent spots suitable for exploratory (not paper-replicated) niche analysis."
            if total_proxy > 0
            else "Few or no V3-proxy spots detected with pixel-heuristic; consistent with VS exclusion."
        ),
        "method_limitation": (
            "Pixel-coordinate heuristic without Allen CCF registration — anatomical claims require "
            "histology-aligned spot mapping."
        ),
    }
    json_path = out_dir / "exploratory_gse287202_v3_spot_audit_summary.json"
    json_path.write_text(json.dumps(summ, indent=2), encoding="utf-8")
    log(f"Wrote {csv_path.name} ({len(audit_df)} samples); total V3-proxy spots={total_proxy}")
    log(f"Wrote {json_path.name}")
