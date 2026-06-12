"""Run R fgsea (multilevel when available) + limma camera + limma fry on pseudobulk (exploratory)."""
from __future__ import annotations

import subprocess
from pathlib import Path

import msigdb_hallmark_gsea


def run_pseudobulk_hallmark_pathway_tests(
    project_root: Path,
    out_dir: Path,
    log,
) -> None:
    mtx = out_dir / "gse188646_pseudobulk_counts.mtx"
    rn = out_dir / "gse188646_pseudobulk_counts_rownames.csv"
    cn = out_dir / "gse188646_pseudobulk_counts_colnames.csv"
    meta = out_dir / "gse188646_pseudobulk_metadata.csv"
    deg = out_dir / "gse188646_young_vs_aged_deg.csv"
    if not all(p.is_file() for p in (mtx, rn, cn, meta, deg)):
        log(
            "\n=== Exploratory: pseudobulk Hallmark fgsea/camera/fry skipped (need .mtx + metadata + DE) ==="
        )
        return
    gmt = msigdb_hallmark_gsea.ensure_hallmark_gmt()
    log("\n=== Exploratory: fgsea on pseudobulk DE ranks vs Hallmark ===")
    fgsea_out = out_dir / "exploratory_fgsea_hallmark_pseudobulk_de_ranks.csv"
    cmd1 = [
        "Rscript",
        str(project_root / "r" / "fgsea_pseudobulk_de_hallmark.R"),
        str(deg.resolve()),
        str(gmt.resolve()),
        str(fgsea_out.resolve()),
        "2000",
    ]
    try:
        p1 = subprocess.run(
            cmd1,
            cwd=str(project_root),
            capture_output=True,
            text=True,
            timeout=3600,
            encoding="utf-8",
            errors="replace",
        )
        if p1.returncode != 0:
            log(f"fgsea pseudobulk DE failed: {p1.stderr[-2000:]}")
        else:
            log(f"Wrote {fgsea_out.name}")
    except FileNotFoundError:
        log("Rscript not found; fgsea pseudobulk skipped.")
    except subprocess.TimeoutExpired:
        log("fgsea pseudobulk timed out.")

    log("\n=== Exploratory: limma camera on pseudobulk log-CPM vs Hallmark ===")
    cam_out = out_dir / "exploratory_camera_hallmark_pseudobulk_young_vs_aged.csv"
    cmd2 = [
        "Rscript",
        str(project_root / "r" / "camera_pseudobulk_hallmark.R"),
        str(mtx.resolve()),
        str(rn.resolve()),
        str(cn.resolve()),
        str(meta.resolve()),
        str(gmt.resolve()),
        str(cam_out.resolve()),
    ]
    try:
        p2 = subprocess.run(
            cmd2,
            cwd=str(project_root),
            capture_output=True,
            text=True,
            timeout=3600,
            encoding="utf-8",
            errors="replace",
        )
        if p2.returncode != 0:
            log(f"camera pseudobulk failed: {p2.stderr[-2500:]}")
        else:
            log(f"Wrote {cam_out.name}")
    except FileNotFoundError:
        log("Rscript not found; camera skipped.")
    except subprocess.TimeoutExpired:
        log("camera timed out.")

    log("\n=== Exploratory: limma fry on pseudobulk log-CPM vs Hallmark ===")
    fry_out = out_dir / "exploratory_fry_hallmark_pseudobulk_young_vs_aged.csv"
    cmd3 = [
        "Rscript",
        str(project_root / "r" / "fry_pseudobulk_hallmark.R"),
        str(mtx.resolve()),
        str(rn.resolve()),
        str(cn.resolve()),
        str(meta.resolve()),
        str(gmt.resolve()),
        str(fry_out.resolve()),
    ]
    try:
        p3 = subprocess.run(
            cmd3,
            cwd=str(project_root),
            capture_output=True,
            text=True,
            timeout=3600,
            encoding="utf-8",
            errors="replace",
        )
        if p3.returncode != 0:
            log(f"fry pseudobulk failed: {p3.stderr[-2500:]}")
        else:
            log(f"Wrote {fry_out.name}")
    except FileNotFoundError:
        log("Rscript not found; fry skipped.")
    except subprocess.TimeoutExpired:
        log("fry timed out.")
