"""Tier-2: R fgsea (Hallmark), decoupler PROGENy + DoRothEA consensus on pseudobulk, GenAge overlap, extra GEO metadata."""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pandas as pd

import decoupler_dorothea_pb
import decoupler_progeny_pb
import extra_geo_context
import genage_context
import msigdb_hallmark_gsea
import overlap_fisher
import pseudobulk_pathway_r


def run_fgsea_r(rank_csv: Path, gmt_path: Path, out_csv: Path, project_root: Path, log) -> None:
    r_script = project_root / "r" / "fgsea_hallmark_prerank.R"
    if not r_script.is_file():
        log("fgsea: R script missing; skipped.")
        return
    if not rank_csv.is_file():
        log("fgsea: ranked CSV missing; skipped.")
        return
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    nperm = 2000
    cmd = ["Rscript", str(r_script), str(rank_csv), str(gmt_path), str(out_csv), str(nperm)]
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(project_root),
            capture_output=True,
            text=True,
            timeout=3600,
            encoding="utf-8",
            errors="replace",
        )
        if proc.returncode != 0:
            log(f"fgsea R failed (exit {proc.returncode}).")
            if proc.stderr:
                log(proc.stderr[-2500:])
            return
        log(f"fgsea: wrote {out_csv.name}")
    except FileNotFoundError:
        log("fgsea: Rscript not on PATH; skipped.")
    except subprocess.TimeoutExpired:
        log("fgsea: Rscript timed out; skipped.")


def run_tier2(
    log,
    project_root: Path,
    out_dir: Path,
    union: set[str],
    uni_mmu: set[str],
    ranked_csv: Path,
) -> None:
    log("\n=== Tier-2: R fgsea (MSigDB Hallmark GMT; complements gseapy prerank) ===")
    gmt = msigdb_hallmark_gsea.ensure_hallmark_gmt()
    run_fgsea_r(
        ranked_csv,
        gmt,
        out_dir / "fgsea_hallmark_prerank_results.csv",
        project_root,
        log,
    )

    log("\n=== Tier-2: GenAge model-organism genes (Mus musculus) vs targets ===")
    try:
        genage = genage_context.mouse_genage_symbols()
        gen_in_u = {g.upper() for g in genage} & {g.upper() for g in uni_mmu}
        res = overlap_fisher.fisher_overlap(union, gen_in_u, uni_mmu)
        pd.Series(res).to_csv(out_dir / "fisher_targets_vs_genage_mouse_models.csv")
        pd.Series(sorted(gen_in_u), name="gene").to_csv(
            out_dir / "genage_mouse_models_in_mirtarbase_universe.csv", index=False
        )
        log(str(res))
        log("See data/provenance/GenAge_HAGR.txt for licensing and citation.")
    except Exception as e:
        log(f"GenAge context skipped: {e}")

    log("\n=== Tier-2: decoupler PROGENy (mouse, academic license) on pseudobulk ===")
    try:
        pr = decoupler_progeny_pb.run_progeny_on_pseudobulk(out_dir)
        if pr is None:
            exp = os.environ.get("GSE188646_EXPORT_COUNTS", "").strip().lower() in ("1", "true", "yes")
            if exp:
                log(
                    "Skipped PROGENy: no pseudobulk .mtx (counts are only written after a successful "
                    "r/pseudobulk_edgeR_gse188646.R run). Set GSE188646_RDS or place the canonical RDS "
                    "under data/, ensure pseudobulk completes, then re-run."
                )
            else:
                log(
                    "Skipped: need outputs/gse188646_pseudobulk_counts.mtx and row/col name CSVs "
                    "(set GSE188646_EXPORT_COUNTS=1 when running pseudobulk)."
                )
        else:
            acts, tests = pr
            acts.to_csv(out_dir / "progeny_pseudobulk_pathway_scores.csv", index=False)
            if tests is not None and len(tests):
                tests.to_csv(out_dir / "progeny_pseudobulk_young_vs_aged_welch.csv", index=False)
                log(tests.sort_values("p_two_sided").head(12).to_string(index=False))
            else:
                log("PROGENy activities written; Young/Aged tests skipped (need >=2 reps per group).")
    except Exception as e:
        log(f"PROGENy / decoupler skipped: {e}")

    log("\n=== Tier-2: decoupler DoRothEA (mouse A+B, mlm; exploratory) on pseudobulk ===")
    try:
        dr = decoupler_dorothea_pb.run_dorothea_on_pseudobulk(out_dir)
        if dr is None:
            log(
                "Skipped DoRothEA: need pseudobulk .mtx + metadata (same as PROGENy; export counts when running pseudobulk)."
            )
        else:
            acts_d, tests_d = dr
            acts_d.to_csv(out_dir / "exploratory_dorothea_pseudobulk_tf_scores_mlm.csv", index=False)
            if tests_d is not None and len(tests_d):
                tests_d.to_csv(out_dir / "exploratory_dorothea_pseudobulk_young_vs_aged_mlm.csv", index=False)
                log(tests_d.nsmallest(12, "p_two_sided").to_string(index=False))
            else:
                log("DoRothEA mlm scores written; Young/Aged tests skipped (need >=2 reps per group).")
    except Exception as e:
        log(f"DoRothEA / decoupler skipped: {e}")

    pseudobulk_pathway_r.run_pseudobulk_hallmark_pathway_tests(project_root, out_dir, log)

    log("\n=== Tier-2: extra GEO sample metadata (context only) ===")
    try:
        paths = extra_geo_context.save_extra_geo_sample_tables(out_dir)
        for p in paths:
            log(f"Wrote {p.name}")
        if not paths:
            log("No extra GEO metadata written (GEOparse fetch failed).")
    except Exception as e:
        log(f"Extra GEO metadata skipped: {e}")
