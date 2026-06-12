"""
Pathway- and TF-activity–level coupling between the miRNA target burden vector and
pseudobulk-derived pathway / regulon readouts (GSE188646).

Uses the same per-gene weighted burden as journal_tier_crossmodal (from mirna_targets_long
+ MOESM weights). For each PROGENy pathway (and each DoRothEA TF with enough targets):

  - mean_burden_pathway = mean weighted burden over genes in the pathway that appear in the
    merged burden table (intersection with DE gene list used in cross-modal CSV)
  - delta_activity = mean_aged − mean_young from decoupler-derived Young vs Aged tests

Then Spearman across pathways (or TFs) between mean_burden_pathway and delta_activity.

Null: permute weighted_burden across genes within the merged table, recompute pathway means,
repeat (destroys gene→burden assignment while preserving the marginal burden vector).

Exploratory: pathway activities are pseudobulk sample-level projections; burden is miRNA-layer.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def _spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    r, p = stats.spearmanr(x, y, nan_policy="omit")
    if not np.isfinite(r):
        return float("nan"), float("nan")
    return float(r), float(p)


def _pathway_mean_burdens(
    burden_by_gene: dict[str, float],
    net: pd.DataFrame,
    *,
    min_genes: int,
) -> tuple[list[str], np.ndarray]:
    """net columns: source (pathway), target (gene symbol)."""
    names: list[str] = []
    means: list[float] = []
    for pw, sub in net.groupby("source", sort=False):
        genes = {str(t).strip().upper() for t in sub["target"] if str(t).strip()}
        vals = [burden_by_gene[g] for g in genes if g in burden_by_gene]
        if len(vals) < min_genes:
            continue
        names.append(str(pw))
        means.append(float(np.mean(vals)))
    return names, np.asarray(means, dtype=float)


def _tf_mean_burdens(
    burden_by_gene: dict[str, float],
    net: pd.DataFrame,
    *,
    min_genes: int,
) -> tuple[list[str], np.ndarray]:
    """DoRothEA net: source = TF, target = gene."""
    names: list[str] = []
    means: list[float] = []
    for tf, sub in net.groupby("source", sort=False):
        genes = {str(t).strip().upper() for t in sub["target"] if str(t).strip()}
        vals = [burden_by_gene[g] for g in genes if g in burden_by_gene]
        if len(vals) < min_genes:
            continue
        names.append(str(tf))
        means.append(float(np.mean(vals)))
    return names, np.asarray(means, dtype=float)


def run_pathway_tf_coupling(
    out_dir: Path,
    log,
    *,
    min_genes_pathway: int = 8,
    min_genes_tf: int = 10,
    n_perm: int = 800,
    seed: int = 45,
) -> None:
    log("\n=== Pathway / TF–level coupling (burden vs PROGENy & DoRothEA deltas) ===")
    bur_path = out_dir / "exploratory_crossmodal_gene_burden_vs_aging_logfc.csv"
    prog_tests = out_dir / "progeny_pseudobulk_young_vs_aged_welch.csv"
    doro_tests = out_dir / "exploratory_dorothea_pseudobulk_young_vs_aged_mlm.csv"
    if not bur_path.is_file():
        log("Pathway/TF coupling: missing exploratory_crossmodal_gene_burden_vs_aging_logfc.csv; skipped.")
        return
    if not prog_tests.is_file() and not doro_tests.is_file():
        log(
            "Pathway/TF coupling: need progeny_pseudobulk_young_vs_aged_welch.csv and/or "
            "exploratory_dorothea_pseudobulk_young_vs_aged_mlm.csv (Tier-2 / pseudobulk counts); skipped."
        )
        return

    try:
        import decoupler as dc
    except ImportError:
        log("Pathway/TF coupling: decoupler not installed; skipped.")
        return

    bdf = pd.read_csv(bur_path)
    if "gene" not in bdf.columns or "weighted_burden" not in bdf.columns:
        log("Pathway/TF coupling: unexpected burden CSV columns; skipped.")
        return
    bdf["gene_u"] = bdf["gene"].astype(str).str.strip().str.upper()
    burden_by_gene = dict(zip(bdf["gene_u"], bdf["weighted_burden"].astype(float)))
    genes_ordered = bdf["gene_u"].tolist()
    w_vec = bdf["weighted_burden"].astype(float).values.copy()

    rng = np.random.default_rng(seed)
    out: dict = {"n_perm": n_perm, "pathway_block": None, "tf_block": None}

    # --- PROGENy ---
    if prog_tests.is_file():
        ptests = pd.read_csv(prog_tests)
        if not {"pathway", "mean_young", "mean_aged"}.issubset(ptests.columns):
            log("Pathway coupling: unexpected PROGENy tests columns; skipped pathway block.")
        else:
            net_p = dc.op.progeny(organism="mouse", license="academic", verbose=False)
            pnames, pmeans = _pathway_mean_burdens(burden_by_gene, net_p, min_genes=min_genes_pathway)
            if len(pnames) < 4:
                log(f"Pathway coupling: too few pathways after min_genes={min_genes_pathway}; skipped.")
            else:
                delta_map = {}
                for _, r in ptests.iterrows():
                    try:
                        my = float(r["mean_young"])
                        ma = float(r["mean_aged"])
                    except (TypeError, ValueError):
                        continue
                    if np.isfinite(my) and np.isfinite(ma):
                        delta_map[str(r["pathway"])] = ma - my
                deltas = np.array([delta_map.get(n, float("nan")) for n in pnames], dtype=float)
                ok = np.isfinite(deltas) & np.isfinite(pmeans)
                if ok.sum() < 4:
                    log("Pathway coupling: insufficient matched pathways; skipped.")
                else:
                    rho, pval = _spearman(pmeans[ok], deltas[ok])
                    null_r = []
                    for _ in range(n_perm):
                        w_vec[:] = rng.permutation(w_vec)
                        tmp = dict(zip(genes_ordered, w_vec))
                        _, pm = _pathway_mean_burdens(tmp, net_p, min_genes=min_genes_pathway)
                        if len(pm) != len(pnames):
                            continue
                        d2 = deltas.copy()
                        r0, _ = _spearman(pm[ok], d2[ok])
                        if np.isfinite(r0):
                            null_r.append(r0)
                    null_r = np.asarray(null_r, dtype=float)
                    perm_p = float(np.mean(np.abs(null_r) >= abs(rho))) if null_r.size and np.isfinite(rho) else float("nan")
                    out["pathway_block"] = {
                        "n_pathways": int(ok.sum()),
                        "min_genes_per_pathway": min_genes_pathway,
                        "spearman_rho_mean_burden_vs_delta_activity": rho,
                        "spearman_p": pval,
                        "perm_p_rho": perm_p,
                    }
                    pd.DataFrame(
                        {"pathway": [pnames[i] for i in range(len(pnames)) if ok[i]],
                         "mean_weighted_burden": pmeans[ok],
                         "delta_activity_aged_minus_young": deltas[ok]}
                    ).to_csv(out_dir / "exploratory_pathway_burden_vs_progeny_delta.csv", index=False)
                    log(json.dumps(out["pathway_block"], indent=2))

    # --- DoRothEA TF ---
    if doro_tests.is_file():
        dtests = pd.read_csv(doro_tests)
        if not {"tf", "mean_young", "mean_aged"}.issubset(dtests.columns):
            log("TF coupling: unexpected DoRothEA tests columns; skipped TF block.")
        else:
            net_d = dc.op.dorothea(
                organism="mouse",
                levels=["A", "B"],
                license="academic",
                verbose=False,
            )
            tnames, tmeans = _tf_mean_burdens(burden_by_gene, net_d, min_genes=min_genes_tf)
            if len(tnames) < 5:
                log(f"TF coupling: too few TFs after min_genes={min_genes_tf}; skipped.")
            else:
                delta_tf = {}
                for _, r in dtests.iterrows():
                    tf = str(r["tf"]).strip()
                    try:
                        my = float(r["mean_young"])
                        ma = float(r["mean_aged"])
                    except (TypeError, ValueError):
                        continue
                    if np.isfinite(my) and np.isfinite(ma):
                        delta_tf[tf] = ma - my
                dvec = np.array([delta_tf.get(n, float("nan")) for n in tnames], dtype=float)
                ok = np.isfinite(dvec) & np.isfinite(tmeans)
                if ok.sum() < 5:
                    log("TF coupling: insufficient matched TFs; skipped.")
                else:
                    rho2, p2 = _spearman(tmeans[ok], dvec[ok])
                    null2 = []
                    w_vec[:] = bdf["weighted_burden"].astype(float).values
                    for _ in range(n_perm):
                        w_vec[:] = rng.permutation(w_vec)
                        tmp = dict(zip(genes_ordered, w_vec))
                        _, tm = _tf_mean_burdens(tmp, net_d, min_genes=min_genes_tf)
                        if len(tm) != len(tnames):
                            continue
                        d2 = dvec.copy()
                        r0, _ = _spearman(tm[ok], d2[ok])
                        if np.isfinite(r0):
                            null2.append(r0)
                    null2 = np.asarray(null2, dtype=float)
                    perm_p2 = float(np.mean(np.abs(null2) >= abs(rho2))) if null2.size and np.isfinite(rho2) else float("nan")
                    out["tf_block"] = {
                        "n_tfs": int(ok.sum()),
                        "min_targets_per_tf": min_genes_tf,
                        "spearman_rho_mean_burden_vs_delta_tf_activity": rho2,
                        "spearman_p": p2,
                        "perm_p_rho": perm_p2,
                    }
                    pd.DataFrame(
                        {"tf": [tnames[i] for i in range(len(tnames)) if ok[i]],
                         "mean_weighted_burden_on_tf_targets": tmeans[ok],
                         "delta_tf_activity_aged_minus_young": dvec[ok]}
                    ).to_csv(out_dir / "exploratory_tf_burden_vs_dorothea_delta.csv", index=False)
                    log(json.dumps(out["tf_block"], indent=2))

    out["methodology_note"] = (
        "Pathway/TF activities are Welch Young vs Aged on pseudobulk-derived decoupler scores; "
        "burden is MOESM-weighted miRTarBase target burden on genes present in exploratory_crossmodal_* CSV."
    )
    out["caveat"] = "Not a generative joint model; permutation null destroys gene-specific burden while preserving its empirical distribution."

    if out["pathway_block"] is None and out["tf_block"] is None:
        log("Pathway/TF coupling: no blocks written (missing inputs or too few features).")
        return

    (out_dir / "exploratory_pathway_tf_coupling_summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    log("Wrote exploratory_pathway_tf_coupling_summary.json")
