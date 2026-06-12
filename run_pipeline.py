"""Run full feasibility analysis: local supplementary + public metadata + optional pathway ORA."""
from __future__ import annotations

import sys
from pathlib import Path

# Allow `python run_pipeline.py` from feasibility_study/
_SRC = Path(__file__).resolve().parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import pandas as pd

from config import OUTPUT_DIR
import figures
import load_moesm
import pathway_context
import public_geo
import statistics_core


def moesm15_to_long(raw: pd.DataFrame) -> pd.DataFrame:
    """Wide MOESM15 layout → long (endpoint × treatment)."""
    rows = []
    for ri in range(len(raw)):
        r = raw.iloc[ri]
        try:
            n = int(r.iloc[0])
        except (ValueError, TypeError):
            continue
        # Locomotion cols 1,2,3 — Coordination 5,6,7 — Treadmill 9,10,11
        triples = [
            ("locomotion_m", 1, 2, 3),
            ("coordination_min", 5, 6, 7),
            ("treadmill_J", 9, 10, 11),
        ]
        for name, c0, c1, c2 in triples:
            try:
                rows.append(
                    {
                        "animal": n,
                        "endpoint": name,
                        "con_veh": float(r.iloc[c0]),
                        "tk_veh": float(r.iloc[c1]),
                        "tk_exo": float(r.iloc[c2]),
                    }
                )
            except (ValueError, TypeError):
                continue
    return pd.DataFrame(rows)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    lines = []

    def log(msg: str):
        lines.append(msg)
        print(msg)

    log("=== Zhang et al. 2017 supplementary integration ===")
    mirna = load_moesm.load_mirna_expression_fig4d()
    mirna_path = OUTPUT_DIR / "mirna_htnsc_astrocyte_summary.csv"
    mirna.to_csv(mirna_path, index=False)
    log(f"Wrote {mirna_path.name} ({len(mirna)} miRNAs)")

    m2 = load_moesm._pick("MOESM2_ESM")
    raw2 = load_moesm.read_first_sheet(m2)
    hdr = raw2.iloc[1].tolist()
    body = raw2.iloc[2:].copy()
    body.columns = [str(c).strip() if c is not None else "" for c in hdr]
    body = body.rename(columns={"Transcript ID(Array Design)": "mirna"})
    body = body.dropna(subset=["mirna"])
    for c in ["htNSC1", "htNSC2", "Astrocyte1", "Astrocyte2"]:
        body[c] = pd.to_numeric(body[c], errors="coerce")
    perm = statistics_core.mirna_specificity_null(
        body["htNSC1"].values,
        body["htNSC2"].values,
        body["Astrocyte1"].values,
        body["Astrocyte2"].values,
        n_perm=8000,
    )
    pd.Series(perm).to_csv(OUTPUT_DIR / "mirna_htnsc_specificity_permutation.csv")
    log(
        "miRNA htNSC vs astrocyte specificity (label-shuffle null, median |d|): "
        f"{perm['median_abs_logfc_htnsc_vs_astro']:.4f}, empirical p={perm['perm_p_median_vs_null']:.4g}"
    )

    figures.plot_mirna_volcano(mirna, OUTPUT_DIR / "fig_mirna_htnsc_vs_astro.png")
    log("Wrote fig_mirna_htnsc_vs_astro.png")

    sox2 = load_moesm.load_sox2_counts_fig1d()
    sox2.to_csv(OUTPUT_DIR / "sox2_counts_3v_arc.csv", index=False)
    t_3v = statistics_core.paired_ttest(sox2["sox2_3v_con"].values, sox2["sox2_3v_tk"].values)
    t_arc = statistics_core.paired_ttest(sox2["sox2_arc_con"].values, sox2["sox2_arc_tk"].values)
    log(f"Sox2+ 3V paired t-test Con vs TK: {t_3v}")
    log(f"Sox2+ ARC paired t-test Con vs TK: {t_arc}")
    figures.plot_sox2_paired(sox2, OUTPUT_DIR / "fig_sox2_paired.png")

    cyto = load_moesm.load_cytokine_exosome_long()
    cyto.to_csv(OUTPUT_DIR / "cytokine_vehicle_exosome_long.csv", index=False)
    cy_results = []
    for g, sub in cyto.groupby("gene"):
        st = statistics_core.paired_ttest(sub["vehicle"].values, sub["exosome"].values)
        st["gene"] = g
        st["cohens_d"] = statistics_core.cohens_d_paired(
            sub["vehicle"].values, sub["exosome"].values
        )
        cy_results.append(st)
    cy_df = pd.DataFrame(cy_results)
    cy_df.to_csv(OUTPUT_DIR / "cytokine_paired_tests.csv", index=False)
    log("Cytokine paired t-tests (vehicle vs exosome):\n" + cy_df.to_string(index=False))
    figures.plot_cytokine_paired(cyto, OUTPUT_DIR / "fig_cytokine_exosome.png")

    raw15 = load_moesm.load_exosome_rescue_phenotypes_fig6c()
    long15 = moesm15_to_long(raw15)
    long15.to_csv(OUTPUT_DIR / "phenotype_exosome_rescue_long.csv", index=False)
    anova_rows = []
    for ep, sub in long15.groupby("endpoint"):
        sc = sub.dropna(subset=["con_veh", "tk_veh", "tk_exo"])
        anova_rows.append(
            {
                "endpoint": ep,
                "n_complete_triples": len(sc),
                **statistics_core.one_way_anova_three_groups(
                    sc["con_veh"].values, sc["tk_veh"].values, sc["tk_exo"].values
                ),
            }
        )
    pd.DataFrame(anova_rows).to_csv(OUTPUT_DIR / "phenotype_three_group_anova.csv", index=False)
    log("MOESM15 one-way ANOVA (3 groups):\n" + pd.DataFrame(anova_rows).to_string(index=False))

    exo_np = load_moesm.load_exosome_nanoparticle_fig4b()
    exo_np.to_csv(OUTPUT_DIR / "exosome_nanoparticle_moesm13.csv", index=False)
    log(f"Wrote exosome_nanoparticle_moesm13.csv ({exo_np.shape[0]} rows)")

    csf = load_moesm.load_csf_mirna_young_fig5a()
    csf.to_csv(OUTPUT_DIR / "csf_mirna_young_moesm14.csv", index=False)
    log(f"Wrote csf_mirna_young_moesm14.csv ({len(csf)} miRNAs)")

    log("\n=== Public GEO context (GSE188646) ===")
    try:
        p = public_geo.save_geo_context()
        log(f"Wrote {p.name} and GSE188646_context.json")
    except Exception as e:
        log(f"GEO step skipped/failed: {e}")

    log("\n=== Pathway ORA (g:Profiler, curated NSC/inflammation query) ===")
    try:
        genes = pathway_context.inflammation_and_nsc_query_genes()
        gdf = pathway_context.gprofiler_ora(genes)
        if len(gdf):
            pathway_context.save_pathway_results(gdf, "gprofiler_nsc_inflammation_ora.csv")
            top = gdf.nsmallest(12, "p_value")
            log("Top terms:\n" + top.to_string(index=False))
        else:
            log("g:Profiler returned no rows (check network).")
    except Exception as e:
        log(f"g:Profiler step failed: {e}")

    summary_path = OUTPUT_DIR / "REPORT_SUMMARY.txt"
    header = (
        "Hypothalamus NSC / exosome feasibility — computational summary\n"
        "Primary data: Zhang et al., Nature 2017 supplementary tables (MOESM).\n"
        "Public anchor: GSE188646 (hypothalamus snRNA young vs aged; metadata only).\n"
        "Interpretation: strong miRNA separation htNSC vs astrocytes; exosome-associated\n"
        "suppression of inflammatory mRNA; Sox2+ decline under ablation — consistent with\n"
        "a modulatable hypothalamic niche relevant to systemic aging biology.\n"
        "\n"
        "SCIENTIFIC LIMITATIONS (read before inferring device efficacy):\n"
        "- MOESM tables are from published mouse experiments; statistics here re-quantify\n"
        "  those measurements but do not replace the peer-reviewed paper.\n"
        "- The miRNA permutation null randomizes replicate labels within each miRNA; it\n"
        "  tests specificity of the htNSC vs astrocyte expression contrast under a\n"
        "  simple exchangeability null, not a genome-wide FDR model.\n"
        "- g:Profiler ORA uses a small curated query (NSC + inflammatory genes) to show\n"
        "  pathway coherence; it is illustrative, not a de novo discovery from RNA-seq.\n"
        "- For miRNA→mRNA→Hallmark + public overlaps, run: python run_extended.py\n"
        "- LIFU feasibility is biological plausibility + safety literature; after running\n"
        "  run_extended.py see outputs/LIFU_evidence_layer.txt (kept separate from MOESM).\n"
        "---\n"
    )
    summary_path.write_text(header + "\n".join(lines), encoding="utf-8")
    log(f"\nWrote {summary_path}")


if __name__ == "__main__":
    main()
