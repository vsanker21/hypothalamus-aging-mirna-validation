"""
Generate Nature Portfolio-style reporting summary checklist from pipeline outputs.

Outputs:
  outputs/manuscript/NATURE_PORTFOLIO_REPORTING_SUMMARY.md
  outputs/manuscript/NATURE_PORTFOLIO_REPORTING_CHECKLIST.csv
"""
from __future__ import annotations

import argparse
from datetime import date
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _ok(p: Path) -> str:
    return "yes" if p.is_file() else "no"


def build_reporting_summary(out_dir: Path) -> tuple[Path, Path]:
    out_dir = out_dir.resolve()
    ms_dir = out_dir / "manuscript"
    ms_dir.mkdir(parents=True, exist_ok=True)

    checks = [
        ("Data availability — processed outputs", out_dir / "exploratory_crossmodal_mirna_aging_summary.json"),
        ("Data availability — GSE188646 DE", out_dir / "gse188646_young_vs_aged_deg.csv"),
        ("Code availability — pipeline entry", PROJECT_ROOT / "run_extended.py"),
        ("Code availability — requirements", PROJECT_ROOT / "requirements.txt"),
        ("Reproducibility — figure manifest", PROJECT_ROOT / "data" / "provenance" / "FIGURE_PANEL_MANIFEST.csv"),
        ("Reproducibility — SA completeness", out_dir / "SA_COMPLETENESS_CHECK.txt"),
        ("External validation — Jin ABC aging scRNA", out_dir / "exploratory_allen_aging_scrna_validation_summary.json"),
        ("External validation — MERFISH spatial", out_dir / "exploratory_merfish_spatial_validation_summary.json"),
        ("External validation — RSTE3 BIL in situ", out_dir / "exploratory_allen_aging_spatial_validation_summary.json"),
        ("External validation — Allen ISH anatomy", out_dir / "exploratory_allen_ish_marker_anatomy_summary.json"),
        ("External validation — GSE287202 Visium audit", out_dir / "exploratory_gse287202_v3_spot_audit_summary.json"),
        ("Small niche strata DE (48/50/57)", out_dir / "exploratory_gse188646_small_niche_strata_de_summary.json"),
        ("Statistical analysis plan", PROJECT_ROOT / "data" / "provenance" / "STATISTICAL_ANALYSIS_PLAN.txt"),
        ("Zenodo deposit manifest", out_dir / "zenodo" / "ZENODO_DEPOSIT_MANIFEST.csv"),
    ]

    rows = [{"item": item, "path": str(path), "present": _ok(path)} for item, path in checks]
    df = pd.DataFrame(rows)
    csv_path = ms_dir / "NATURE_PORTFOLIO_REPORTING_CHECKLIST.csv"
    df.to_csv(csv_path, index=False)

    n_yes = int((df["present"] == "yes").sum())
    md_lines = [
        "# Nature Portfolio reporting summary (auto-generated checklist)",
        "",
        f"Generated: {date.today().isoformat()}",
        "",
        f"Automated artifact coverage: **{n_yes}/{len(df)}** items present.",
        "",
        "## Computational reproducibility",
        "",
        "- Public reanalysis of GSE188646/GSE87102 with external Allen/BIL spatial validation.",
        "- V3 niche: Jin ABC aging scRNA, MERFISH CCF, RSTE3 in situ Resolve, GSE287202 Visium audit.",
        "",
        "## Ethics",
        "",
        "- Public data only; no new human or animal experiments.",
        "",
        "## Manual items before submission",
        "",
        "- [ ] GitHub repository URL",
        "- [ ] Zenodo DOI",
        "- [ ] Author list, funding, competing interests",
        "",
        "## Checklist",
        "",
    ]
    for _, r in df.iterrows():
        mark = "x" if r["present"] == "yes" else " "
        md_lines.append(f"- [{mark}] {r['item']}")

    md_path = ms_dir / "NATURE_PORTFOLIO_REPORTING_SUMMARY.md"
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    return md_path, csv_path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs-dir", type=Path, default=PROJECT_ROOT / "outputs")
    args = ap.parse_args()
    md, csv = build_reporting_summary(args.outputs_dir)
    print(f"Wrote {md}")
    print(f"Wrote {csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
