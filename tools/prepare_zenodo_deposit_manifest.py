"""
Prepare Zenodo deposit manifest (file list + sizes) for large artifacts excluded from Git.

Outputs:
  outputs/zenodo/ZENODO_DEPOSIT_MANIFEST.csv
  outputs/zenodo/ZENODO_DEPOSIT_README.txt
"""
from __future__ import annotations

import argparse
import hashlib
from datetime import date
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]

ZENODO_INCLUDE = [
    ("data/GSE188646_hypo.integrated.final.20210719.RDS", "GSE188646 Seurat object (publisher)"),
    ("data/references/bil_rstE3/RSTE3_metadata.csv", "Jin RSTE3 BIL metadata"),
    ("data/references/bil_rstE3/RSTE3_cellxgene.csv", "Jin RSTE3 Resolve expression"),
    ("outputs/gse188646_young_vs_aged_deg.csv", "Primary pseudobulk DE"),
    ("outputs/gse188646_pseudobulk_counts.mtx", "Pseudobulk counts matrix"),
    ("outputs/exploratory_crossmodal_mirna_aging_summary.json", "Cross-modal summary"),
    ("outputs/exploratory_allen_aging_spatial_validation_summary.json", "RSTE3 validation summary"),
    ("outputs/figures/sa_bundle/", "Figure PNG bundle directory"),
    ("outputs/figures/sa_bundle_pdf/", "Figure PDF bundle (300 dpi embedding)"),
    ("requirements.txt", "Python dependencies"),
    ("run_extended.py", "Main reproducibility entry point"),
]


def _md5(path: Path, max_bytes: int = 50_000_000) -> str:
    if not path.is_file() or path.stat().st_size > max_bytes:
        return ""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_manifest(out_dir: Path) -> tuple[Path, Path]:
    out_dir = out_dir.resolve()
    zen_dir = out_dir / "zenodo"
    zen_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for rel, desc in ZENODO_INCLUDE:
        p = (PROJECT_ROOT / rel).resolve()
        if p.is_dir():
            for f in sorted(p.rglob("*")):
                if f.is_file():
                    rows.append(
                        {
                            "relative_path": str(f.relative_to(PROJECT_ROOT)),
                            "description": desc,
                            "bytes": f.stat().st_size,
                            "md5_partial": _md5(f),
                            "present": "yes",
                        }
                    )
            continue
        rows.append(
            {
                "relative_path": rel,
                "description": desc,
                "bytes": p.stat().st_size if p.is_file() else 0,
                "md5_partial": _md5(p) if p.is_file() else "",
                "present": "yes" if p.is_file() else "no",
            }
        )
    df = pd.DataFrame(rows)
    csv_path = zen_dir / "ZENODO_DEPOSIT_MANIFEST.csv"
    df.to_csv(csv_path, index=False)
    readme = zen_dir / "ZENODO_DEPOSIT_README.txt"
    readme.write_text(
        "\n".join(
            [
                f"Zenodo deposit manifest — {date.today().isoformat()}",
                "",
                "Upload PROJECT_ROOT files listed in ZENODO_DEPOSIT_MANIFEST.csv.",
                "Large RDS/BIL files are intentionally excluded from Git.",
                "",
                f"Present: {(df['present']=='yes').sum()}/{len(df)} entries",
                "",
                "Suggested Zenodo metadata:",
                "  Title: Hypothalamic miRNA targetome × aging transcriptome integrative reanalysis",
                "  Resource type: Dataset + Software",
                "  License: CC-BY-4.0 (adjust per journal policy)",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return csv_path, readme


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs-dir", type=Path, default=PROJECT_ROOT / "outputs")
    args = ap.parse_args()
    csv, readme = build_manifest(args.outputs_dir)
    print(f"Wrote {csv}")
    print(f"Wrote {readme}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
