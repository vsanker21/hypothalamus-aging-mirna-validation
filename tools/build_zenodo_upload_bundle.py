"""
Build a Zenodo-ready upload zip from ZENODO_DEPOSIT_MANIFEST.csv (present files only).

Large RDS is included when present; use --exclude-rds to split uploads.

Outputs:
  outputs/zenodo/hypothalamus-aging-mirna-validation_YYYYMMDD.zip
  outputs/zenodo/ZENODO_UPLOAD_INSTRUCTIONS.txt
"""
from __future__ import annotations

import argparse
import zipfile
from datetime import date
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def build_zip(
    out_dir: Path,
    *,
    exclude_rds: bool = False,
    max_file_bytes: int | None = None,
) -> Path | None:
    out_dir = out_dir.resolve()
    manifest_csv = out_dir / "zenodo" / "ZENODO_DEPOSIT_MANIFEST.csv"
    if not manifest_csv.is_file():
        raise FileNotFoundError(f"Missing {manifest_csv}; run prepare_zenodo_deposit_manifest.py first.")

    df = pd.read_csv(manifest_csv)
    present = df[df["present"].astype(str).str.lower() == "yes"].copy()
    stamp = date.today().strftime("%Y%m%d")
    zip_path = out_dir / "zenodo" / f"hypothalamus-aging-mirna-validation_{stamp}.zip"

    included: list[str] = []
    skipped: list[str] = []
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for _, row in present.iterrows():
            rel = str(row["relative_path"]).replace("\\", "/")
            src = (PROJECT_ROOT / rel).resolve()
            if exclude_rds and rel.lower().endswith(".rds"):
                skipped.append(f"{rel} (exclude-rds)")
                continue
            if not src.is_file():
                skipped.append(f"{rel} (missing)")
                continue
            nbytes = src.stat().st_size
            if max_file_bytes is not None and nbytes > max_file_bytes:
                skipped.append(f"{rel} (>{max_file_bytes} bytes)")
                continue
            zf.write(src, arcname=rel)
            included.append(rel)

    instr = out_dir / "zenodo" / "ZENODO_UPLOAD_INSTRUCTIONS.txt"
    instr.write_text(
        "\n".join(
            [
                f"Zenodo upload bundle — {date.today().isoformat()}",
                "",
                f"Archive: {zip_path.name}",
                f"Files included: {len(included)}",
                f"Files skipped: {len(skipped)}",
                "",
                "Steps:",
                "  1. Create deposit at https://zenodo.org/deposit/new",
                "  2. Upload this zip (and RDS separately if excluded or >50 GB limit)",
                "  3. Metadata (suggested):",
                "       Title: Hypothalamic miRNA targetome × aging transcriptome integrative reanalysis",
                "       Resource type: Dataset",
                "       License: CC-BY-4.0",
                "       Related identifier: https://github.com/vsanker21/hypothalamus-aging-mirna-validation",
                "  4. Publish and paste DOI into manuscript + NATURE_PORTFOLIO_REPORTING_SUMMARY.md",
                "",
                "Included:",
                *[f"  + {p}" for p in included[:200]],
                *(["  ..."] if len(included) > 200 else []),
                "",
                "Skipped:",
                *[f"  - {s}" for s in skipped[:50]],
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {zip_path} ({len(included)} files)")
    print(f"Instructions: {instr}")
    return zip_path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs-dir", type=Path, default=PROJECT_ROOT / "outputs")
    ap.add_argument("--exclude-rds", action="store_true", help="Omit GSE188646 RDS from zip")
    ap.add_argument("--max-file-bytes", type=int, default=None)
    args = ap.parse_args()
    try:
        build_zip(args.outputs_dir, exclude_rds=args.exclude_rds, max_file_bytes=args.max_file_bytes)
    except FileNotFoundError as e:
        print(e)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
