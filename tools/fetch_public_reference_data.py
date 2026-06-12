"""
CLI entry point for programmatic fetch of public reference artifacts.

Delegates to src/data_acquisition.py (same code path as run_extended.py hooks).

Usage (from feasibility_study/):
  python tools/fetch_public_reference_data.py --public-gmts
  python tools/fetch_public_reference_data.py --gse188646-rds
  python tools/fetch_public_reference_data.py --hypomap-h5ad
  python tools/fetch_public_reference_data.py --all
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))


def main() -> int:
    p = argparse.ArgumentParser(description="Fetch public reference data (GEO, GMTs, CELLxGENE).")
    p.add_argument("--public-gmts", action="store_true", help="miRTarBase + MSigDB Hallmark GMTs to data/cache")
    p.add_argument("--gse188646-rds", action="store_true", help="GSE188646 integrated Seurat RDS to data/")
    p.add_argument("--hypomap-h5ad", action="store_true", help="HypoMap CELLxGENE h5ad to data/references/")
    p.add_argument("--all", action="store_true", help="All of the above")
    args = p.parse_args()
    if not (args.all or args.public_gmts or args.gse188646_rds or args.hypomap_h5ad):
        p.print_help()
        return 2

    import data_acquisition

    def log(msg: str) -> None:
        print(msg, flush=True)

    ok = True
    if args.all or args.public_gmts:
        try:
            data_acquisition.ensure_public_gmts(log)
        except Exception as e:
            log(f"ERROR: public GMTs: {e}")
            ok = False
    if args.all or args.gse188646_rds:
        os.environ["AUTO_FETCH_GSE188646_RDS"] = "1"
        if not data_acquisition.ensure_gse188646_rds_download(log):
            ok = False
    if args.all or args.hypomap_h5ad:
        os.environ["AUTO_FETCH_HYPOMAP_H5AD"] = "1"
        if not data_acquisition.ensure_hypomap_h5ad_download(log):
            ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
