"""
Download the curated HypoMap AnnData from CZ CELLxGENE Discover (stable dataset asset URL).

Source (collection + citation):
  https://cellxgene.cziscience.com/collections/d86517f0-fa7e-4266-b82e-a521350d6d36
  DOI: 10.1038/s42255-022-00657-y

Default asset (CELLxGENE API, 2025-10 revision):
  https://datasets.cellxgene.cziscience.com/d3be7423-d664-4913-89a9-a506cae4c28f.h5ad
  Expected size: 3814481296 bytes (~3.55 GiB)

CELLxGENE schema: primary matrix `.X` is the normalized / analysis matrix; raw counts are in
`adata.raw` when `raw_data_location` is `raw.X` (see collection metadata). For cross-study
Spearman vs GSE188646 Seurat `RNA$data`, use `.X` (default in build_hypomap_ref_mean_matrix.py).

Usage (from feasibility_study/):
  python tools/download_hypomap_cellxgene_h5ad.py
  python tools/download_hypomap_cellxgene_h5ad.py --output data/references/cellxgene_hypomap.h5ad
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import requests

DEFAULT_URL = (
    "https://datasets.cellxgene.cziscience.com/"
    "d3be7423-d664-4913-89a9-a506cae4c28f.h5ad"
)
EXPECTED_BYTES = 3_814_481_296


def download_stream(url: str, dest: Path, expected: int | None, chunk: int = 4 << 20) -> int:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    n = 0
    last_gb = -1
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        cl = r.headers.get("Content-Length")
        exp = expected
        if exp is None and cl and cl.isdigit():
            exp = int(cl)
        with open(tmp, "wb") as f:
            for block in r.iter_content(chunk_size=chunk):
                if block:
                    f.write(block)
                    n += len(block)
                    gb = n // (1 << 30)
                    if gb > last_gb:
                        last_gb = gb
                        print(f"  downloaded {n / (1 << 30):.2f} GiB...", flush=True)
        if exp is not None and n != exp:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(f"Size mismatch: got {n} bytes, expected {exp}")
    tmp.replace(dest)
    return n


def main() -> int:
    p = argparse.ArgumentParser(description="Download HypoMap h5ad from CELLxGENE Discover.")
    p.add_argument("--url", type=str, default=DEFAULT_URL, help="Dataset asset URL (.h5ad)")
    p.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "references" / "cellxgene_hypomap.h5ad",
        help="Destination .h5ad path",
    )
    p.add_argument(
        "--no-size-check",
        action="store_true",
        help="Do not require exact byte match to EXPECTED_BYTES (not recommended).",
    )
    args = p.parse_args()
    out: Path = args.output.resolve()
    if out.is_file() and out.stat().st_size >= EXPECTED_BYTES - 10_000_000:
        print(f"Already present ({out.stat().st_size} bytes): {out}", flush=True)
        return 0

    free = shutil.disk_usage(out.parent).free
    if free < EXPECTED_BYTES + (500 << 20):
        print(
            f"Insufficient disk space: {free / (1 << 30):.2f} GiB free; need roughly "
            f"{EXPECTED_BYTES / (1 << 30):.2f} GiB + buffer.",
            file=sys.stderr,
        )
        return 2

    print(f"Downloading:\n  {args.url}\n-> {out}", flush=True)
    n = download_stream(args.url, out, None if args.no_size_check else EXPECTED_BYTES)
    print(f"Done: {out} ({n} bytes)", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        raise SystemExit(130) from None
