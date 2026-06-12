"""
Download and decompress the official GSE188646 integrated Seurat object from NCBI GEO FTP.

Source (canonical NLM HTTPS mirror of GEO supplementary):
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188646/suppl/
  GSE188646_hypo.integrated.final.20210719.RDS.gz

GEO distributes the object gzip-compressed; Seurat readRDS() expects the decompressed .RDS.
This script streams the .gz, verifies size against the server Content-Length when provided,
decompresses atomically to --output, then optionally runs R readRDS() as a smoke test.

Rationale: use the official FTP hostname and series-level path (GSE188nnn/GSE188646) so the
artifact is unambiguously the dataset publisher file, not a third-party mirror.
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import requests

GSE188646_RDS_GZ_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188646/suppl/"
    "GSE188646_hypo.integrated.final.20210719.RDS.gz"
)


def _disk_free_bytes(path: Path) -> int:
    return shutil.disk_usage(path).free


def download_stream(url: str, dest_gz: Path, chunk: int = 4 << 20) -> int:
    dest_gz.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest_gz.with_suffix(dest_gz.suffix + ".partial")
    expected: int | None = None
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        cl = r.headers.get("Content-Length")
        if cl and cl.isdigit():
            expected = int(cl)
        n = 0
        last_100mb = -1
        with open(tmp, "wb") as f:
            for block in r.iter_content(chunk_size=chunk):
                if block:
                    f.write(block)
                    n += len(block)
                    m100 = n // (100 << 20)
                    if m100 > last_100mb:
                        last_100mb = m100
                        print(f"  downloaded {n / (1 << 30):.2f} GiB...", flush=True)
        if expected is not None and n != expected:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(f"Size mismatch: got {n} bytes, expected {expected}")
    tmp.replace(dest_gz)
    return n


def gunzip_to_rds(gz_path: Path, rds_path: Path) -> None:
    rds_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = rds_path.with_suffix(rds_path.suffix + ".partial")
    try:
        with gzip.open(gz_path, "rb") as f_in, open(tmp, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out, length=16 << 20)
        tmp.replace(rds_path)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise


def smoke_read_rds(rds_path: Path) -> None:
    """Full-object readRDS; use only if you have sufficient RAM (multi-GB Seurat object)."""
    cmd = [
        "Rscript",
        "-e",
        f"invisible(readRDS('{rds_path.as_posix()}')); message('readRDS_OK')",
    ]
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=7200,
        encoding="utf-8",
        errors="replace",
    )
    out = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0 or "readRDS_OK" not in out:
        raise RuntimeError(f"R readRDS smoke test failed (exit {proc.returncode}): {out[-2000:]}")


def main() -> int:
    p = argparse.ArgumentParser(description="Download GSE188646 integrated Seurat RDS from GEO FTP.")
    p.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "data"
        / "GSE188646_hypo.integrated.final.20210719.RDS",
        help="Decompressed .RDS destination path",
    )
    p.add_argument(
        "--keep-gz",
        action="store_true",
        help="Keep the .gz after successful decompression (default: delete to save space)",
    )
    p.add_argument(
        "--r-smoke-test",
        action="store_true",
        help="After decompress, load the entire RDS in R (readRDS). Very slow and RAM-heavy; off by default.",
    )
    p.add_argument(
        "--min-free-gib",
        type=float,
        default=8.0,
        help="Abort if less than this many GiB free on destination volume (compressed + decompressed peaks).",
    )
    args = p.parse_args()
    out_rds: Path = args.output.resolve()
    out_gz = out_rds.with_suffix(out_rds.suffix + ".gz")

    if out_rds.is_file() and out_rds.stat().st_size > 100_000_000:
        print(f"Already present ({out_rds.stat().st_size} bytes): {out_rds}", flush=True)
        return 0

    free = _disk_free_bytes(out_rds.parent)
    need = int(args.min_free_gib * (1 << 30))
    if free < need:
        print(
            f"Insufficient disk space: {free / (1 << 30):.1f} GiB free, "
            f"recommended >= {args.min_free_gib} GiB for download + decompress.",
            file=sys.stderr,
        )
        return 2

    print(f"Downloading:\n  {GSE188646_RDS_GZ_URL}\n-> {out_gz}", flush=True)
    n = download_stream(GSE188646_RDS_GZ_URL, out_gz)
    print(f"Downloaded {n} bytes compressed.", flush=True)

    print(f"Decompressing gzip -> {out_rds}", flush=True)
    gunzip_to_rds(out_gz, out_rds)
    if not args.keep_gz:
        out_gz.unlink(missing_ok=True)

    rds_size = out_rds.stat().st_size
    if rds_size < 500_000_000:
        raise RuntimeError(
            f"Decompressed RDS unexpectedly small ({rds_size} bytes); "
            "decompression may have failed or the upstream file changed."
        )

    if args.r_smoke_test:
        print("Running R readRDS smoke test (RAM-heavy; full object load)...", flush=True)
        smoke_read_rds(out_rds)
        print("R readRDS OK.", flush=True)

    print(f"Done: {out_rds} ({rds_size} bytes)", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        raise SystemExit(130) from None
