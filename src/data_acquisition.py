"""
Programmatic acquisition of public reference data used by run_extended.py.

Rationale (see data/provenance/DATA_ACQUISITION_AND_QC_RATIONALE.txt):
  - Prefer publisher / NLM mirrors (GEO FTP HTTPS) with Content-Length verification where available.
  - Small GMTs: Enrichr text endpoints (miRTarBase, MSigDB Hallmark) with post-download sanity checks.
  - Large artifacts (GSE188646 RDS, CELLxGENE h5ad): opt-in via env flags to avoid surprise multi-GB pulls.

Environment (optional):
  AUTO_FETCH_PUBLIC_DATA=1   — refresh / ensure miRTarBase + MSigDB Hallmark GMTs under data/cache.
  AUTO_FETCH_GSE188646_RDS=1 — run tools/download_gse188646_rds.py if canonical RDS missing or truncated.
  AUTO_FETCH_HYPOMAP_H5AD=1  — run tools/download_hypomap_cellxgene_h5ad.py if reference h5ad missing.

GSE188646_RDS resolution order (no env required if file is already in place):
  1) GSE188646_RDS if set and the path exists.
  2) Else data/GSE188646_hypo.integrated.final.20210719.RDS if present and >= 500 MiB (sanity on truncated downloads).
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from config import DATA_DIR, PROJECT_ROOT

CANONICAL_GSE188646_RDS = DATA_DIR / "GSE188646_hypo.integrated.final.20210719.RDS"
CANONICAL_HYPOMAP_H5AD = DATA_DIR / "references" / "cellxgene_hypomap.h5ad"
MIN_RDS_BYTES = 500_000_000
# Must match tools/download_hypomap_cellxgene_h5ad.py EXPECTED_BYTES (CELLxGENE asset revision).
EXPECTED_HYPOMAP_H5AD_BYTES = 3_814_481_296


def _truthy_env(name: str) -> bool:
    return os.environ.get(name, "").strip().lower() in ("1", "true", "yes")


def validate_gmt_file(path: Path, label: str) -> None:
    """Public wrapper: reject truncated or HTML error responses masquerading as GMT."""
    _validate_gmt_tabular(path, label)


def _validate_gmt_tabular(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label}: missing {path}")
    n = path.stat().st_size
    if n < 200:
        raise RuntimeError(f"{label}: file unexpectedly small ({n} bytes): {path}")
    head = path.read_text(encoding="utf-8", errors="replace")[:8000]
    if "\t" not in head:
        raise RuntimeError(f"{label}: GMT at {path} has no tab-separated fields (possible HTML error page)")


def ensure_public_gmts(log=print) -> None:
    """Download miRTarBase + MSigDB Hallmark GMTs if missing; validate non-degenerate content."""
    import mirna_target_union
    import msigdb_hallmark_gsea

    p1 = mirna_target_union.ensure_mirtarbase_gmt()
    validate_gmt_file(p1, "miRTarBase_2017")
    p2 = msigdb_hallmark_gsea.ensure_hallmark_gmt()
    validate_gmt_file(p2, "MSigDB_Hallmark_2020")
    log(f"Public GMTs OK: {p1.name}, {p2.name} (under data/cache).")


def validate_existing_public_gmts(log=print) -> None:
    """If GMT caches already exist, verify they are not HTML error pages (common on network failures)."""
    import mirna_target_union
    import msigdb_hallmark_gsea

    p1 = mirna_target_union.GMT_CACHE
    p2 = msigdb_hallmark_gsea.HALLMARK_GMT_CACHE
    checked = False
    if p1.is_file():
        validate_gmt_file(p1, "miRTarBase_2017")
        checked = True
    if p2.is_file():
        validate_gmt_file(p2, "MSigDB_Hallmark_2020")
        checked = True
    if checked:
        log("Existing public GMT caches validated (tabular GMT structure OK).")
    else:
        log("No GMT caches on disk yet; miRTarBase/Hallmark will download on first use.")


def run_tool_script(rel_tool: str, argv: list[str], log=print, timeout: int = 86_400) -> int:
    script = PROJECT_ROOT / rel_tool
    if not script.is_file():
        log(f"Acquisition: tool script not found: {script}")
        return 1
    cmd = [sys.executable, str(script), *argv]
    log(f"Acquisition: running {' '.join(cmd)}")
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            timeout=timeout,
            encoding="utf-8",
            errors="replace",
        )
    except subprocess.TimeoutExpired:
        log(f"Acquisition: timeout after {timeout}s: {rel_tool}")
        return 124
    except Exception as e:
        log(f"Acquisition: failed to run {rel_tool}: {e}")
        return 1
    if proc.returncode != 0:
        log(f"Acquisition: {rel_tool} exit {proc.returncode}")
        tail = (proc.stderr or proc.stdout or "")[-4000:]
        if tail.strip():
            log(tail)
        return int(proc.returncode)
    if proc.stdout and proc.stdout.strip():
        for line in proc.stdout.strip().splitlines()[-3:]:
            log(line)
    return 0


def ensure_gse188646_rds_download(log=print) -> bool:
    """
    If AUTO_FETCH_GSE188646_RDS=1 and canonical RDS missing or too small, stream from NCBI GEO FTP.
    Returns True if a usable file exists after the attempt.
    """
    if not _truthy_env("AUTO_FETCH_GSE188646_RDS"):
        return CANONICAL_GSE188646_RDS.is_file() and CANONICAL_GSE188646_RDS.stat().st_size >= MIN_RDS_BYTES
    if CANONICAL_GSE188646_RDS.is_file() and CANONICAL_GSE188646_RDS.stat().st_size >= MIN_RDS_BYTES:
        log("GSE188646 RDS already present; skipping download.")
        return True
    rc = run_tool_script(
        "tools/download_gse188646_rds.py",
        ["--output", str(CANONICAL_GSE188646_RDS)],
        log=log,
        timeout=86_400,
    )
    ok = rc == 0 and CANONICAL_GSE188646_RDS.is_file() and CANONICAL_GSE188646_RDS.stat().st_size >= MIN_RDS_BYTES
    if not ok:
        log("GSE188646 RDS download did not produce a valid file; set GSE188646_RDS manually or fix disk/network.")
    return ok


def ensure_hypomap_h5ad_download(log=print) -> bool:
    """If AUTO_FETCH_HYPOMAP_H5AD=1 and canonical h5ad missing, download from CELLxGENE asset URL."""
    if not _truthy_env("AUTO_FETCH_HYPOMAP_H5AD"):
        return CANONICAL_HYPOMAP_H5AD.is_file()
    if CANONICAL_HYPOMAP_H5AD.is_file() and CANONICAL_HYPOMAP_H5AD.stat().st_size >= EXPECTED_HYPOMAP_H5AD_BYTES - 10_000_000:
        log("HypoMap CELLxGENE h5ad already present; skipping download.")
        return True
    CANONICAL_HYPOMAP_H5AD.parent.mkdir(parents=True, exist_ok=True)
    rc = run_tool_script(
        "tools/download_hypomap_cellxgene_h5ad.py",
        ["--output", str(CANONICAL_HYPOMAP_H5AD)],
        log=log,
        timeout=86_400,
    )
    ok = rc == 0 and CANONICAL_HYPOMAP_H5AD.is_file()
    if ok:
        log(f"HypoMap h5ad OK: {CANONICAL_HYPOMAP_H5AD} ({CANONICAL_HYPOMAP_H5AD.stat().st_size} bytes).")
    return ok


def run_prefetch_hooks(log=print) -> None:
    """Called at the start of run_extended.main() when hooks are enabled."""
    if _truthy_env("AUTO_FETCH_PUBLIC_DATA"):
        try:
            ensure_public_gmts(log)
        except Exception as e:
            log(f"AUTO_FETCH_PUBLIC_DATA failed: {e}")
    if _truthy_env("AUTO_FETCH_GSE188646_RDS"):
        ensure_gse188646_rds_download(log)
    if _truthy_env("AUTO_FETCH_HYPOMAP_H5AD"):
        ensure_hypomap_h5ad_download(log)


def resolve_gse188646_rds_path() -> Path | None:
    """
    Path to GSE188646 Seurat object for R pseudobulk:
    env GSE188646_RDS if valid, else canonical data/ RDS if present and large enough.
    """
    raw = os.environ.get("GSE188646_RDS", "").strip()
    if raw:
        p = Path(raw).expanduser()
        if p.is_file():
            return p.resolve()
        return None
    if CANONICAL_GSE188646_RDS.is_file() and CANONICAL_GSE188646_RDS.stat().st_size >= MIN_RDS_BYTES:
        return CANONICAL_GSE188646_RDS.resolve()
    return None


def describe_rds_resolution(log=print) -> None:
    """Log how GSE188646 RDS was resolved (transparent provenance)."""
    raw = os.environ.get("GSE188646_RDS", "").strip()
    if raw:
        p0 = Path(raw).expanduser()
        if not p0.is_file():
            log(f"GSE188646_RDS is set but file not found (pseudobulk will be skipped): {raw}")
            return
    p = resolve_gse188646_rds_path()
    if p is None:
        return
    if raw and Path(raw).expanduser().resolve() == p:
        log(f"GSE188646 RDS: using GSE188646_RDS={p}")
    elif not raw:
        log(
            "GSE188646 RDS: GSE188646_RDS unset; using canonical local file "
            f"(place publisher .RDS under data/ or set env): {p}"
        )
