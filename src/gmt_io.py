"""Download and parse Enrichr-style GMT (tab-separated: term, blank, gene, gene, ...)."""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Callable, Iterator

import requests

from config import DATA_DIR

CACHE_DIR = DATA_DIR / "cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)


def download_text(url: str, dest: Path, timeout: int = 300) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    f.write(chunk)
    # Ensure trailing newline so R fgsea::gmtPathways does not warn on "incomplete final line"
    with open(dest, "rb+") as f:
        f.seek(0, 2)
        if f.tell() == 0:
            return
        f.seek(-1, 2)
        last = f.read(1)
        if last != b"\n":
            f.write(b"\n")


def iter_gmt_lines(path: Path) -> Iterator[tuple[str, list[str]]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    mode = "rt" if str(path).endswith(".gz") else "r"
    with opener(path, mode, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            term = parts[0].strip()
            genes = [p.strip() for p in parts[2:] if p.strip()]
            yield term, genes


def build_mirna_to_targets_mmu(
    gmt_path: Path, mirna_filter: Callable[[str], bool]
) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for term, genes in iter_gmt_lines(gmt_path):
        if not term.startswith("mmu-miR") and not term.startswith("mmu-let"):
            continue
        if not mirna_filter(term):
            continue
        out.setdefault(term, set()).update(g for g in genes if g)
    return out
