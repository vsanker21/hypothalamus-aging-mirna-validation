"""Load Zhang et al. Nature 2017 supplementary Excel tables from the workspace."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from config import MOESM_GLOB, WORKSPACE_ROOT


def list_moesm_files(root: Optional[Path] = None) -> list[Path]:
    root = root or WORKSPACE_ROOT
    return sorted(root.glob(MOESM_GLOB))


def read_first_sheet(path: Path) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name=0, header=None, engine=None)


def _pick(pattern_substr: str) -> Path:
    root = WORKSPACE_ROOT
    hits = [p for p in root.glob(MOESM_GLOB) if pattern_substr in p.name]
    if not hits:
        raise FileNotFoundError(f"No file matching {pattern_substr}")
    return hits[0]


def load_mirna_expression_fig4d(path: Optional[Path] = None) -> pd.DataFrame:
    """
    MOESM2: miRNA expression (htNSC vs hippocampal NSC vs astrocyte replicates).
    """
    path = path or _pick("MOESM2_ESM")
    raw = read_first_sheet(path)
    header = raw.iloc[1].tolist()
    df = raw.iloc[2:].copy()
    df.columns = [str(c).strip() if c is not None else "" for c in header]
    if "Transcript ID(Array Design)" in df.columns:
        df = df.rename(columns={"Transcript ID(Array Design)": "mirna"})
    else:
        df = df.rename(columns={df.columns[0]: "mirna"})
    df = df.dropna(subset=["mirna"])
    df["mirna"] = df["mirna"].astype(str).str.strip()
    for c in ["htNSC1", "htNSC2", "hippoNSC1", "hippoNSC2", "Astrocyte1", "Astrocyte2"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    htnsc = df[["htNSC1", "htNSC2"]].mean(axis=1)
    hippo = df[["hippoNSC1", "hippoNSC2"]].mean(axis=1)
    astro = df[["Astrocyte1", "Astrocyte2"]].mean(axis=1)
    return pd.DataFrame(
        {
            "mirna": df["mirna"],
            "htnsc_mean": htnsc,
            "hippo_nsc_mean": hippo,
            "astrocyte_mean": astro,
            "logfc_htnsc_vs_astro": htnsc - astro,
            "logfc_htnsc_vs_hippo": htnsc - hippo,
        }
    ).reset_index(drop=True)


def load_exosome_nanoparticle_fig4b(path: Optional[Path] = None) -> pd.DataFrame:
    """MOESM13: nanoparticle / exosome-related particle measures by cell type."""
    path = path or _pick("MOESM13_ESM")
    raw = read_first_sheet(path)
    header_row = 2
    cols = raw.iloc[header_row].tolist()
    body = raw.iloc[header_row + 1 :].copy()
    body.columns = [str(c).strip() if c is not None else f"c{i}" for i, c in enumerate(cols)]
    body = body.dropna(how="all")
    return body.reset_index(drop=True)


def load_csf_mirna_young_fig5a(path: Optional[Path] = None) -> pd.DataFrame:
    """MOESM14: CSF miRNA (young cohort)."""
    path = path or _pick("MOESM14_ESM")
    raw = read_first_sheet(path)
    sub = raw.iloc[4:].copy()
    sub = sub.dropna(how="all")
    mirna = sub.iloc[:, 0].astype(str)
    vals = sub.iloc[:, 1:7]
    vals.columns = [f"rep{i}" for i in range(vals.shape[1])]
    out = pd.concat([mirna.rename("mirna"), vals], axis=1)
    out = out[out["mirna"].str.contains(r"^[0-9A-Za-z\-]+", regex=True, na=False)]
    for c in out.columns:
        if c.startswith("rep"):
            out[c] = pd.to_numeric(out[c], errors="coerce")
    repcols = [c for c in out.columns if c.startswith("rep")]
    out["mean_expr"] = out[repcols].mean(axis=1)
    return out.reset_index(drop=True)


def load_sox2_counts_fig1d(path: Optional[Path] = None) -> pd.DataFrame:
    """MOESM10: Sox2+ counts in 3V wall and ARC (control vs TK)."""
    path = path or _pick("MOESM10_ESM")
    raw = read_first_sheet(path)
    rows = []
    for i in range(4, 9):
        r = raw.iloc[i]
        rows.append(
            {
                "animal_idx": int(r[0]) if pd.notna(r[0]) else i - 3,
                "sox2_3v_con": float(r[1]),
                "sox2_3v_tk": float(r[2]),
                "sox2_arc_con": float(r[3]),
                "sox2_arc_tk": float(r[4]),
            }
        )
    return pd.DataFrame(rows)


def load_exosome_rescue_phenotypes_fig6c(path: Optional[Path] = None) -> pd.DataFrame:
    """MOESM15: multi-endpoint phenotypes (control vs TK vs TK+exosome)."""
    path = path or _pick("MOESM15_ESM")
    raw = read_first_sheet(path)
    header = [str(h).strip() if h is not None else "" for h in raw.iloc[2].tolist()]
    body = raw.iloc[3:28].copy()
    body.columns = header[: body.shape[1]]
    body = body.loc[:, [c for c in body.columns if c]]
    return body.reset_index(drop=True)


def load_cytokine_exosome_long(path: Optional[Path] = None) -> pd.DataFrame:
    """MOESM22: cytokine mRNA vehicle vs exosome (long tidy). Column layout is fixed in source file."""
    path = path or _pick("MOESM22_ESM")
    raw = read_first_sheet(path)
    # Row index 3..7 animals; genes at fixed pairs (Veh, Exo)
    blocks = [
        ("TNF", 1, 2),
        ("IL6", 4, 5),
        ("IL1B", 7, 8),
        ("NFKB1", 10, 11),
        ("CXCL10", 13, 14),
    ]
    rows = []
    for ri in range(3, 8):
        r = raw.iloc[ri]
        if pd.isna(r[0]):
            continue
        aid = int(r[0])
        for gene, c0, c1 in blocks:
            rows.append(
                {
                    "animal": aid,
                    "gene": gene,
                    "vehicle": float(r[c0]),
                    "exosome": float(r[c1]),
                }
            )
    return pd.DataFrame(rows)

