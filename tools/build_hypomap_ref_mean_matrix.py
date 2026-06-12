"""
Build a wide gene × cell_type mean expression matrix from an AnnData / h5ad reference
(e.g. HypoMap after you download it locally), for use with:

  set GSE188646_HYPOMAP_REF_EXPR_CSV=...   (Windows)
  Rscript r/gse188646_hypomap_reference_mapping.R <seurat.rds>

The R script expects CSV columns: `gene` (symbol) + one numeric column per cell type
(mean expression in the same scale as GSE188646 RNA `data` layer is ideal; Spearman
is scale-invariant but cross-study shifts still matter).

HypoMap (CELLxGENE Discover) — recommended workflow
  1) python tools/download_hypomap_cellxgene_h5ad.py
  2) python tools/build_hypomap_ref_mean_matrix.py inspect --h5ad data/references/cellxgene_hypomap.h5ad
  3) Build with harmonized C66 taxonomy, symbols from feature_name, means of normalized X:
     python tools/build_hypomap_ref_mean_matrix.py build \\
       --h5ad data/references/cellxgene_hypomap.h5ad \\
       --celltype-column C66_named --var-gene-col feature_name \\
       --output data/references/hypomap_cellxgene_C66_named_mean_X.csv \\
       --min-cells-per-type 50 --max-types 80
  4) Finer C185 resolution (≥200 cells per type, 170 types): \\
       --celltype-column C185_named \\
       --output data/references/hypomap_cellxgene_C185_named_mean_X_min200.csv \\
       --min-cells-per-type 200 --max-types 200
     Rationale: data/provenance/HypoMap_CELLxGENE_REFERENCE_MATRIX.txt

Other h5ad files:

  # Inspect obs/var without loading the full matrix into RAM aggressively
  python tools/build_hypomap_ref_mean_matrix.py inspect --h5ad D:\\refs\\hypomap.h5ad

  # Use a specific expression layer if present
  python tools/build_hypomap_ref_mean_matrix.py build --h5ad path.h5ad \\
    --celltype-column subclass --layer logcounts --output out.csv

  # Tiny synthetic CSV to test the R mapping path (no h5ad required)
  python tools/build_hypomap_ref_mean_matrix.py demo --output data/references/demo_ref_means.csv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def _safe_colname(s: str) -> str:
    s = str(s).strip()
    s = re.sub(r"[^\w]+", "_", s, flags=re.UNICODE)
    s = re.sub(r"_+", "_", s).strip("_")
    if not s or s[0].isdigit():
        s = "t_" + s
    return s[:120] if len(s) > 120 else s


def cmd_inspect(h5ad: Path) -> int:
    import anndata as ad

    adata = ad.read_h5ad(h5ad, backed="r")
    print("=== obs columns (first 40) ===")
    print(list(adata.obs.columns)[:40])
    print("\n=== var columns (first 20) ===")
    print(list(adata.var.columns)[:20])
    print("\n=== shape ===", adata.shape)
    print("\n=== X type ===", type(adata.X))
    if adata.layers:
        print("\n=== layers ===", list(adata.layers.keys())[:30])
    if adata.raw is not None:
        print("\n=== raw is attached (often raw counts); primary X is usually normalized ===")
    print("\n=== candidate obs columns for --celltype-column (5 <= nunique <= 250) ===")
    for c in adata.obs.columns:
        try:
            nu = adata.obs[c].astype(str).nunique()
        except Exception:
            continue
        if 5 <= nu <= 250:
            print(f"  {c!r}  nunique={nu}")
    print("\nTip: prefer harmonized biology labels (e.g. cell_type); use finer columns only if")
    print("      n_cells per type stays high enough for stable means (see --min-cells-per-type).")
    return 0


def cmd_demo(output: Path) -> int:
    import numpy as np
    import pandas as pd

    root = Path(__file__).resolve().parents[1]
    deg = root / "outputs" / "gse188646_young_vs_aged_deg.csv"
    if deg.is_file():
        g = pd.read_csv(deg, nrows=400)["gene"].astype(str).tolist()
    else:
        g = [f"Gm{i}" for i in range(300)]
    rng = np.random.default_rng(42)
    types = ["RefA", "RefB", "RefC", "RefD", "RefE"]
    df = pd.DataFrame({"gene": g})
    for t in types:
        df[t] = rng.normal(scale=0.5, size=len(g))
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    print("Wrote", output.resolve(), "n_genes=", len(g))
    return 0


def cmd_build(
    h5ad: Path,
    celltype_col: str,
    output: Path,
    layer: str | None,
    var_gene_col: str | None,
    min_cells: int,
    max_types: int,
    gzip_out: bool,
) -> int:
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import issparse

    adata = ad.read_h5ad(h5ad, backed="r")
    if celltype_col not in adata.obs.columns:
        print(f"ERROR: obs column not found: {celltype_col!r}", file=sys.stderr)
        print("Available:", list(adata.obs.columns), file=sys.stderr)
        return 1

    ct = adata.obs[celltype_col].astype(str)
    vc = ct.value_counts()
    vc = vc[vc >= min_cells].sort_values(ascending=False)
    types = list(vc.index[:max_types])
    if len(types) < 2:
        print("ERROR: fewer than 2 cell types after filters.", file=sys.stderr)
        return 1

    genes = adata.var_names.astype(str)
    if var_gene_col and var_gene_col in adata.var.columns:
        genes = adata.var[var_gene_col].astype(str)

    def get_x(sub):
        if layer and layer in sub.layers:
            return sub.layers[layer]
        return sub.X

    means: dict[str, np.ndarray] = {}
    for t in types:
        idx = np.flatnonzero((ct == t).values)
        if idx.size == 0:
            continue
        sub = adata[idx].to_memory()
        X = get_x(sub)
        if issparse(X):
            mu = np.asarray(X.mean(axis=0)).ravel()
        else:
            mu = np.asarray(np.mean(X, axis=0)).ravel()
        if mu.shape[0] != len(genes):
            print("ERROR: mean length mismatch with var.", file=sys.stderr)
            return 1
        means[_safe_colname(t)] = mu.astype(np.float64)
        del sub

    mat = pd.DataFrame(means, index=genes.astype(str))
    mat.index.name = "gene"
    mat = mat.reset_index()
    gcol = mat.columns[0]
    rest = [c for c in mat.columns if c != gcol]
    if mat[gcol].duplicated().any():
        mat = mat.groupby(gcol, as_index=False)[rest].mean(numeric_only=True)
    mat = mat.rename(columns={gcol: "gene"})

    output.parent.mkdir(parents=True, exist_ok=True)
    if gzip_out:
        out_path = output.with_suffix(output.suffix + ".gz") if not str(output).endswith(".gz") else output
        mat.to_csv(out_path, index=False, compression="gzip")
        print("Wrote", out_path.resolve(), "shape", mat.shape)
    else:
        mat.to_csv(output, index=False)
        print("Wrote", output.resolve(), "shape", mat.shape)
    print("Columns:", list(mat.columns)[:12], "..." if mat.shape[1] > 12 else "")
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description="Build gene × cell_type mean matrix from h5ad for HypoMap-style mapping.")
    sub = p.add_subparsers(dest="cmd", required=True)

    pi = sub.add_parser("inspect", help="List obs/var columns and shape (backed read).")
    pi.add_argument("--h5ad", type=Path, required=True)

    pdemo = sub.add_parser("demo", help="Write a tiny random CSV for testing the R pipeline.")
    pdemo.add_argument("--output", type=Path, default=Path("data/references/demo_ref_celltype_means.csv"))

    pb = sub.add_parser("build", help="Aggregate mean expression per cell type.")
    pb.add_argument("--h5ad", type=Path, required=True)
    pb.add_argument("--celltype-column", type=str, required=True, help="adata.obs column (categorical or string).")
    pb.add_argument("--output", type=Path, required=True, help="Output CSV path (wide: gene + one column per type).")
    pb.add_argument("--layer", type=str, default=None, help="Use adata.layers[layer]; default is .X")
    pb.add_argument(
        "--var-gene-col",
        type=str,
        default=None,
        help="If set, take gene symbols from adata.var[this] instead of var_names.",
    )
    pb.add_argument("--min-cells-per-type", type=int, default=50)
    pb.add_argument("--max-types", type=int, default=60, help="Keep the largest N types by cell count after min filter.")
    pb.add_argument("--gzip", action="store_true", help="Write .csv.gz")

    args = p.parse_args()
    if args.cmd == "inspect":
        return cmd_inspect(args.h5ad)
    if args.cmd == "demo":
        return cmd_demo(args.output)
    if args.cmd == "build":
        return cmd_build(
            args.h5ad,
            args.celltype_column,
            args.output,
            args.layer,
            args.var_gene_col,
            args.min_cells_per_type,
            args.max_types,
            args.gzip,
        )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
