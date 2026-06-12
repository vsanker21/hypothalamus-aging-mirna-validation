# Hypothalamus aging × miRNA targetome — computational validation

Integrative reanalysis linking htNSC-biased miRNA target predictions to hypothalamic aging
transcriptomics (GSE188646/GSE87102), with external third-ventricle validation (Allen ABC aging
scRNA, MERFISH CCF anatomy, Jin RSTE3 BIL in situ spatial, GSE287202 Visium audit, Allen ISH).

**Repository:** https://github.com/vsanker21/hypothalamus-aging-mirna-validation  
**Author:** [vsanker21](https://github.com/vsanker21)

## Quick start

```powershell
cd feasibility_study
pip install -r requirements.txt
python run_pipeline.py          # MOESM / miRNA prerequisites
python run_extended.py          # full journal-tier + external validation spine
```

Key outputs land in `outputs/` (gitignored). Regenerate figures:

```powershell
python figures/build_sa_figure_bundle.py --outputs-dir outputs
python tools/export_figure_bundle_pdf.py --outputs-dir outputs
```

Submission helpers (after a local run):

```powershell
python tools/prepare_zenodo_deposit_manifest.py --outputs-dir outputs
python tools/build_zenodo_upload_bundle.py --outputs-dir outputs
python tools/build_nature_portfolio_reporting_checklist.py --outputs-dir outputs
```

## Large data (not in Git)

Download programmatically before `run_extended.py`:

| Asset | Tool / env |
|-------|------------|
| GSE188646 Seurat RDS (~1.4 GB) | `python tools/download_gse188646_rds.py` or `AUTO_FETCH_GSE188646_RDS=1` |
| RSTE3 BIL (doi:10.35077/g.1157) | auto-fetched by `bil_rstE3_fetch.py` on first run |
| Allen ABC MERFISH metadata | auto-fetched by `abc_atlas_fetch.py` |
| HypoMap h5ad (~3.6 GB) | `tools/download_hypomap_cellxgene_h5ad.py` (optional) |

Zenodo deposit manifest: `outputs/zenodo/ZENODO_DEPOSIT_MANIFEST.csv` (after one local run).

## Repository layout

- `src/` — analysis modules (cross-modal nulls, external validation, RSTE3, ISH, etc.)
- `r/` — GSE188646 pseudobulk edgeR / stratified DE
- `tools/` — downloads, manuscript build, GitHub/Zenodo helpers
- `figures/` — publication figure builders
- `data/provenance/` — SAP, figure manifest, submission readiness notes

## Citation

Public datasets: GSE188646, GSE87102, Allen Brain Cell Atlas, Brain Image Library g.1157 (RSTE3),
GSE287202, Allen Mouse Brain Atlas ISH.

## License

Code: MIT (adjust before release). Processed public data remain under their respective archive terms.
