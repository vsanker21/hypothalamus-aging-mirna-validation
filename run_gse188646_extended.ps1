# Run extended pipeline with GSE188646: RDS download (optional), full pseudobulk from RDS,
# exported counts (Tier-2), stratified pseudobulk by default, then full run_extended.py.
#
# Typical run (from feasibility_study/):
#   .\run_gse188646_extended.ps1
#
# First run: downloads the official RDS (large), verifies size (>= 1 GiB), re-runs pseudobulk,
# exports counts (.mtx + sidecar CSVs), runs stratified DE under outputs/gse188646_strata/,
# runs python run_extended.py, then optional cluster marker mapping (outputs/gse188646_cluster_annotation/).
#
# Usage (from feasibility_study/):
#   .\run_gse188646_extended.ps1 -Download          # force re-fetch .RDS.gz from GEO
#   .\run_gse188646_extended.ps1 -RdsPath "D:\path\to\GSE188646_hypo.integrated.final.20210719.RDS"
#   .\run_gse188646_extended.ps1 -NoStratified       # skip per-cluster pseudobulk
#   .\run_gse188646_extended.ps1 -StratumCol ""     # same as -NoStratified
#   .\run_gse188646_extended.ps1 -NoAutoDownload   # fail if RDS missing (no implicit fetch)
#   .\run_gse188646_extended.ps1 -StratumCol "RNA_snn_res.0.5"   # if seurat_clusters absent
#   .\run_gse188646_extended.ps1 -ForcePseudobulkFromRds:$false # reuse existing DE if present
#   .\run_gse188646_extended.ps1 -NoClusterMap              # skip r/gse188646_cluster_marker_mapping.R
#   .\run_gse188646_extended.ps1 -NoHypomapMap              # skip GSE208355 axis correlation mapping
#   .\run_gse188646_extended.ps1 -HypomapRefExprCsv "data\references\hypomap_cellxgene_C185_named_mean_X_min200.csv"
#   (C66 overview: data\references\hypomap_cellxgene_C66_named_mean_X.csv)
#       sets GSE188646_HYPOMAP_REF_EXPR_CSV so r/gse188646_hypomap_reference_mapping.R can emit hypomap_custom_ref_spearman.csv
#
# Defaults: if RDS is absent, downloads from official GEO supplementary URL (same as tools\download_gse188646_rds.py).
#           GSE188646_FORCE_PSEUDOBULK=1 so pseudobulk always re-runs from RDS when this launcher is used.
#           GSE188646_STRATUM_COL=seurat_clusters for r/pseudobulk_stratified_edgeR_gse188646.R
#           (override if your object uses another column, e.g. RNA_snn_res.0.5).
#           run_extended.py runs GSE87102 cohort2 + DL meta and the SA figure PNG bundle by default
#           (opt-out: $env:GSE188646_SKIP_COHORT2_META="1"; $env:SKIP_SA_FIGURE_BUNDLE="1").

param(
    [string] $RdsPath = $(Join-Path $PSScriptRoot "data\GSE188646_hypo.integrated.final.20210719.RDS"),
    [switch] $Download,
    [string] $StratumCol = "seurat_clusters",
    [bool] $AutoDownloadIfMissing = $true,
    [bool] $ForcePseudobulkFromRds = $true,
    [switch] $NoStratified,
    [switch] $NoAutoDownload,
    [switch] $NoClusterMap,
    [switch] $NoHypomapMap,
    [string] $HypomapRefExprCsv = ""
)

$ErrorActionPreference = "Stop"

$env:GSE188646_RDS = $RdsPath
$env:GSE188646_EXPORT_COUNTS = "1"
if ($ForcePseudobulkFromRds) {
    $env:GSE188646_FORCE_PSEUDOBULK = "1"
}
else {
    Remove-Item Env:\GSE188646_FORCE_PSEUDOBULK -ErrorAction SilentlyContinue
}

if ($NoStratified -or [string]::IsNullOrWhiteSpace($StratumCol)) {
    Remove-Item Env:\GSE188646_STRATUM_COL -ErrorAction SilentlyContinue
}
else {
    $env:GSE188646_STRATUM_COL = $StratumCol
}

Set-Location $PSScriptRoot
Write-Host "GSE188646_RDS=$env:GSE188646_RDS"
Write-Host "GSE188646_EXPORT_COUNTS=$env:GSE188646_EXPORT_COUNTS"
Write-Host "GSE188646_FORCE_PSEUDOBULK=$($env:GSE188646_FORCE_PSEUDOBULK)"
Write-Host "GSE188646_STRATUM_COL=$($env:GSE188646_STRATUM_COL)"

$needDownload = -not (Test-Path -LiteralPath $RdsPath)
if ($needDownload) {
    $doFetch = $Download -or ($AutoDownloadIfMissing -and -not $NoAutoDownload)
    if (-not $doFetch) {
        Write-Host "RDS not found at $RdsPath. Re-run with -Download, or omit -NoAutoDownload to fetch automatically." -ForegroundColor Red
        exit 1
    }
    Write-Host "Downloading GSE188646 RDS from NCBI GEO (official supplementary)..." -ForegroundColor Cyan
    $dl = Join-Path $PSScriptRoot "tools\download_gse188646_rds.py"
    python $dl --output $RdsPath
    if ($LASTEXITCODE -ne 0) {
        Write-Host "download_gse188646_rds.py failed (exit $LASTEXITCODE)." -ForegroundColor Red
        exit $LASTEXITCODE
    }
}

if (-not (Test-Path -LiteralPath $RdsPath)) {
    Write-Host "RDS still missing at $RdsPath after download attempt." -ForegroundColor Red
    exit 1
}

$len = (Get-Item -LiteralPath $RdsPath).Length
if ($len -lt 1GB) {
    Write-Host "RDS file is unexpectedly small ($len bytes); expected a multi-GiB Seurat object. Aborting." -ForegroundColor Red
    exit 1
}
Write-Host "RDS OK: $RdsPath ($([math]::Round($len/1GB, 2)) GiB)" -ForegroundColor Green

python run_extended.py
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

if (-not $NoClusterMap) {
    $mapR = Join-Path $PSScriptRoot "r\gse188646_cluster_marker_mapping.R"
    Write-Host "Cluster marker mapping (outputs/gse188646_cluster_annotation/) ..." -ForegroundColor Cyan
    Rscript $mapR $RdsPath
    if ($LASTEXITCODE -ne 0) {
        Write-Host "gse188646_cluster_marker_mapping.R failed (exit $LASTEXITCODE)." -ForegroundColor Red
        exit $LASTEXITCODE
    }
}

if (-not $NoHypomapMap) {
    if (-not [string]::IsNullOrWhiteSpace($HypomapRefExprCsv)) {
        $refAbs = if ([System.IO.Path]::IsPathRooted($HypomapRefExprCsv)) {
            $HypomapRefExprCsv
        } else {
            Join-Path $PSScriptRoot $HypomapRefExprCsv
        }
        if (-not (Test-Path -LiteralPath $refAbs)) {
            Write-Host "Hypomap ref CSV not found: $refAbs" -ForegroundColor Red
            exit 1
        }
        $env:GSE188646_HYPOMAP_REF_EXPR_CSV = (Resolve-Path -LiteralPath $refAbs).Path
        Write-Host "GSE188646_HYPOMAP_REF_EXPR_CSV=$env:GSE188646_HYPOMAP_REF_EXPR_CSV" -ForegroundColor Cyan
    } else {
        Remove-Item Env:\GSE188646_HYPOMAP_REF_EXPR_CSV -ErrorAction SilentlyContinue
    }
    $hypoR = Join-Path $PSScriptRoot "r\gse188646_hypomap_reference_mapping.R"
    Write-Host "HypoMap GSE208355 axis mapping (outputs/gse188646_hypomap_mapping/) ..." -ForegroundColor Cyan
    Rscript $hypoR $RdsPath
    if ($LASTEXITCODE -ne 0) {
        Write-Host "gse188646_hypomap_reference_mapping.R failed (exit $LASTEXITCODE)." -ForegroundColor Red
        exit $LASTEXITCODE
    }
}

exit 0
