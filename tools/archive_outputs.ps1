# Copy key pipeline outputs to outputs/archive/YYYY-MM-DD/ for publication freezes.
# Run from feasibility_study/:  powershell -File tools/archive_outputs.ps1

$root = Split-Path -Parent $PSScriptRoot
$dest = Join-Path $root "outputs\archive\$(Get-Date -Format 'yyyy-MM-dd')"
New-Item -ItemType Directory -Force -Path $dest | Out-Null
$files = @(
    "outputs\EXTENDED_REPORT.txt",
    "outputs\REPORT_SUMMARY.txt",
    "outputs\fisher_targets_vs_gse188646_de.csv",
    "outputs\fisher_targets_vs_hypomap_geo_deseq_union.csv",
    "outputs\gse188646_young_vs_aged_deg.csv",
    "outputs\cohort2_GSE87102_C57_hypothalamus_aged_vs_young_limma.csv",
    "outputs\exploratory_meta_DE_two_cohort_DL.csv",
    "outputs\VERSION_SNAPSHOT.txt",
    "outputs\gse188646_cluster_annotation\cluster_putative_labels.csv",
    "outputs\gse188646_hypomap_mapping\cluster_hypomap_axes_combined.csv",
    "outputs\gse188646_hypomap_mapping\hypomap_axis_spearman.csv",
    "outputs\gse188646_hypomap_mapping\hypomap_custom_ref_spearman.csv"
)
foreach ($rel in $files) {
    $p = Join-Path $root $rel
    if (Test-Path -LiteralPath $p) {
        Copy-Item -LiteralPath $p -Destination $dest -Force
    }
}
Write-Host "Archived to $dest"
