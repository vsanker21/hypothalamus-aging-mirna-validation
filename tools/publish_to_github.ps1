# Publish feasibility_study to https://github.com/vsanker21
# Usage (from feasibility_study/):
#   .\tools\publish_to_github.ps1
#   .\tools\publish_to_github.ps1 -RepoName hypothalamus-aging-mirna-validation -Private

param(
    [string]$RepoName = "hypothalamus-aging-mirna-validation",
    [switch]$Private,
    [string]$GitExe = "C:\Program Files\Git\bin\git.exe",
    [string]$GhExe = "C:\Program Files\GitHub CLI\gh.exe"
)

$ErrorActionPreference = "Stop"
$Root = Split-Path $PSScriptRoot -Parent
Set-Location $Root

if (-not (Test-Path $GitExe)) { throw "git not found at $GitExe" }
if (-not (Test-Path $GhExe)) { throw "gh not found at $GhExe" }

function Invoke-Git { & $GitExe @Args }
function Invoke-Gh { & $GhExe @Args }

Write-Host "Checking GitHub CLI authentication..."
& $GhExe auth status 2>&1 | Out-String | Write-Host
if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "Not logged in. Run once interactively:" -ForegroundColor Yellow
    Write-Host "  gh auth login -h github.com -p https -w"
    Write-Host "Or set GH_TOKEN with repo scope, then re-run this script."
    exit 1
}

if (-not (Test-Path ".git")) {
    Write-Host "Initializing git repository..."
    Invoke-Git init -b main
}

Invoke-Git add -A
$status = Invoke-Git status --porcelain
if ($status) {
    Invoke-Git commit -m @"
Add hypothalamus aging miRNA integrative pipeline and V3 external validation.

Includes Allen/Jin/MERFISH/RSTE3/ISH modules, reproducibility docs, and submission artifacts.
Large reference files excluded per .gitignore; use download tools and Zenodo manifest.
"@
} else {
    Write-Host "Nothing to commit."
}

$visibility = if ($Private) { "--private" } else { "--public" }
$remote = "https://github.com/vsanker21/$RepoName.git"

if (-not (Invoke-Gh repo view "vsanker21/$RepoName" 2>$null)) {
    Write-Host "Creating GitHub repo vsanker21/$RepoName ..."
    Invoke-Gh repo create $RepoName $visibility --source=. --remote=origin --description "Hypothalamic miRNA targetome x aging transcriptomics with V3 niche external validation"
} else {
    Write-Host "Repo exists; ensuring remote origin..."
    $hasOrigin = Invoke-Git remote 2>$null | Select-String -Pattern "^origin$"
    if (-not $hasOrigin) {
        Invoke-Git remote add origin $remote
    }
}

Write-Host "Pushing to $remote ..."
Invoke-Git push -u origin main
Write-Host "Done: https://github.com/vsanker21/$RepoName" -ForegroundColor Green
