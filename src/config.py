"""Paths and constants for the feasibility computational study."""
from pathlib import Path

# feasibility_study/ directory (parent of this file's package)
PROJECT_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = PROJECT_ROOT.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
DATA_DIR = PROJECT_ROOT / "data"

MOESM_GLOB = "41586_2017_BFnature23282_MOESM*.xlsx"

# snRNA-seq hypothalamus young vs aged female mice (Nature Aging companion dataset)
GEO_HYPOTHALAMUS_AGING = "GSE188646"

GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
