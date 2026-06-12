#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_knn_smooth_from_rds.sh TYPE K [options passed to R script]

Examples:
  bash scripts/run_knn_smooth_from_rds.sh Microglia 16

  bash run_knn_smooth_from_rds.sh 14mWT 16 \
    --rds states_with_plaque_info.rds \
    --d 20 \
    --outdir knn_smooth_outputs

Notes:
  This shell script orchestrates the run:
    1. R reads the Seurat RDS and exports totalRNA/counts for one metadata type.
    2. Python runs knn_smooth.py on the exported genes x cells matrix.

  Defaults are defined in scripts/run_knn_smooth_from_rds.R:
    --meta-col type
    --assay totalRNA
    --slot counts
    --python python3
EOF
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 2 ]]; then
  usage >&2
  exit 1
fi

TYPE_VALUE="$1"
K_VALUE="$2"
shift 2

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [[ -f "${REPO_DIR}/knn_smooth.py" ]]; then
  KNN_SCRIPT="${REPO_DIR}/knn_smooth.py"
elif [[ -f "${SCRIPT_DIR}/knn_smooth.py" ]]; then
  KNN_SCRIPT="${SCRIPT_DIR}/knn_smooth.py"
elif [[ -f "./knn_smooth.py" ]]; then
  KNN_SCRIPT="$(cd "$(dirname "./knn_smooth.py")" && pwd)/knn_smooth.py"
else
  echo "Error: could not find knn_smooth.py." >&2
  echo "Put knn_smooth.py next to this script, in the parent directory, or pass --knn-script PATH to the R script directly." >&2
  exit 1
fi

Rscript "${SCRIPT_DIR}/run_knn_smooth_from_rds.R" \
  --type "${TYPE_VALUE}" \
  --k "${K_VALUE}" \
  --knn-script "${KNN_SCRIPT}" \
  "$@"
