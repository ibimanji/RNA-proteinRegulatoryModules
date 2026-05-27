#!/bin/bash
#SBATCH -J Stitch_IF_Registered
#SBATCH -o stitch_logs/log_Stitch_IF_Registered_%x_%A.out
#SBATCH -e stitch_logs/log_Stitch_IF_Registered_%x_%A.err
#SBATCH -p C64M512G
#SBATCH --qos=normal
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --export=ALL

start_time=$(date +%s)
echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"

# sbatch run_stitch_registered_channels.sh \
#     /gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mFAD/04_stitch/registered_tilecoord.csv \
#     /gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mFAD/02_registration/IF \
#     /gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mFAD/04_stitch/IF_registered_channels \
#     90 \
#     PI X34

source /gpfs/share/home/2301920002/software/miniconda3/etc/profile.d/conda.sh
conda activate pystar310

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
SCRIPT="${SCRIPT_DIR}/stitch_registered_channels.py"

REGISTERED_TILECOORD=$1
IF_DIR=$2
OUTPUT_DIR=$3
ROTATION=$4
shift 4
CHANNELS=("$@")

if [[ -z "${REGISTERED_TILECOORD}" || -z "${IF_DIR}" || -z "${OUTPUT_DIR}" || -z "${ROTATION}" || ${#CHANNELS[@]} -eq 0 ]]; then
    echo "Usage: sbatch run_stitch_registered_channels.sh REGISTERED_TILECOORD IF_DIR OUTPUT_DIR ROTATION CHANNEL [CHANNEL ...]" >&2
    exit 2
fi

mkdir -p "${OUTPUT_DIR}"

echo "registered_tilecoord: ${REGISTERED_TILECOORD}"
echo "IF_DIR:               ${IF_DIR}"
echo "OUTPUT_DIR:           ${OUTPUT_DIR}"
echo "rotation clockwise:   ${ROTATION}"
echo "channels:             ${CHANNELS[*]}"

python "${SCRIPT}" \
    --registered_tilecoord "${REGISTERED_TILECOORD}" \
    --if_dir "${IF_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --rotation "${ROTATION}" \
    --channels "${CHANNELS[@]}"

echo "==============================================="
end_time=$(date +%s)
echo "End time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Elapsed time: $(($end_time - $start_time)) seconds"
