#!/bin/bash

#SBATCH -p C64M256G
#SBATCH -o final_%A.out
#SBATCH -e final_%A.err
#SBATCH -N 1
#SBATCH --no-requeue
#SBATCH -c 16

SAMPLE="14mWT"
STITCH_DIR="/gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mWT/04_stitch"
RESULTS_DIR="/gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mWT/05_final"
DAPI_MASK="${RESULTS_DIR}/maskdapi_segmentation.tif"
OUTPUT_TAG="20260527"

source /gpfs/share/home/2301920002/software/miniconda3/etc/profile.d/conda.sh
conda activate /gpfs/share/home/2301920002/software/miniconda3/envs/sopa

SCRIPT_DIR="/gpfs/share/home/2301920002/labShare/2301920002/ADdecon/RNA-proteinRegulatoryModules/pystar/04_final_output"

python "${SCRIPT_DIR}/final_mask2h5ad.py" \
  --sample "${SAMPLE}" \
  --stitch-dir "${STITCH_DIR}" \
  --results-dir "${RESULTS_DIR}" \
  --dapi-mask "${DAPI_MASK}" \
  --output-tag "${OUTPUT_TAG}"
