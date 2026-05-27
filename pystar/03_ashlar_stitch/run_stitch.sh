#!/bin/bash
#SBATCH -J Stitch_Alshlar
#SBATCH -o stitch_logs/log_Stitch_Alshlar_%x_%A.out
#SBATCH -e stitch_logs/log_Stitch_Alshlar_%x_%A.err
#SBATCH -p C64M512G
#SBATCH --qos=normal

#SBATCH -n 1
#SBATCH -c 60

#SBATCH --mem=64G
#SBATCH --time=24:00:00

#SBATCH --no-requeue
#SBATCH --export=ALL


start_time=$(date +%s)
echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"

# sbatch run_stitch.sh \
#     /gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mFAD/03_segmentation \
#     /gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mFAD/04_stitch \
#     /gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/output_2025-11-01_20251101A2B2_14mWT_14mAD.maf \
#     0.0946 \
#     ch00 \
#     90 \
#     /gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mFAD/02_registration/IF/PI \
#     128 \
#     254 \
#     0 \
#     8

source /gpfs/share/home/2301920002/software/miniconda3/etc/profile.d/conda.sh
conda activate pystar310
export JAVA_HOME="${CONDA_PREFIX}"
export PATH="${JAVA_HOME}/bin:${PATH}"
JVM_PATH="$(find "${JAVA_HOME}" -name libjvm.so -type f | head -n 1)"
if [[ -z "${JVM_PATH}" ]]; then
    echo "ERROR: libjvm.so not found under JAVA_HOME=${JAVA_HOME}" >&2
    echo "Run this in the pystar310 environment to locate it: find \"\$CONDA_PREFIX\" -name libjvm.so -type f" >&2
    exit 2
fi
export JVM_PATH

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
ASHLAR_SCRIPT="${SCRIPT_DIR}/main.py"


echo "============= SLURM Job Info =================="
echo "Job ID:          $SLURM_JOB_ID"
echo "Job Name:        $SLURM_JOB_NAME"
echo "User:            $SLURM_JOB_USER"
echo "Submit Host:     $SLURM_SUBMIT_HOST"
echo "Submit Directory:$SLURM_SUBMIT_DIR"
echo "Node List:       $SLURM_NODELIST"
echo "Job Node:        $SLURMD_NODENAME"
echo "Number of Nodes: $SLURM_JOB_NUM_NODES"
echo "Partition:       $SLURM_JOB_PARTITION"

echo "============= Allocated CPUs Info ============="
echo "CPUs per task:   $SLURM_CPUS_PER_TASK"
echo "Allocated CPUs:  $SLURM_JOB_CPUS_PER_NODE"

echo "Tasks per node:  $SLURM_NTASKS_PER_NODE"
echo "Total Tasks:     $SLURM_NTASKS"
echo "Memory per node: $SLURM_MEM_PER_NODE MB"
echo "==============================================="



INPUT_DIR=$1
OUTPUT_DIR=$2
MAF_FILE=$3
PIXEL_SIZE=$4
REF_CHANNEL=$5
ROTATION=$6 
IMAGE_DIR="${7:-}"
MAF_START="${8:-}"
MAF_END="${9:-}"
SPOT_ROTATION="${10:-}" 
LOCAL_QC_COUNT="${11:-8}"

EXTRA_ARGS=()
if [[ -n "${IMAGE_DIR}" ]]; then
    EXTRA_ARGS+=(--image_dir "${IMAGE_DIR}")
fi
if [[ -n "${MAF_START}" || -n "${MAF_END}" ]]; then
    if [[ -z "${MAF_START}" || -z "${MAF_END}" ]]; then
        echo "ERROR: MAF_START and MAF_END must be provided together." >&2
        exit 2
    fi
    EXTRA_ARGS+=(--maf_start "${MAF_START}" --maf_end "${MAF_END}")
fi
if [[ -n "${SPOT_ROTATION}" ]]; then
    EXTRA_ARGS+=(--spot_rotation "${SPOT_ROTATION}")
fi
if [[ -n "${LOCAL_QC_COUNT}" ]]; then
    EXTRA_ARGS+=(--local_qc_count "${LOCAL_QC_COUNT}")
fi

python "$ASHLAR_SCRIPT" \
    --input_dir "$INPUT_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --maf_file "$MAF_FILE" \
    --pixelsize ${PIXEL_SIZE} \
    --ref_channel "$REF_CHANNEL" \
    --rotation "$ROTATION" \
    "${EXTRA_ARGS[@]}"


echo "==============================================="

end_time=$(date +%s)
echo "End time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Elapsed time: $(($end_time - $start_time)) seconds"
