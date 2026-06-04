#!/bin/bash

# =================================================================
# 自动读取 Config，计算 Array 长度，提交任务
# 用法:
#   bash scripts/run_pystar.sh [config_path] [main|if|all]
# =================================================================

set -euo pipefail

CONFIG_FILE="${1:-config/experiment_config.yaml}"
RUN_MODE="${2:-all}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# 你的 conda 安装位置
CONDA_ROOT="$HOME/software/miniconda3"
CONDA_ENV_NAME="pystar310"

if [[ "$RUN_MODE" != "main" && "$RUN_MODE" != "if" && "$RUN_MODE" != "all" ]]; then
    echo "Error: RUN_MODE must be one of: main, if, all"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file not found at $CONFIG_FILE"
    exit 1
fi

CONFIG_FILE="$(cd "$(dirname "$CONFIG_FILE")" && pwd)/$(basename "$CONFIG_FILE")"

if [ ! -f "${CONDA_ROOT}/etc/profile.d/conda.sh" ]; then
    echo "Error: conda.sh not found at ${CONDA_ROOT}/etc/profile.d/conda.sh"
    exit 1
fi

echo "--- PyStar Launcher ---"
echo "Reading config: $CONFIG_FILE"
echo "Run mode: $RUN_MODE"
echo "Using conda env: $CONDA_ENV_NAME"
echo "Repo root: $REPO_ROOT"

module purge
module load matlab/2023a

source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"

NUM_JOBS=$(python -c "
import yaml
try:
    with open('$CONFIG_FILE', encoding='utf-8') as f:
        data = yaml.safe_load(f)
    fovs = data['dataset']['fov_list']
    if isinstance(fovs, str) and '-' in fovs:
        start, end = map(int, fovs.split('-'))
        print(end - start + 1)
    elif isinstance(fovs, str) and ',' in fovs:
        print(len([x for x in fovs.split(',') if x.strip()]))
    elif isinstance(fovs, list):
        print(len(fovs))
    else:
        print(1)
except Exception:
    print(0)
")

if [ "$NUM_JOBS" -eq "0" ]; then
    echo "Error: Failed to parse fov_list from yaml."
    exit 1
fi

echo "Detected $NUM_JOBS FOVs to process."

mkdir -p logs/pystar

if [ "$RUN_MODE" = "main" ]; then
    CPUS_PER_FOV=8
    Batch_FOV=127
else
    CPUS_PER_FOV=4
    Batch_FOV=127
fi

JOB_ID=$(sbatch << EOF | awk '{print $4}'
#!/bin/bash
#SBATCH -J pystar_${RUN_MODE}
#SBATCH -o logs/pystar/%x.%A_%a.out
#SBATCH -e logs/pystar/%x.%A_%a.err
#SBATCH -p C64M512G
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS_PER_FOV}
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --array=1-${NUM_JOBS}%${Batch_FOV}
#SBATCH --no-requeue
#SBATCH --export=ALL

echo "Running on node: \$(hostname)"
echo "Slurm Task ID: \$SLURM_ARRAY_TASK_ID"
echo "Run mode: ${RUN_MODE}"
echo "Config file: ${CONFIG_FILE}"

cd "${REPO_ROOT}"

source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"

export OMP_NUM_THREADS=${CPUS_PER_FOV}
export MKL_NUM_THREADS=${CPUS_PER_FOV}
export OPENBLAS_NUM_THREADS=${CPUS_PER_FOV}
export PYTHONPATH="${REPO_ROOT}:\${PYTHONPATH:-}"

python "${REPO_ROOT}/01_decode/batch_pystar.py" --config "$CONFIG_FILE" --task_id "\$SLURM_ARRAY_TASK_ID" --mode "${RUN_MODE}"

EOF
)

echo "Job submitted! ID: $JOB_ID"
echo "Monitor with: squeue -j $JOB_ID"
echo "Stage logs:"
echo "  main/decode: <export_directory>/PositionXXX/decode.log"
echo "  IF:          <export_directory>/PositionXXX/if.log"
