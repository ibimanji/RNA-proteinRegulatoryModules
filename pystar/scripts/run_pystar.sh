#!/bin/bash

# =================================================================
# 自动读取 Config，计算 Array 长度，提交任务
# 用法: bash scripts/run_pystar.sh [config_path]
# =================================================================

set -euo pipefail

CONFIG_FILE="${1:-config/experiment_config.yaml}"

# 你的 conda 安装位置
CONDA_ROOT="$HOME/software/miniconda3"
CONDA_ENV_NAME="pystar310"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file not found at $CONFIG_FILE"
    exit 1
fi

if [ ! -f "${CONDA_ROOT}/etc/profile.d/conda.sh" ]; then
    echo "Error: conda.sh not found at ${CONDA_ROOT}/etc/profile.d/conda.sh"
    exit 1
fi

echo "--- PyStar Launcher ---"
echo "Reading config: $CONFIG_FILE"
echo "Using conda env: $CONDA_ENV_NAME"

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
#mkdir -p logs/pystar

CPUS_PER_FOV=8
Batch_FOV=127

JOB_ID=$(sbatch << EOF | awk '{print $4}'
#!/bin/bash
#SBATCH -J pystar_batch
#SBATCH -o /dev/null
#SBATCH -e /dev/null
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

TASK_LOG_FILE=\$(python -c "
from pathlib import Path
import yaml

config_file = '$CONFIG_FILE'
task_id = int('\$SLURM_ARRAY_TASK_ID')

with open(config_file, encoding='utf-8') as f:
    data = yaml.safe_load(f)

fovs = data['dataset']['fov_list']
if isinstance(fovs, str) and '-' in fovs:
    start, end = map(int, fovs.split('-'))
    fov_list = list(range(start, end + 1))
elif isinstance(fovs, str) and ',' in fovs:
    fov_list = [int(x.strip()) for x in fovs.split(',') if x.strip()]
elif isinstance(fovs, list):
    fov_list = [int(x) for x in fovs]
else:
    fov_list = [int(fovs)]

output_cfg = data['pipeline']['output']
export_base = output_cfg.get('export_directory') or output_cfg['directory']
digits = int(output_cfg.get('export_fov_digits', 3))
log_name = output_cfg.get('export_log_name', 'log.out')

fov_id = fov_list[task_id - 1]
log_dir = Path(export_base) / f'Position{fov_id:0{digits}d}'
log_dir.mkdir(parents=True, exist_ok=True)
print(log_dir / log_name)
")

exec > "\$TASK_LOG_FILE" 2>&1

echo "Running on node: \$(hostname)"
echo "Slurm Task ID: \$SLURM_ARRAY_TASK_ID"

source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"

export OMP_NUM_THREADS=${CPUS_PER_FOV}
export MKL_NUM_THREADS=${CPUS_PER_FOV}
export OPENBLAS_NUM_THREADS=${CPUS_PER_FOV}

python scripts/batch_pystar.py --config "$CONFIG_FILE" --task_id "\$SLURM_ARRAY_TASK_ID"

EOF
)

echo "Job submitted! ID: $JOB_ID"
echo "Monitor with: squeue -j $JOB_ID"
