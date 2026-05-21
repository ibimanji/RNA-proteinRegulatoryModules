#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# run_cm_repair_missing_tiles_and_stitch.sh
#
# 作用：
#   1) 扫描每个 POS×param 的 WORK_DIR，读取 manifest.json 的 n_tiles
#   2) 检查 tile_results 下缺失/损坏的 tile_XXXXX.pkl（不存在/空文件/joblib load失败/None）
#   3) 只对缺的 tile 用 sbatch array 补提交（patch tiles）
#   4) patch tiles 完成后提交 stitch（新 job，不再依赖旧 rescue 链，杜绝 NeverSatisfied）
#
# 用法：
#   bash run_cm_repair_missing_tiles_and_stitch.sh [PARAM_TSV] [POS_SPEC] [--dry] [--wait] [--stitch-only]
#   bash run_cm_repair_missing_tiles_and_stitch.sh params.tsv "Position001-Position127"
#   bash run_cm_repair_missing_tiles_and_stitch.sh params.tsv "Position010,Position040,Position117" --dry
#   --dry         只打印会提交什么，不真的 sbatch
#   --wait        在脚本里等待 patch tiles 结束后再 stitch（更稳，但 login 上会挂着等）
#   --stitch-only 不补 tile，只对未 done 的直接提交 stitch（适合你确认 tile 已齐）
# =========================================================

PY="/gpfs/share/home/2301920002/software/miniconda3/envs/clustermaptest/bin/python"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT_DIR}/new_clustermap_fast.py"

PARAM_TSV="${1:-${SCRIPT_DIR}/params.tsv}"
POS_SPEC="${2:-Position010,Position040,Position117}"

BASE="/gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mWT"
Z_NUM=42
XY_SIZE=3072
GENE_CSV="${BASE}/01_data/genes.csv"
IDN="4"
ROUND_DAPI="IF"
DAPI_ROUND_NUM="1"

LOG_ROOT="logs"
WORK_ROOT_BASE="${BASE}/03_segmentation/clustermap_tmp"
OUT_ROOT_BASE="${BASE}/03_segmentation"

ARRAY_CONCURRENCY=60
PARTITION="C64M512G"
QOS="normal"

# tiles patch 资源（用你原来的）
TILE_CPUS=8
TILE_MEM="64G"
TILE_TIME="24:00:00"
PATCH_CONCURRENCY=12

# stitch 资源（用你原来的）
STITCH_CPUS=16
STITCH_MEM="192G"
STITCH_TIME="12:00:00"

# sbatch 重试（用你原来的）
SBATCH_RETRIES="${SBATCH_RETRIES:-6}"
SBATCH_RETRY_SLEEP="${SBATCH_RETRY_SLEEP:-2}"

mkdir -p "$LOG_ROOT"
RUN_ID="$(date +%Y%m%d_%H%M%S)_$$"

log(){ echo "$@" >&2; }
need_cmd(){ command -v "$1" >/dev/null 2>&1 || { log "ERROR: missing command: $1"; exit 1; }; }

need_cmd sbatch
need_cmd squeue
need_cmd awk
need_cmd sed
need_cmd date

[[ -x "$PY" ]] || { log "ERROR: PY not executable: $PY"; exit 1; }
[[ -f "$SCRIPT" ]] || { log "ERROR: SCRIPT not found: $SCRIPT"; exit 1; }
[[ -f "$PARAM_TSV" ]] || { log "ERROR: PARAM_TSV not found: $PARAM_TSV"; exit 1; }

# ---------------- args flags ----------------
DRY_RUN=0
WAIT_PATCH=0
STITCH_ONLY=0
for a in "${@:3}"; do
  case "$a" in
    --dry) DRY_RUN=1 ;;
    --wait) WAIT_PATCH=1 ;;
    --stitch-only) STITCH_ONLY=1 ;;
    *) log "[WARN] unknown arg ignored: $a" ;;
  esac
done

log "[INFO] PY=$PY"
log "[INFO] SCRIPT=$SCRIPT"
log "[INFO] PARAM_TSV=$PARAM_TSV"
log "[INFO] POS_SPEC=$POS_SPEC"
log "[INFO] LOG_ROOT=$LOG_ROOT"
log "[INFO] DRY_RUN=$DRY_RUN WAIT_PATCH=$WAIT_PATCH STITCH_ONLY=$STITCH_ONLY"
log "[INFO] SBATCH_RETRIES=$SBATCH_RETRIES SBATCH_RETRY_SLEEP=$SBATCH_RETRY_SLEEP"

# ---------- 原子写文件 ----------
write_atomic() {
  local dst="${1:?write_atomic: missing destination path}"
  local tmp
  tmp="$(mktemp "${dst}.tmp.XXXXXX")"
  cat > "$tmp"
  chmod +x "$tmp"
  mv -f "$tmp" "$dst"
}

# ---------- sbatch 重试 ----------
sbatch_retry() {
  local attempt=1
  local sleep_s="${SBATCH_RETRY_SLEEP}"
  local out
  while (( attempt <= SBATCH_RETRIES )); do
    set +e
    out="$(sbatch --parsable "$@" 2>&1)"
    local rc=$?
    set -e
    if (( rc == 0 )) && [[ "$out" =~ ^[0-9]+ ]]; then
      echo "$out"
      return 0
    fi
    log "[WARN] sbatch failed (attempt ${attempt}/${SBATCH_RETRIES}, rc=${rc})"
    log "[WARN] sbatch output: ${out}"
    if (( attempt == SBATCH_RETRIES )); then
      log "[ERROR] sbatch still failing after retries."
      return 99
    fi
    log "[RETRY] sleep ${sleep_s}s then retry..."
    sleep "${sleep_s}"
    sleep_s=$(( sleep_s * 2 ))
    attempt=$(( attempt + 1 ))
  done
  return 99
}

# ---------- Position 展开 ----------
expand_pos_spec() {
  local spec="$1"
  local tok
  spec="${spec//,/ }"
  for tok in $spec; do
    if [[ "$tok" =~ ^[Pp]osition([0-9]+)-([Pp]osition)?([0-9]+)$ ]]; then
      local a=$((10#${BASH_REMATCH[1]})) b=$((10#${BASH_REMATCH[3]}))
      if (( a <= b )); then for ((i=a;i<=b;i++)); do printf "Position%03d\n" "$i"; done
      else for ((i=a;i>=b;i--)); do printf "Position%03d\n" "$i"; done; fi
    elif [[ "$tok" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      local a=$((10#${BASH_REMATCH[1]})) b=$((10#${BASH_REMATCH[2]}))
      if (( a <= b )); then for ((i=a;i<=b;i++)); do printf "Position%03d\n" "$i"; done
      else for ((i=a;i>=b;i--)); do printf "Position%03d\n" "$i"; done; fi
    else
      if [[ "$tok" =~ ^[Pp]osition[0-9]+$ ]]; then
        local n="${tok#[Pp]osition}"
        printf "Position%03d\n" "$((10#$n))"
      else
        log "ERROR: cannot parse position token: $tok"
        exit 2
      fi
    fi
  done
}
mapfile -t POS_LIST < <(expand_pos_spec "$POS_SPEC" | awk '!seen[$0]++')
log "[INFO] Expanded positions: ${#POS_LIST[@]} -> ${POS_LIST[*]}"

# ---------- params TSV ----------
mapfile -t PARAM_LINES < <(awk 'BEGIN{FS="[ \t]+"} !/^#/ && NF>=6 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' "$PARAM_TSV")
log "[INFO] Params lines: ${#PARAM_LINES[@]}"
(( ${#PARAM_LINES[@]} > 0 )) || { log "ERROR: PARAM_TSV has no valid lines (need: tag IC ID IPF ICR_xy ICR_z)."; exit 3; }

# ---------- done 判定（沿用你原来的严格标准） ----------
is_done_strict() {
  local out_dir="$1"
  local work_dir="$2"
  [[ -s "${out_dir}/model_final.pkl" ]] || return 1
  [[ -s "${out_dir}/remain_reads_raw.csv" || -s "${out_dir}/remain_reads.csv" ]] || return 1
  if [[ -e "${out_dir}/cell_center.csv" && ! -s "${out_dir}/cell_center.csv" ]]; then
    return 1
  fi
  if compgen -G "${work_dir}/tile_results/tile_*.pkl" > /dev/null; then
    local newest_tile
    newest_tile="$(ls -1t "${work_dir}/tile_results"/tile_*.pkl 2>/dev/null | head -n 1 || true)"
    if [[ -n "${newest_tile}" ]]; then
      if [[ "${out_dir}/model_final.pkl" -ot "${newest_tile}" ]]; then
        return 1
      fi
    fi
  fi
  return 0
}

# ---------- 写模板（patch tiles + stitch） ----------
TPL_DIR="${LOG_ROOT}/_repair_templates_${RUN_ID}"
mkdir -p "$TPL_DIR"

PATCH_TILES_SBATCH="${TPL_DIR}/sbatch_patch_tiles.sh"
PATCH_STITCH_SBATCH="${TPL_DIR}/sbatch_patch_stitch.sh"

# ---- PATCH TILES ----
write_atomic "$PATCH_TILES_SBATCH" <<EOF
#!/bin/bash
#SBATCH -p ${PARTITION}
#SBATCH --qos=${QOS}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${TILE_CPUS}
#SBATCH --mem=${TILE_MEM}
#SBATCH --time=${TILE_TIME}
#SBATCH -o ${LOG_ROOT}/%x.%A_%a.out
#SBATCH -e ${LOG_ROOT}/%x.%A_%a.err
set -euo pipefail

ts(){ date +"%F %T"; }
echo "[\$(ts)] [BEGIN] stage=patch_tile jobid=\${SLURM_JOB_ID} array=\${SLURM_ARRAY_JOB_ID}.\${SLURM_ARRAY_TASK_ID} node=\${SLURMD_NODENAME:-NA} POS=\${POS} TAG=\${TAG}" >&2

: "\${PY:?}" "\${SCRIPT:?}" "\${BASE:?}" "\${POS:?}" "\${TAG:?}" "\${Z_NUM:?}" "\${XY_SIZE:?}" "\${GENE_CSV:?}" "\${DAPI_DIR:?}" "\${GOOD_POINTS:?}" "\${IDN:?}" "\${WORK_DIR:?}" "\${OUT_DIR:?}"
: "\${IC:?}" "\${ID:?}" "\${IPF:?}" "\${ICR_XY:?}" "\${ICR_Z:?}"
: "\${SHARED_PREP_DIR:?}" "\${DAPI_ROUND_NUM:?}"

# tiles link (必须存在)
if [[ ! -e "\${WORK_DIR}/tiles" ]]; then
  ln -sfn "\${SHARED_PREP_DIR}/tiles" "\${WORK_DIR}/tiles"
fi

MANIFEST="\${WORK_DIR}/manifest.json"
[[ -f "\${MANIFEST}" ]] || { echo "ERROR: manifest missing: \${MANIFEST}" >&2; exit 2; }

N_TILES=\$("\${PY}" - <<PY
import json
print(json.load(open("\${MANIFEST}"))["n_tiles"])
PY
)

TID=\${SLURM_ARRAY_TASK_ID}
[[ "\${TID}" -lt "\${N_TILES}" ]] || { echo "[\$(ts)] [END] stage=patch_tile status=SKIP_out_of_range tid=\${TID} n=\${N_TILES}" >&2; exit 0; }

thr=\${SLURM_CPUS_PER_TASK:-1}
(( thr > 8 )) && thr=8
export OMP_NUM_THREADS=\$thr MKL_NUM_THREADS=\$thr OPENBLAS_NUM_THREADS=\$thr NUMEXPR_NUM_THREADS=\$thr

OUTPKL="\${WORK_DIR}/tile_results/tile_\$(printf "%05d" "\${TID}").pkl"

ICR="\${ICR_XY},\${ICR_Z}"

"\${PY}" "\${SCRIPT}" \\
  -IP "\${POS}" -IZ "\${Z_NUM}" -IEP F \\
  -IDir "\${BASE}" \\
  -Igenecsv "\${GENE_CSV}" \\
  -IDapi_path "\${DAPI_DIR}" \\
  -IXY "\${XY_SIZE}" -IC "\${IC}" -ID "\${ID}" -IDR "\${DAPI_ROUND_NUM}" -IPF "\${IPF}" -ICR "\${ICR}" \\
  -OP "\${OUT_DIR}" \\
  -Igood_points_max3d "\${GOOD_POINTS}" \\
  -IDN "\${IDN}" \\
  --stage tile \\
  --work_dir "\${WORK_DIR}" \\
  --tile_id "\${TID}"

# sanity check
if [[ ! -s "\${OUTPKL}" ]]; then
  echo "[\$(ts)] [END] stage=patch_tile status=FAIL no_output tid=\${TID}" >&2
  exit 5
fi

"\${PY}" - <<PY
import joblib, sys
p=r"""\${OUTPKL}"""
try:
    obj=joblib.load(p)
    sys.exit(0 if obj is not None else 6)
except Exception:
    sys.exit(7)
PY

echo "[\$(ts)] [END] stage=patch_tile status=OK tid=\${TID}" >&2
EOF

# ---- PATCH STITCH ----
write_atomic "$PATCH_STITCH_SBATCH" <<EOF
#!/bin/bash
#SBATCH -p ${PARTITION}
#SBATCH --qos=${QOS}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${STITCH_CPUS}
#SBATCH --mem=${STITCH_MEM}
#SBATCH --time=${STITCH_TIME}
#SBATCH -o ${LOG_ROOT}/%x.%j.out
#SBATCH -e ${LOG_ROOT}/%x.%j.err
set -euo pipefail

ts(){ date +"%F %T"; }
echo "[\$(ts)] [BEGIN] stage=patch_stitch jobid=\${SLURM_JOB_ID} node=\${SLURMD_NODENAME:-NA} POS=\${POS} TAG=\${TAG}" >&2

: "\${PY:?}" "\${SCRIPT:?}" "\${BASE:?}" "\${POS:?}" "\${TAG:?}" "\${Z_NUM:?}" "\${XY_SIZE:?}" "\${GENE_CSV:?}" "\${DAPI_DIR:?}" "\${GOOD_POINTS:?}" "\${IDN:?}" "\${WORK_DIR:?}" "\${OUT_DIR:?}"
: "\${IC:?}" "\${ID:?}" "\${IPF:?}" "\${ICR_XY:?}" "\${ICR_Z:?}" "\${DAPI_ROUND_NUM:?}"

thr=\${SLURM_CPUS_PER_TASK:-1}
(( thr > 16 )) && thr=16
export OMP_NUM_THREADS=\$thr MKL_NUM_THREADS=\$thr OPENBLAS_NUM_THREADS=\$thr NUMEXPR_NUM_THREADS=\$thr

ICR="\${ICR_XY},\${ICR_Z}"

"\${PY}" "\${SCRIPT}" \\
  -IP "\${POS}" -IZ "\${Z_NUM}" -IEP F \\
  -IDir "\${BASE}" \\
  -Igenecsv "\${GENE_CSV}" \\
  -IDapi_path "\${DAPI_DIR}" \\
  -IXY "\${XY_SIZE}" -IC "\${IC}" -ID "\${ID}" -IDR "\${DAPI_ROUND_NUM}" -IPF "\${IPF}" -ICR "\${ICR}" \\
  -OP "\${OUT_DIR}" \\
  -Igood_points_max3d "\${GOOD_POINTS}" \\
  -IDN "\${IDN}" \\
  --stage stitch \\
  --work_dir "\${WORK_DIR}"

touch "\${OUT_DIR}/.cm_done"
echo "[\$(ts)] [END] stage=patch_stitch status=OK" >&2
EOF

# ---------- helper: 计算缺失 tiles 列表 ----------
missing_tiles_csv() {
  local work_dir="$1"
  local manifest="$work_dir/manifest.json"
  [[ -f "$manifest" ]] || { echo ""; return 0; }
  WORK_DIR="$work_dir" "$PY" - <<'PY'
import os, json, joblib
work=os.environ["WORK_DIR"]
man=os.path.join(work,"manifest.json")
try:
    n=json.load(open(man))["n_tiles"]
except Exception:
    print("")
    raise SystemExit(0)

tile_out=os.path.join(work,"tile_results")
miss=[]
for t in range(int(n)):
    p=os.path.join(tile_out, f"tile_{t:05d}.pkl")
    if (not os.path.exists(p)) or os.path.getsize(p)==0:
        miss.append(t); continue
    try:
        obj=joblib.load(p)
        if obj is None:
            miss.append(t)
    except Exception:
        miss.append(t)
print(",".join(map(str, miss)))
PY
}

# ---------- 主流程：扫 -> patch -> stitch ----------
log ""
log "[INFO] Repairing pipelines: positions(${#POS_LIST[@]}) × params(${#PARAM_LINES[@]})"
log "[INFO] Templates: $TPL_DIR"
log ""

TOTAL=0
PATCHED=0
STITCHED=0

for POS in "${POS_LIST[@]}"; do
  log "#################### POSITION: ${POS} ####################"

  local_DAPI_DIR="${BASE}/02_registration/${ROUND_DAPI}/PI/${POS}.tif"
  local_GOOD_POINTS="${BASE}/02_registration/${POS}/goodPoints_max3d.csv"
  [[ -f "$local_DAPI_DIR" ]] || { log "[WARN] skip ${POS}: DAPI missing: $local_DAPI_DIR"; continue; }
  [[ -f "$local_GOOD_POINTS" ]] || { log "[WARN] skip ${POS}: goodPoints missing: $local_GOOD_POINTS"; continue; }

  POS_ROOT_WORK="${WORK_ROOT_BASE}/${POS}"
  POS_ROOT_OUT="${OUT_ROOT_BASE}/${POS}"
  SHARED_PREP_DIR="${POS_ROOT_WORK}/_shared_prep"
  [[ -f "${SHARED_PREP_DIR}/manifest.json" ]] || { log "[WARN] skip ${POS}: shared_prep manifest missing: ${SHARED_PREP_DIR}/manifest.json"; continue; }
  [[ -d "${SHARED_PREP_DIR}/tiles" ]] || { log "[WARN] skip ${POS}: shared_prep tiles dir missing: ${SHARED_PREP_DIR}/tiles"; continue; }

  for line in "${PARAM_LINES[@]}"; do
    TOTAL=$((TOTAL+1))
    IFS=$'\t' read -r TAG IC ID IPF ICR_XY ICR_Z <<< "$line"

    WORK_DIR="${POS_ROOT_WORK}/${TAG}"
    OUT_DIR="${POS_ROOT_OUT}"
    mkdir -p "$WORK_DIR/tile_results" "$OUT_DIR"

    if is_done_strict "$OUT_DIR" "$WORK_DIR"; then
      log "[SKIP] DONE(strict) POS=${POS} TAG=${TAG}"
      continue
    fi

    # 缺 tile?
    MISS="$(missing_tiles_csv "$WORK_DIR" || true)"
    if [[ -z "$MISS" ]]; then
      log "[OK] tiles complete POS=${POS} TAG=${TAG} -> submit stitch"
      PATCH_JOBID=""
    else
      log "[MISS] POS=${POS} TAG=${TAG} missing tiles: ${MISS}"
      PATCHED=$((PATCHED+1))

      if (( STITCH_ONLY == 1 )); then
        log "[STITCH_ONLY] skip patch tiles POS=${POS} TAG=${TAG}"
        PATCH_JOBID=""
      else
        EXPORTS="ALL,PY=${PY},SCRIPT=${SCRIPT},BASE=${BASE},POS=${POS},TAG=${TAG},Z_NUM=${Z_NUM},XY_SIZE=${XY_SIZE},GENE_CSV=${GENE_CSV},DAPI_DIR=${local_DAPI_DIR},GOOD_POINTS=${local_GOOD_POINTS},IDN=${IDN},WORK_DIR=${WORK_DIR},OUT_DIR=${OUT_DIR},IC=${IC},ID=${ID},IPF=${IPF},ICR_XY=${ICR_XY},ICR_Z=${ICR_Z},SHARED_PREP_DIR=${SHARED_PREP_DIR},DAPI_ROUND_NUM=${DAPI_ROUND_NUM}"

        if (( DRY_RUN == 1 )); then
          log "[DRY] sbatch patch tiles: cm_patch_tile_${POS}_${TAG} --array ${MISS}%${PATCH_CONCURRENCY}"
          PATCH_JOBID="DRY"
        else
          PATCH_JOBID="$(sbatch_retry \
            --job-name="cm_patch_tile_${POS}_${TAG}" \
            -p "${PARTITION}" --qos="${QOS}" \
            --array="${MISS}%${PATCH_CONCURRENCY}" \
            --cpus-per-task="${TILE_CPUS}" --mem="${TILE_MEM}" --time="${TILE_TIME}" \
            --export="${EXPORTS}" \
            "${PATCH_TILES_SBATCH}")"
          log "[JOB] PATCH_TILES POS=${POS} TAG=${TAG} -> ${PATCH_JOBID}"
        fi

        if (( WAIT_PATCH == 1 )) && [[ "${PATCH_JOBID}" != "DRY" ]]; then
          log "[WAIT] waiting patch tiles job ${PATCH_JOBID} ..."
          while squeue -j "${PATCH_JOBID}" -h 2>/dev/null | grep -q .; do
            sleep 30
          done
          # 再算一遍缺失，确保补齐
          MISS2="$(missing_tiles_csv "$WORK_DIR" || true)"
          if [[ -n "$MISS2" ]]; then
            log "[WARN] still missing after patch job ${PATCH_JOBID}: ${MISS2}"
          else
            log "[OK] patch tiles completed and verified POS=${POS} TAG=${TAG}"
          fi
        fi
      fi
    fi

    # 提交 stitch（依赖 patch tiles（如果有且没 wait））
    EXPORTS2="ALL,PY=${PY},SCRIPT=${SCRIPT},BASE=${BASE},POS=${POS},TAG=${TAG},Z_NUM=${Z_NUM},XY_SIZE=${XY_SIZE},GENE_CSV=${GENE_CSV},DAPI_DIR=${local_DAPI_DIR},GOOD_POINTS=${local_GOOD_POINTS},IDN=${IDN},WORK_DIR=${WORK_DIR},OUT_DIR=${OUT_DIR},IC=${IC},ID=${ID},IPF=${IPF},ICR_XY=${ICR_XY},ICR_Z=${ICR_Z},DAPI_ROUND_NUM=${DAPI_ROUND_NUM}"

    DEP_ARGS=()
    if (( WAIT_PATCH == 0 )) && [[ -n "${PATCH_JOBID:-}" && "${PATCH_JOBID}" != "DRY" ]]; then
      DEP_ARGS+=(--dependency="afterany:${PATCH_JOBID}")
    fi

    if (( DRY_RUN == 1 )); then
      log "[DRY] sbatch stitch: cm_patch_stitch_${POS}_${TAG} (dep=${DEP_ARGS[*]:-none})"
      STITCHED=$((STITCHED+1))
      continue
    fi

    STITCH_JOBID="$(sbatch_retry \
      --job-name="cm_patch_stitch_${POS}_${TAG}" \
      -p "${PARTITION}" --qos="${QOS}" \
      --cpus-per-task="${STITCH_CPUS}" --mem="${STITCH_MEM}" --time="${STITCH_TIME}" \
      "${DEP_ARGS[@]}" \
      --export="${EXPORTS2}" \
      "${PATCH_STITCH_SBATCH}")"

    log "[JOB] STITCH POS=${POS} TAG=${TAG} -> ${STITCH_JOBID}"
    echo -e "${POS}\t${TAG}\t${IC}\t${ID}\t${IPF}\t${ICR_XY}\t${ICR_Z}\t${MISS:-}\t${STITCH_JOBID}\t${PATCH_JOBID:-}" >> "${LOG_ROOT}/repair_submitted.tsv"
    STITCHED=$((STITCHED+1))
  done

  log ""
done

log "✅ scanned combos=${TOTAL} patched=${PATCHED} stitched_submitted=${STITCHED}"
log "Track: ${LOG_ROOT}/repair_submitted.tsv"
log "Watch: squeue -u \$USER | egrep 'cm_patch_(tile|stitch)_'"
