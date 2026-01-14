#!/bin/bash
#SBATCH --job-name=ReadCounts_Parallel
#SBATCH --output=%j_ReadCounts_Parallel.out
#SBATCH --error=%j_ReadCounts_Parallel.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --time=80:00:00

export PATH="your_conda_path/bin:$PATH"

# 创建临时文件来跟踪处理状态
TEMP_DIR=$(mktemp -d)
FAILED_TASKS="${TEMP_DIR}/failed_tasks.txt"
COMPLETED_TASKS="${TEMP_DIR}/completed_tasks.txt"
touch "${FAILED_TASKS}" "${COMPLETED_TASKS}"

# 激活环境
source activate hmmcopy
if [ $? -ne 0 ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Failed to activate conda environment"  >&2
    exit 1
fi

# 检查 parallel
if ! command -v parallel &> /dev/null; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] GNU Parallel not found. Installing via conda..."
    conda install -y -c conda-forge parallel
    if [ $? -ne 0 ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Failed to install GNU Parallel"  >&2
        exit 1
    fi
fi

# 设置变量
READ_COUNTER="your_Read_Counter_path" # hmmcopy_utils/bin/readCounter工具
BAM_DIR="your_bam_dir"
OUTPUT_DIR="your_output_dir"
WINDOW_SIZES=(your_window_size)

mkdir -p "${OUTPUT_DIR}"

# === 定义标准染色体列表 ===
# 这一步非常重要，它保证了 ReadCounts WIG 的行数和顺序与 GC/Map WIG 完全一致
# 从而避免了后续 R 脚本中 wigsToRangedData 报错
export CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX"

# 修改获取BAM文件的方式，使用find命令递归查找所有子目录中的BAM文件
mapfile -t SAMPLE_BAMS < <(find -L "${BAM_DIR}" -type f -name "*.bam")
TOTAL_SAMPLES=${#SAMPLE_BAMS[@]}

if [ ${TOTAL_SAMPLES} -eq 0 ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: No BAM files found in ${BAM_DIR}"  >&2
    exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found ${TOTAL_SAMPLES} BAM files"

# 改进的处理函数
process_readcounts() {
    local BAM_FILE=$1
    local WINDOW=$2
    local RETRY_COUNT=3
    local retry=0

    # 1. 样本ID清洗：确保只提取核心ID，方便后续匹配
    # 假设 BAM 名为 CL100164113_L01_1.sorted.rmdup...bam -> 提取 CL100164113_L01_1
    SAMPLE_ID=$(basename "${BAM_FILE}" | sed -E 's/(\.sorted)?(\.rmdup)?(\.realign)?(\.BQSR)?\.bam$//')

    WINDOW_DIR="${OUTPUT_DIR}/window_${WINDOW}"
    mkdir -p "${WINDOW_DIR}"

    # 输出文件名：带上 _20000bp 后缀，与 R 脚本预期一致
    READCOUNTER_OUTPUT="${WINDOW_DIR}/${SAMPLE_ID}_${WINDOW}bp.readcounts.wig"
    TASK_ID="${SAMPLE_ID}_${WINDOW}"

    if [ -f "${READCOUNTER_OUTPUT}" ] && [ -s "${READCOUNTER_OUTPUT}" ]; then
        echo "[SKIP] ${TASK_ID} exists."
        echo "${TASK_ID}" >> "${COMPLETED_TASKS}"
        return 0
    fi

    while [ $retry -lt $RETRY_COUNT ]; do
        # === 关键修改 ===
        # 1. 增加 -c "${CHROMS}" 参数：只输出指定染色体，顺序固定
        # 2. 移除 sed '/^fixedStep chrom=chrM/,$d'：
        #    因为我们已经用 -c 指定了 chr1..chrY (不含M)，工具本身就不会输出 chrM，
        #    这样比 sed 更安全且高效。

        ${READ_COUNTER} -w "${WINDOW}" -q 20 -c "${CHROMS}" "${BAM_FILE}" > "${READCOUNTER_OUTPUT}.tmp"

        if [ $? -eq 0 ] && [ -s "${READCOUNTER_OUTPUT}.tmp" ]; then
            mv "${READCOUNTER_OUTPUT}.tmp" "${READCOUNTER_OUTPUT}"
            echo "[DONE] ${TASK_ID}"
            echo "${TASK_ID}" >> "${COMPLETED_TASKS}"
            return 0
        else
            echo "[FAIL] Attempt $((retry+1)) for ${TASK_ID}" >&2
            rm -f "${READCOUNTER_OUTPUT}.tmp"
            retry=$((retry+1))
            sleep 1
        fi
    done

    echo "[ERROR] Final fail for ${TASK_ID}" >&2
    echo "${TASK_ID}" >> "${FAILED_TASKS}"
    return 1
}

export -f process_readcounts
# 务必导出 CHROMS 变量供函数内部使用
export READ_COUNTER OUTPUT_DIR COMPLETED_TASKS FAILED_TASKS CHROMS

# 使用parallel执行任务，并设置超时 (自行修改timeout时间和并行jobs数量)
parallel --timeout 3600 --jobs 5 process_readcounts {} ::: "${SAMPLE_BAMS[@]}" ::: "${WINDOW_SIZES[@]}"  # --timeout 3600 表示设置超时时间为3600秒  # --jobs 5 表示同时运行5个任务

# 验证处理结果
TOTAL_EXPECTED=$((TOTAL_SAMPLES * ${#WINDOW_SIZES[@]}))
COMPLETED_COUNT=$(wc -l < "${COMPLETED_TASKS}")
FAILED_COUNT=$(wc -l < "${FAILED_TASKS}")

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing summary:"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Expected tasks: ${TOTAL_EXPECTED}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed tasks: ${COMPLETED_COUNT}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Failed tasks: ${FAILED_COUNT}"

if [ -s "${FAILED_TASKS}" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Failed tasks list:"
    cat "${FAILED_TASKS}" | while read -r task; do
        echo "  - ${task}"
    done
fi

# 检查是否所有任务都已完成
if [ "$((COMPLETED_COUNT + FAILED_COUNT))" -ne "${TOTAL_EXPECTED}" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Some tasks were not processed. Please check the logs."  >&2
    exit 1
fi

# 清理临时文件
rm -rf "${TEMP_DIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Process completed"