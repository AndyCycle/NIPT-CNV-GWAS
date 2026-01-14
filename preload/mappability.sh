#!/bin/bash
#SBATCH --job-name=HMMCopy_Mappability
#SBATCH --output=HMMCopy_Mappability_%j.out
#SBATCH --error=HMMCopy_Mappability_%j.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10G
#SBATCH --time=00:30:00

echo "Mappability processing started at: $(date)"

export PATH="your_conda_path/bin:$PATH"

# 加载必要的模块或激活环境
source activate hmmcopy  # 激活 Conda 环境

# 设置变量
MAP_COUNTER=mapCounter

# 输入文件路径
MAPPABILITY_BW="k100.Umap.MultiTrackMappability.bw"

# 输出目录
OUTPUT_DIR="preload/mappability"  # 替换为实际的输出目录
mkdir -p "${OUTPUT_DIR}"

# === 定义标准染色体列表 ===
# 必须不带空格，逗号分隔
CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX"

# 定义窗口大小数组（以 bp 为单位）
WINDOW_SIZES=(50000 20000 10000 5000)

# 循环处理每个窗口大小
for WINDOW in "${WINDOW_SIZES[@]}"; do
    echo "Processing mappability for window size: ${WINDOW} bp"
    WINDOW_DIR="${OUTPUT_DIR}/window_${WINDOW}"
    mkdir -p "${WINDOW_DIR}"
    MAPPABILITY_OUTPUT="${WINDOW_DIR}/k100.Umap.MultiTrackMappability_${WINDOW}bp.wig"

    echo "Running mapCounter with window size ${WINDOW}..."

    # === 关键修改：增加 -c 参数 ===
    ${MAP_COUNTER} -w "${WINDOW}" -c "${CHROMS}" "${MAPPABILITY_BW}" > "${MAPPABILITY_OUTPUT}"

    if [ $? -eq 0 ]; then
        echo "Success: ${MAPPABILITY_OUTPUT}"
    else
        echo "Error in mapCounter" >&2
        exit 1
    fi
done

echo "All mappability tasks completed successfully."
echo "Process completed at: $(date)"
