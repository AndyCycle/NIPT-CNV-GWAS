#!/bin/sh
#SBATCH --job-name=HMMCopy_GCContent
#SBATCH --output=HMMCopy_GCContent_%j.out
#SBATCH --error=HMMCopy_GCContent_%j.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10G
#SBATCH --time=01:00:00

echo "GC Content processing started at: $(date)"

export PATH=/share/home/lsy_chenyanchao/software/miniconda3/bin:$PATH

source activate hmmcopy  # 激活 Conda 环境

# 设置变量
GC_COUNTER=gcCounter

# 输入文件路径
FASTA_FILE="/share/home/lsy_liusiyang/20220708_Alignment/hg38/Homo_sapiens_assembly38.fasta"

# 输出目录
OUTPUT_DIR="/share/home/lsy_student/chenyanchao/GDM/CNV/HMMcopy/preload/gc"  # 替换为实际的输出目录
mkdir -p "${OUTPUT_DIR}"

# === 定义标准染色体列表 ===
# 必须不带空格，逗号分隔
CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX"

# 定义窗口大小数组（以 bp 为单位）
WINDOW_SIZES=(50000 20000 10000 5000)

# 生成 FASTA 索引文件（如果未存在）
FASTA_INDEX="${FASTA_FILE}.fai"
if [ ! -f "${FASTA_INDEX}" ]; then
    echo "Generating FASTA index..."
    samtools faidx "${FASTA_FILE}"
    if [ $? -ne 0 ]; then
        echo "Error generating FASTA index." >&2
        exit 1
    fi
fi

# 循环处理每个窗口大小
for WINDOW in "${WINDOW_SIZES[@]}"; do
    echo "Processing GC content for window size: ${WINDOW} bp"
    WINDOW_DIR="${OUTPUT_DIR}/window_${WINDOW}"
    mkdir -p "${WINDOW_DIR}"
    GC_OUTPUT="${WINDOW_DIR}/Homo_sapiens_assembly38_${WINDOW}bp.gc.wig"

    echo "Running gcCounter with window size ${WINDOW}..."

    # === 关键修改：增加 -c 参数 ===
    ${GC_COUNTER} -w "${WINDOW}" -c "${CHROMS}" "${FASTA_FILE}" > "${GC_OUTPUT}"

    if [ $? -eq 0 ]; then
        echo "Success: ${GC_OUTPUT}"
    else
        echo "Error in gcCounter" >&2
        exit 1
    fi
done

echo "All GC content tasks completed successfully."
echo "Process completed at: $(date)"

