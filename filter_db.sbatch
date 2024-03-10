#!/bin/bash -l
#SBATCH -J filtered-bt2
#SBATCH --time 120:00:00
#SBATCH --mem 1200gb
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --partition highmem
#SBATCH --output %x-%A.out
#SBATCH --error %x-%A.err

mamba activate qiime2-2023.5

set -x 
set -e

echo $TMPDIR

if [[ -z ${FILE} ]]; then
    FILE=${BASE_LABEL}.fna
fi
mamba activate bowtie2

python exhaustive.py filter-db \
    --database-fasta ${DB_FNA_STAGED} \
    --contaminated-fasta ${FILE} \
    --output-fasta ${DB_FNA_STAGED}.masked \
    --output-bt2 ${BASE_LABEL}-bt2/${BASE_LABEL}-bt2 \
    --threads 32

# python filter_fna.py \
#     --database-fasta ${DB_FNA_STAGED} \
#     --output-fasta ${DB_FNA_STAGED}.masked \
#     --contaminated-fasta ${FILE}

# mkdir -p ${BASE_LABEL}-bt2
# bowtie2-build --version
# bowtie2-build \
#     --seed 42 \
#     --threads 32 ${DB_FNA_STAGED}.masked \
#     ${BASE_LABEL}-bt2/${BASE_LABEL}-bt2
