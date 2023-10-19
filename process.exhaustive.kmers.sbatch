#!/bin/bash -l
#SBATCH -J exhaustive
#SBATCH --time 120:00:00
#SBATCH -N 1
#SBATCH --exclude=brncl-[23-31]
#SBATCH --array 0-15
#SBATCH --output %x-%A_%a.out
#SBATCH --error %x-%A_%a.err

mamba activate bowtie2

set -x 
set -e
set -o pipefail


if [[ -z ${FILES} ]]; then
    echo "FILES not set"
    exit 1
fi

if [[ -z ${EX_OUT} ]]; then
    echo "FILES not set"
    exit 1
fi

if [[ -z ${DB_BT2} ]]; then
    echo "DB_BT2 not set"
    exit 1
fi

if [[ ! -f ${FILES} ]]; then
    echo "FILES not found"
    exit 1
fi

mkdir -p ${EX_OUT}

bt2_threads=$((${SLURM_CPUS_PER_TASK} - 2))
for f in $(cat ${FILES})
do
    zcat ${f}
done | \
    python kmer.py ${K} ${J} ${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_TASK_COUNT} | \
    bowtie2 -p ${bt2_threads} \
        -x ${DB_BT2} \
        -f - \
        --seed 42 \
        --very-sensitive \
        -k 32 \
        --np 1 \
        --mp "1,1" \
        --rdg "0,1" \
        --rfg "0,1" \
        --score-min "L,0,-0.05" \
        --no-head \
        --no-unal | \
            xz -1 -z -c > ${EX_OUT}/${SLURM_ARRAY_TASK_ID}.sam.xz
