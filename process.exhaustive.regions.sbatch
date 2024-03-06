#!/bin/bash -l

#SBATCH -J contig_seqs
#SBATCH --time=100:00:00        # Walltime
#SBATCH --mem 120g     
#SBATCH -c 64       # number of cores per task
#SBATCH -N 1               # number of nodes
#SBATCH --output %x-%A.out
#SBATCH --error %x-%A.err

mamba activate human-depletion

if [[ -z ${TMPDIR} ]]; then
    echo "Must set TMPDIR"
    exit 1
fi

set -x 
set -e
set -o pipefail


#-----------------------------

#Create bed path based on file input
export LC_ALL=C
xzcat -T 7 ${EX_OUT}/*.sam.xz | \
    awk -F'\t' -v OFS='\t' '{sum = $4 + length($10); print $3, $4, sum, $1, $10}' | \
    xz -1 -z -c -T 8 > ${EX_OUT}/awk.aggregated.xz

xzcat -T 4 ${EX_OUT}/awk.aggregated.xz | \
    sort \
        --parallel 8 \
        --buffer-size=100g \
        -k1,1 \
        -k2,2n | \
        python deduplicate.py | \
        xz -1 -z -c -T 4 > ${EX_OUT}/awk.aggregated.sorted.xz

xzcat -T 15 ${EX_OUT}/awk.aggregated.sorted.xz | \
    bedtools merge \
        -i stdin \
        -c 5 \
        -o count > ${EX_OUT}/exhaustive.bed


OUTPUT=${EX_OUT}/exhaustive.fna


python trim-to-max-length.py \
    ${DB_FNA_STAGED} \
    ${EX_OUT}/exhaustive.bed \
    ${EX_OUT}/exhaustive.trimmed.bed

bedtools getfasta \
    -fi ${DB_FNA_STAGED} \
    -bed ${EX_OUT}/exhaustive.trimmed.bed \
    -fo ${OUTPUT}
