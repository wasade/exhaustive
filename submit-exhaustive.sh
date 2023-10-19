#!/bin/bash

if [[ $# -eq 4 ]]; then
    BASE_LABEL=${1}  # a label for the run
    DB_FNA=${2}  # the database as fasta
    DB_BT2=${3}  # the database to map against
    FILES=${4}  # genomes to generate substrings from
else
    echo "bad arguments"
    exit 1
fi

export TMPDIR=$(mktemp -d)
echo ${TMPDIR}

extra_args=""

K=150
JUMP=75
EX_OUT=${PANFS}/${BASE_LABEL}-exhaustive-${K}-${JUMP}
DB_FNA_STAGED=${TMPDIR}/concat.fna

jd=$(sbatch \
     -J ${BASE_LABEL}-stage \
     --parsable \
     --mem 8G \
     -N 1 \
     -c 1 \
     ${extra_args} \
     --export DB_FNA=${DB_FNA},DB_FNA_STAGED=${DB_FNA_STAGED},TMPDIR=${TMPDIR} \
     stage.export_db.sbatch)

# exhaustive
jex=$(sbatch \
      -J ${BASE_LABEL}-ex \
      --mem 100G \
      -N 1 \
      -c 16 \
      ${extra_args} \
      --parsable \
      --export EX_OUT=${EX_OUT},DB_BT2=${DB_BT2},K=${K},J=${JUMP} \
      process.exhaustive.kmers.sbatch)

jex=$(sbatch \
      --parsable \
      -J ${BASE_LABEL}-exr \
      --export EX_OUT=${EX_OUT},DB_FNA_STAGED=${DB_FNA_STAGED},TMPDIR=${TMPDIR} \
      ${extra_args} \
      --dependency=afterok:${jex}:${jd} \
      process.exhaustive.regions.sbatch)

sbatch \
    -J ${BASE_LABEL}-bt2 \
    --dependency=afterok:${jex} \
    --export TMPDIR=${TMPDIR},DB_FNA_STAGED=${DB_FNA_STAGED},BASE_LABEL=${BASE_LABEL},FILE=${EX_OUT}/exhaustive.fna \
    filter_db.sbatch
