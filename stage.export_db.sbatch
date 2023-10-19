#!/bin/bash -l
#SBATCH -J stage
#SBATCH --time 24:00:00
#SBATCH --mem 8G
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --output %x-%A.out
#SBATCH --error %x-%A.err

function logger () { 
    echo "$(date) :: ${@}"; 
    echo "$(date) :: ${@}" 1>&2; 
}

set -x 
set -e
set -o pipefail

echo $TMPDIR


if [[ "${DB_FNA}" == *".fna.xz" ]]; then
    name=$(basename ${DB_FNA} .xz)
    xzcat ${DB_FNA} > ${DB_FNA_STAGED}
else
    cp ${DB_FNA} ${DB_FNA_STAGED}
fi
