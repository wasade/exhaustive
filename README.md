Exhaustive
----------

Perform exhaustive short read mapping against a reference database. 
Specifically, we generate all successive `k` length substrings, jumping `j`
characters, across the set of input sequences. These substrings are then
mapped against a database to determine regions of high identity. Mapped 
coordinates are transformed into contiguous regions, and a fasta file
containing the observed sequences is produced. Finally, the database is masked
such that regions of similarity are converted to N characters, and a new
bowtie2 index is then constructed.

Install
-------

WYSIWYG for the time being. The scripts are designed around a SLURM submission
system. The included Python scripts should generalize across Python versions
with the only dependencies being Click and NumPy. 

Running
-------

`submit-exhaustive.sh` requires a label, the database as fasta or xz compressed
fasta, a bowtie2 index to map against, and a file containing the paths to the
sequence files to generate short reads from. The submission will issue calls
to `sbatch` and setup dependencies for subsequent work.
