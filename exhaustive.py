#!/usr/bin/env python3

import os
import sys
import click
import subprocess
from pathlib import Path
from collections import OrderedDict

ROOT_DIR = Path(__file__).parent.resolve()

DEFAULTS = {"kmer": 150, "jump": 75, "threads": 8}


shared_options = {
    "db_fna": click.option(
        "--db-fna",
        "-d",
        help="the database in fasta format (fna, fna.gz, or fna.xz) to map against",
        required=True,
        type=click.Path(exists=True, dir_okay=False, path_type=Path),
    ),
    "db_bt2": click.option(
        "--db-bt2",
        "-b",
        help="the bowtie2 index of the database to map against",
        required=True,
        type=click.Path(path_type=Path),
    ),
    "files": click.option(
        "--files",
        "-f",
        help="the paths of genomes to generate substrings from",
        type=click.Path(exists=True, dir_okay=False, path_type=Path),
    ),
    "kmer": click.option(
        "--kmer", "-k", help="kmer length", type=int, default=DEFAULTS["kmer"]
    ),
    "jump": click.option(
        "--jump", "-j", help="jump length", type=int, default=DEFAULTS["jump"]
    ),
    "threads": click.option(
        "--threads",
        "-t",
        help="number of threads to use",
        type=int,
        default=DEFAULTS["threads"],
    ),
}


class CliOptions:
    def __init__(self, options):
        self.options = options
        for name, option in options.items():
            setattr(self, name, option)


shared_opts = CliOptions(shared_options)


def parse_base_label(base_label):
    base_label = Path(base_label)
    out_dir = base_label.parent.resolve()
    out_dir.mkdir(exist_ok=True, parents=True)
    base_label = base_label.name
    return base_label, out_dir


def _export_db(db_fna, db_fna_staged, build_index=False):

    if db_fna.suffix == ".gz":
        db_fna_staged = db_fna.with_suffix("")
        os.system(f"gunzip -c {db_fna}.gz > {db_fna_staged}")
    elif db_fna.suffix == ".xz":
        db_fna_staged = db_fna.with_suffix("")
        os.system(f"xzcat {db_fna} > {db_fna_staged}")
    elif db_fna.suffix == ".fna":
        os.system(f"cp {db_fna} {db_fna_staged}")
    else:
        raise ValueError(f"Unsupported database format: {db_fna}")

    one_line_fna = db_fna_staged.with_suffix(".one-line.fna")
    with (
        open(db_fna_staged) as f_in,
        open(db_fna_staged.with_suffix(".one-line.fna"), "w") as f_out,
    ):
        is_first = True
        for line in f_in:
            if line.startswith(">"):
                if is_first:
                    f_out.write(line)
                    is_first = False
                else:
                    f_out.write(f"\n{line}")

            else:
                f_out.write(line.strip())
    db_fna_staged = one_line_fna.rename(db_fna_staged)
    if build_index:
        db_bt2 = db_fna_staged.parent / f"{db_fna_staged.stem}-bt2"
        p_bt2 = subprocess.run(["bowtie2-build", str(db_fna_staged), str(db_bt2)])
        if p_bt2.returncode:
            raise RuntimeError("bowtie2-build failed")
    else:
        db_bt2 = None
    return db_fna_staged, db_bt2


def _process_kmers(kmer, jump, files, task_id, task_count, threads, db_bt2, kmers_sam):
    if kmers_sam.is_file():
        print(f"Skipping {kmers_sam} as it already exists")
    else:
        genomes = []
        with open(files) as fp:
            for line in fp:
                file = line.strip()
                genomes.append(file)
        p_cat = subprocess.Popen(["cat"] + genomes, stdout=subprocess.PIPE)
        p_kmer = subprocess.Popen(
            [
                "python",
                str(ROOT_DIR / "kmer.py"),
                str(kmer),
                str(jump),
                str(task_id),
                str(task_count),
            ],
            stdin=p_cat.stdout,
            stdout=subprocess.PIPE,
        )
        if p_kmer.returncode:
            raise RuntimeError("kmer.py failed")

        print("Mapping kmer substrings to database")
        p_bt2 = subprocess.run(
            [
                "bowtie2",
                "-p",
                str(threads),
                "-x",
                str(db_bt2),
                "-f",
                "-",
                "--seed",
                "42",
                "--very-sensitive",
                "-k",
                "32",
                "--np",
                "1",
                "--mp",
                "1,1",
                "--rdg",
                "0,1",
                "--rfg",
                "0,1",
                "--score-min",
                "L,0,-0.05",
                "--no-head",
                "--no-unal",
                "-S",
                str(kmers_sam),
            ],
            stdin=p_kmer.stdout,
        )
        print(" ".join(p_bt2.args))
        if p_bt2.returncode:
            raise RuntimeError("bowtie2 failed")
    return kmers_sam


def _process_regions(db_fna_staged, ex_out, output):
    awk_tab = ex_out / "awk.aggregated"
    for kmers_sam in ex_out.glob("*.sam"):
        with open(kmers_sam) as f_in, open(awk_tab, "w") as f_out:
            for line in f_in:
                line = line.strip().split("\t")
                f_out.write(
                    f"{line[2]}\t{line[3]}\t{int(line[3]) + len(line[9])}\t{line[0]}\t{line[9]}\n"
                )

    p_sort = subprocess.Popen(
        [
            "sort",
            "--parallel",
            "8",
            "--buffer-size=100g",
            "-k1,1",
            "-k2,2n",
            awk_tab,
        ],
        stdout=subprocess.PIPE,
    )
    sort_tab = awk_tab.with_suffix(".aggregated.sorted")
    with open(sort_tab, "w") as fp:
        p_dedup = subprocess.run(
            ["python", str(ROOT_DIR / "deduplicate.py")],
            stdin=p_sort.stdout,
            stdout=fp,
        )
    if p_dedup.returncode:
        raise RuntimeError("deduplicate.py failed")
    p_cat = subprocess.Popen(["cat", sort_tab], stdout=subprocess.PIPE)
    merge_bed = ex_out / "exhaustive.bed"
    with open(merge_bed, "w") as fp:
        p_merge = subprocess.run(
            ["bedtools", "merge", "-i", "stdin", "-c", "5", "-o", "count"],
            stdin=p_cat.stdout,
            stdout=fp,
        )
    if p_merge.returncode:
        raise RuntimeError("bedtools merge failed")

    trim_bed = ex_out / "exhaustive.trimmed.bed"
    p_trim = subprocess.run(
        [
            "python",
            str(ROOT_DIR / "trim-to-max-length.py"),
            db_fna_staged,
            merge_bed,
            trim_bed,
        ]
    )
    if p_trim.returncode:
        raise RuntimeError("trim-to-max-length.py failed")
    p_get = subprocess.run(
        [
            "bedtools",
            "getfasta",
            "-fi",
            db_fna_staged,
            "-bed",
            trim_bed,
            "-fo",
            output,
        ]
    )
    if p_get.returncode:
        raise RuntimeError("bedtools getfasta failed")

    return output


def _filter_db(database_fasta, contaminated_fasta, output_fasta, output_bt2, threads):

    p_filter = subprocess.run(
        [
            "python",
            str(ROOT_DIR / "filter_fna.py"),
            "--database-fasta",
            str(database_fasta),
            "--output-fasta",
            output_fasta,
            "--contaminated-fasta",
            contaminated_fasta,
        ]
    )
    if p_filter.returncode:
        raise RuntimeError("filter_fna.py failed")

    idx_path = output_bt2
    idx_path.parent.mkdir(exist_ok=True)
    p_build = subprocess.run(
        [
            "bowtie2-build",
            "--seed",
            "42",
            "--threads",
            str(threads),
            str(output_fasta),
            str(idx_path),
        ]
    )
    if p_build.returncode:
        raise RuntimeError("bowtie2-build failed")

    return output_fasta, idx_path


class OrderedGroup(click.Group):
    def __init__(self, name=None, commands=None, **attrs):
        super(OrderedGroup, self).__init__(name, commands, **attrs)
        #: the registered subcommands by their exported names.
        self.commands = commands or OrderedDict()

    def list_commands(self, ctx):
        return self.commands


@click.group(cls=OrderedGroup)
def cli():
    pass


@click.command(help="Export staged database fasta")
@shared_opts.db_fna
@click.option(
    "--db-fna-staged",
    "-s",
    help="the staged database in fasta format",
    required=True,
    type=click.Path(exists=False, path_type=Path),
)
@click.option(
    "--build_index",
    "-b",
    help="build bowtie2 index for staged database",
    is_flag=True,
)
def export_db(db_fna, db_fna_staged, build_index):
    _export_db(db_fna, db_fna_staged, build_index)


@click.command(help="Mapping kmer substrings to database")
@shared_opts.kmer
@shared_opts.jump
@shared_opts.files
@shared_opts.threads
@shared_opts.db_bt2
@click.option(
    "--ex_out",
    "-o",
    help="output directory for the sam files",
    type=click.Path(path_type=Path),
)
@click.option(
    "--task_id", "-i", help="task id for slurm array only", type=int, default=0
)
@click.option(
    "--task_count", "-n", help="task count for slurm array only", type=int, default=1
)
def process_kmers(kmer, jump, files, task_id, task_count, threads, db_bt2, kmers_sam):
    _process_kmers(kmer, jump, files, task_id, task_count, threads, db_bt2, kmers_sam)


@click.command(help="Process regions")
@click.option(
    "--db-fna-staged",
    "-d",
    help="the staged database in fasta format",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--ex-out",
    "-e",
    help="directory containing the sam files and saving the output files",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--output",
    "-o",
    help="output fasta for the contaminated regions",
    type=click.Path(path_type=Path),
)
def process_regions(db_fna_staged, ex_out, output):
    _process_regions(db_fna_staged, ex_out, output)


@click.command(
    help="Mask database fasta with contaminated regions and build bowtie2 index"
)
@click.option(
    "--database-fasta",
    "-d",
    help="the database in fasta format to filter",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--contaminated-fasta",
    "-c",
    help="contaminated fasta file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-fasta",
    "-o",
    help="output fasta file",
    required=True,
    type=click.Path(path_type=Path),
)
@click.option(
    "--output-bt2",
    "-m",
    help="output bowtie2 index",
    required=True,
    type=click.Path(path_type=Path),
)
@shared_opts.threads
def filter_db(database_fasta, contaminated_fasta, output_fasta, output_bt2, threads):
    _filter_db(database_fasta, contaminated_fasta, output_fasta, output_bt2, threads)


@click.command(
    help="Run exhaustive cleaning (export-db -> process-kmers -> process-regions -> filter-db)",
    context_settings={"show_default": True},
)
@click.option(
    "--base_label",
    "-l",
    help="base label for output files",
    type=click.Path(path_type=Path),
    default="exhaustive",
)
@shared_opts.db_fna
@click.option(
    "--db-bt2",
    "-b",
    help="the bowtie2 index of the database to map against, it will be built if not provided",
)
@shared_opts.files
@shared_opts.kmer
@shared_opts.jump
@shared_opts.threads
@click.option("--clean", "-c", help="remove temporary files", is_flag=True)
def submit(base_label, db_fna, db_bt2, files, kmer, jump, threads, clean):
    print("Starting exhaustive cleaning")
    base_label, out_dir = parse_base_label(base_label)
    ex_out = out_dir / f"{base_label}_{kmer}_{jump}"
    ex_out.mkdir(exist_ok=True, parents=True)

    db_fna_staged = ex_out / "concat.fna"
    db_fna_staged, db_bt2_new = _export_db(db_fna, db_fna_staged, db_bt2 is None)
    print(f"Staged database fasta: {db_fna_staged}")
    if db_bt2:
        print(f"Provided database index: {db_bt2}")
    else:
        print(f"Staged database index: {db_bt2_new}")
        db_bt2 = db_bt2_new

    kmers_sam = ex_out / "kmers.sam"
    kmers_sam = _process_kmers(kmer, jump, files, 0, 1, threads, db_bt2, kmers_sam)
    print(f"Kmer substrings mapped to database: {kmers_sam}")
    if os.stat(kmers_sam).st_size == 0:
        sys.exit(f"No kmer substrings mapped to database")

    contam_fna = _process_regions(
        db_fna_staged, ex_out, output=out_dir / f"{base_label}.contaminated.fna"
    )
    print(f"Contaminated regions: {contam_fna}")

    db_fna_masked, db_bt2_masked = _filter_db(
        db_fna_staged,
        contam_fna,
        out_dir / f"{base_label}.masked.fna",
        out_dir / f"{base_label}-bt2" / f"{base_label}-bt2",
        threads,
    )
    print(f"Masked db fasta: {db_fna_masked}")
    print(f"Masked db index: {db_bt2_masked}")

    if clean:
        print("Cleaning up temporary files")
        os.system(f"rm -r {ex_out}")
    print("Exhaustive done.")


cli.add_command(export_db)
cli.add_command(process_kmers)
cli.add_command(process_regions)
cli.add_command(filter_db)
cli.add_command(submit)

if __name__ == "__main__":
    cli()
