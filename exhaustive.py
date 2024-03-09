#!/usr/bin/env python3

import os
import sys
import click
import subprocess
from pathlib import Path

ROOT_DIR = Path(__file__).parent.resolve()


@click.command(
    help="Run exhaustive mapping of kmer substrings to a database.",
    context_settings={"show_default": True},
)
@click.option(
    "--base_label",
    "-l",
    help="base label for output files",
    type=click.Path(path_type=Path),
    default="exhaustive",
)
@click.option(
    "--db_fna",
    "-d",
    help="the database in fasta format to map against",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
)
@click.option(
    "--db_bt2",
    "-b",
    help="the bowite2 index of the database to map against, it will be built if not provided",
)
@click.option(
    "--files",
    "-f",
    help="the paths of genomes to generate substrings from",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option("--kmer", "-k", help="kmer length", type=int, default=150)
@click.option("--jump", "-j", help="jump length", type=int, default=75)
@click.option("--threads", "-t", help="number of threads to use", type=int, default=8)
@click.option("--clean", "-c", help="remove temporary files", is_flag=True)
def main(base_label, db_fna, db_bt2, files, kmer, jump, threads, clean):
    base_label = Path(base_label)
    out_dir = base_label.parent.resolve()
    out_dir.mkdir(exist_ok=True, parents=True)
    base_label = base_label.name

    tmp_dir = out_dir / f"{base_label}_{kmer}_{jump}"
    tmp_dir.mkdir(exist_ok=True)

    db_fna_stage = tmp_dir / "concat.fna"
    if db_fna.suffix == ".gz":
        db_fna_stage = db_fna.with_suffix("")
        os.system(f"gunzip -c {db_fna}.gz > {db_fna_stage}")
    elif db_fna.suffix == ".xz":
        db_fna_stage = db_fna.with_suffix("")
        os.system(f"xzcat {db_fna} > {db_fna_stage}")
    else:
        os.system(f"cp {db_fna} {db_fna_stage}")
    one_line_fna = db_fna_stage.with_suffix(".one-line.fna")
    with open(db_fna_stage) as f_in, open(
        db_fna_stage.with_suffix(".one-line.fna"), "w"
    ) as f_out:
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
    db_fna_stage = one_line_fna.rename(db_fna_stage)
    if db_bt2 is None:
        print("Building bowtie2 index for database")
        db_bt2 = db_fna_stage.with_suffix("")
        p_bt2 = subprocess.run(["bowtie2-build", str(db_fna_stage), str(db_bt2)])
    else:
        print(f"Using existing bowtie2 index {db_bt2}")

    ex_sam = tmp_dir / "exhaustive.sam"
    if ex_sam.is_file():
        print(f"Skipping {ex_sam} as it already exists")
    else:
        genomes = []
        with open(files) as fp:
            for line in fp:
                file = line.strip()
                genomes.append(file)
        p_cat = subprocess.Popen(["cat"] + genomes, stdout=subprocess.PIPE)
        p_kmer = subprocess.Popen(
            ["python", str(ROOT_DIR / "kmer.py"), str(kmer), str(jump), "0", "1"],
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
                str(ex_sam),
            ],
            stdin=p_kmer.stdout,
        )
        print(" ".join(p_bt2.args))
        if p_bt2.returncode:
            raise RuntimeError("bowtie2 failed")

    if os.stat(ex_sam).st_size == 0:
        sys.exit(f"No kmer substrings mapped to database")

    awk_tab = tmp_dir / "awk.aggregated"
    with open(ex_sam) as f_in, open(awk_tab, "w") as f_out:
        for line in f_in:
            line = line.strip().split("\t")
            f_out.write(
                f"{line[2]}\t{line[3]}\t{int(line[3]) + len(line[9])}\t{line[0]}\t{line[9]}\n"
            )

    p_sort = subprocess.Popen(
        ["sort", "--parallel", "8", "--buffer-size=100g", "-k1,1", "-k2,2n", awk_tab],
        stdout=subprocess.PIPE,
    )
    sort_tab = awk_tab.with_suffix(".aggregated.sorted")
    with open(sort_tab, "w") as fp:
        p_dedup = subprocess.run(
            ["python", str(ROOT_DIR / "deduplicate.py")], stdin=p_sort.stdout, stdout=fp
        )
    if p_dedup.returncode:
        raise RuntimeError("deduplicate.py failed")
    p_cat = subprocess.Popen(["cat", sort_tab], stdout=subprocess.PIPE)
    merge_bed = tmp_dir / "exhaustive.bed"
    with open(merge_bed, "w") as fp:
        p_merge = subprocess.run(
            ["bedtools", "merge", "-i", "stdin", "-c", "5", "-o", "count"],
            stdin=p_cat.stdout,
            stdout=fp,
        )
    if p_merge.returncode:
        raise RuntimeError("bedtools merge failed")

    contaminated = out_dir / f"{base_label}.contaminated.fna"
    trim_bed = tmp_dir / "exhaustive.trimmed.bed"
    p_trim = subprocess.run(
        [
            "python",
            str(ROOT_DIR / "trim-to-max-length.py"),
            db_fna_stage,
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
            db_fna_stage,
            "-bed",
            trim_bed,
            "-fo",
            contaminated,
        ]
    )
    if p_get.returncode:
        raise RuntimeError("bedtools getfasta failed")

    db_fna_masked = out_dir / f"{base_label}.masked.fna"

    p_filter = subprocess.run(
        [
            "python",
            str(ROOT_DIR / "filter_fna.py"),
            "--database-fasta",
            str(db_fna_stage),
            "--output-fasta",
            db_fna_masked,
            "--contaminated-fasta",
            contaminated,
        ]
    )
    if p_filter.returncode:
        raise RuntimeError("filter_fna.py failed")

    idx_path = out_dir / f"{base_label}-bt2" / f"{base_label}-bt2"
    idx_path.parent.mkdir(exist_ok=True)
    p_build = subprocess.run(
        [
            "bowtie2-build",
            "--seed",
            "42",
            "--threads",
            str(threads),
            str(db_fna_masked),
            str(idx_path),
        ]
    )
    if p_build.returncode:
        raise RuntimeError("bowtie2-build failed")

    print(f"Masked db fasta: {db_fna_masked}")
    print(f"Masked db index: {idx_path}")
    if clean:
        print("Cleaning up temporary files")
        os.system(f"rm -r {tmp_dir}")
    print("Exhaustive done.")


if __name__ == "__main__":
    main()
