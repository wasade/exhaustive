import re
import skbio
import click
import sys
import io
import numpy as np
import numpy.testing as npt
from filter_sam import gather_coordinates


def mask(seq, coordinates):
    seq = np.array(list(seq))
    for start, stop in zip(*coordinates):
        seq[start:stop] = 'N'
    return ''.join(seq)


def test_mask():
    tests = [("AAAAAAAAA", (np.array([0, 3]), np.array([1, 5])), 'NAANNAAAA'),
             ("AAAAAAAAA", (np.array([0, 3]), np.array([3, 5])), 'NNNNNAAAA'),
             ("AAAAAAAAA", (np.array([0, 3]), np.array([1, 9])), 'NAANNNNNN')]

    for seq, coords, exp in tests:
        obs = mask(seq, coords)
        assert obs == exp


test_mask()


@click.command()
@click.option('--database-fasta', type=click.Path(exists=True), required=True)
@click.option('--output-fasta', type=click.Path(exists=False), required=True)
@click.option('--contaminated-fasta', type=click.Path(exists=True), required=True)
def main(database_fasta, output_fasta, contaminated_fasta):
    coordinates = gather_coordinates(open(contaminated_fasta))

    with open(output_fasta, 'w') as out_fp:
        with open(database_fasta) as in_fp:
            ids = iter(in_fp)
            seqs = iter(in_fp)
            for id_, seq in zip(ids, seqs):
                id_ = id_.strip()
                seq = seq.strip()

                assert id_.startswith('>')
                id_ = id_[1:]

                if id_ in coordinates:
                    seq = mask(seq, coordinates[id_])

                out_fp.write(f">{id_}\n{seq}\n")


if __name__ == '__main__':
    main()
