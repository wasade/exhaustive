import re
import click
import sys
import io
import numpy as np
import numpy.testing as npt


def mask(seq, coordinates):
    seq = np.array(list(seq))
    for start, stop in zip(*coordinates):
        seq[start:stop] = 'N'
    return ''.join(seq)


def gather_coordinates(recs):
    ex = re.compile(r"^>(\S+):(\d+)-(\d+)")
    coords = {}

    for r in recs:
        if r.startswith('>'):
            try:
                gid, start, stop = ex.match(r).groups()
            except:
                raise ValueError(f"Bad ID: {r}")

            if gid not in coords:
                coords[gid] = [[], []]

            coords[gid][0].append(int(start))
            coords[gid][1].append(int(stop))

    return {gid: (np.array(a), np.array(b))
            for gid, (a, b) in coords.items()}


def test_gather_coordinates():
    def eq(d1, d2):
        assert d1.keys() == d2.keys()
        for k in d1:
            npt.assert_equal(d1[k][0], d2[k][0])
            npt.assert_equal(d1[k][1], d2[k][1])

    tests = [
        (io.StringIO(">x:1-10\nAAAA\n>x:200-300\nTTTT\n>y:30-40\nGGG\n"),
         {'x': (np.array([1, 200]), np.array([10, 300])),
          'y': (np.array([30, ]), np.array([40, ]))}),
    ]

    exceptions = [
        (io.StringIO(">foobar\nAATT\n"), ValueError),
        (io.StringIO(">x:10-20\nAA\n>y:5\nGGG\n"), ValueError)
    ]

    for test, exp in tests:
        obs = gather_coordinates(test)
        eq(obs, exp)


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
                id_ = id_[1:].split(" ")[0]
                if id_ in coordinates:
                    seq = mask(seq, coordinates[id_])

                out_fp.write(f">{id_}\n{seq}\n")


if __name__ == '__main__':
    main()
