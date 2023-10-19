"""Enumerate kmers

This script generates all successive kmers, with a controllable offset,
from a given set of sequences.

If we specify k=6, and j=3, the following kmers will be emitted

------
   ------
      ------
         ------
            -----
AAAAAAAAAAAAAAAAA

Note that 3' terminal substrings will be emitted if the are > j.

The parallel mode of operation for this script breaks up an input sequence
into "number_of_task" chunks of approximately equal size. Kmers within each
chunk are generated. Each chunk region is extended by j to overlap chunks.

Input sequence can be single or multiline FASTA, and read from stdin.

Sequences obtained are sent to stdout.

Unit tests can be executed by calling the program without arguments.

Author: Daniel McDonald (d3mcdonald@eng.ucsd.edu)
License: BSD-3
"""

import sys
import datetime
import logging
import os
import io

# modified from
# https://stackoverflow.com/a/2135920/19741
def split(a, n):
    """get n roughly even start/stop boundaries for size a"""
    assert a > n
    # k -> quotient
    # m -> remainder
    k, m = divmod(a, n)

    # (start, stop)
    # for start, move to the opening position of the ith
    # quotient bin. we right shift by one, until we exhaust our remainder
    # the effect being that the first m bins contain an additional element
    # for stop, repeat but assuming the next ith position
    return [((i * k) + min(i, m),
            ((i + 1) * k) + min(i + 1, m))
            for i in range(n)]


def check_write(id_, rec, k, j, task, ntask, buf, max_n=100):
    if id_ is None:
        return

    rec = ''.join(rec)
    idfmt = id_ + '_%d\n'

    # j -> jump
    # k -> length
    # >>> a
    # 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    # >>> k = 15
    # >>> j = 8
    # >>> for i in range(0, len(a) - j + 1, j):
    # ...   print(i, a[i:i+k])
    # ...
    # 0 abcdefghijklmno
    # 8 ijklmnopqrstuvw
    # 16 qrstuvwxyzABCDE
    # 24 yzABCDEFGHIJKLM
    # 32 GHIJKLMNOPQRSTU
    # 40 OPQRSTUVWXYZ
    for idx, (rec_start, rec_end) in enumerate(split(len(rec), ntask)):
        if idx % ntask == task:
            task_substring = rec[rec_start:rec_end + j]  # overlap to cover edges
            for start in range(0, len(task_substring) - j + 1, j):
                s = task_substring[start:start + k]
                if s.count('N') >= max_n:
                    continue
                if len(s) <= j:
                    continue

                buf.write(idfmt % (start + rec_start))
                buf.write(s)
                buf.write('\n')


def test_check_write():
    tests = [(('>foo', 'AATTGGCC', 2, 1, 0, 3),
              ['>foo_0', 'AA', '>foo_1', 'AT', '>foo_2', 'TT']),
             (('>bar', 'ATNNNTTGC', 2, 1, 0, 3),
              ['>bar_0', 'AT', '>bar_1', 'TN']),
             (('>bar', 'ATNNNTTGC', 2, 1, 0, 1),
              ['>bar_0', 'AT', '>bar_1', 'TN', '>bar_4', 'NT',
               '>bar_5', 'TT', '>bar_6', 'TG', '>bar_7', 'GC'])]
    for args, exp in tests:
        buf = io.StringIO()
        check_write(*args, buf, max_n=2)
        buf.seek(0)
        obs = [l.strip() for l in buf.readlines()]
        assert obs == exp


def test_split():
    tests = [((10, 5), [(0, 2), (2, 4), (4, 6), (6, 8), (8, 10)]),
             ((11, 5), [(0, 3), (3, 5), (5, 7), (7, 9), (9, 11)]),
             ((11, 2), [(0, 6), (6, 11)])]
    for (a, k), exp in tests:
        obs = split(a, k)
        assert obs == exp


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print("testing...")
        test_check_write()
        test_split()
        sys.exit(0)
    elif len(sys.argv) != 5:
        print("usage: python kmer.py <k> <j> <task> <number_of_tasks>")
        print("")
        print("\tk: kmer size to emit")
        print("\tj: size of jump to perform")
        print("\ttask: the task index (e.g., SLURM_ARRAY_TASK_ID)")
        print("\tnumber_of_tasks: the total number of tasks (e.g., SLURM_ARRAY_TASK_COUNT)")
        sys.exit(1)

    n = datetime.datetime.now()
    nfmt = n.strftime("%Y.%m.%d")
    sjob = os.environ.get('SLURM_ARRAY_JOB_ID', 'no-slurm-id')
    stask = os.environ.get('SLURM_ARRAY_TASK_ID', 'no-slurm-array')

    logging.basicConfig(filename=f'kmer-{nfmt}-{sjob}_{stask}.log',
                        encoding='utf-8',
                        level=logging.DEBUG,
                        format='%(asctime)s:::%(message)s',
                        datefmt="%Y.%m.%d %H:%M:%S")

    k = int(sys.argv[1])
    j = int(sys.argv[2])
    task = int(sys.argv[3])
    ntask = int(sys.argv[4])

    id_ = None
    rec = []

    buf = sys.stdout
    for line in sys.stdin:
        if id_ is None and line.startswith('>'):
            id_ = line.strip().split(' ')[0]
        elif line.startswith('>'):
            logging.info(f"Emitting on: {id_}")
            check_write(id_, rec, k, j, task, ntask, buf)
            id_ = line.strip().split(' ')[0]
            rec = []
        else:
            rec.append(line.strip())

    logging.info(f"Emitting on: {id_}")
    check_write(id_, rec, k, j, task, ntask, buf)
