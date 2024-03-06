import sys
import io

def dedup(in_, out_):
    last = None
    for line in in_:
        gid, start, stop, remainder = line.split(b'\t', 3)
        if (gid, start, stop) != last:
            out_.write(line)
            last = (gid, start, stop)


def test_dedup():
    in_ = io.BytesIO(
b"""a\tb\tc\tstuff\tfoo
a\tb\tc\tcool\tbar
d\te\tf\ta
e\td\tf\tb
e\td\tf\tc
e\td\tf\td
e\td\tf\te
x\ty\tz\taaa\tyyy\tbbb
""")
    out_ = io.BytesIO()
    exp = b"""a\tb\tc\tstuff\tfoo
d\te\tf\ta
e\td\tf\tb
x\ty\tz\taaa\tyyy\tbbb
"""

    dedup(in_, out_)
    out_.seek(0)
    obs = out_.read()
    assert obs == exp


if __name__ == '__main__':
    test_dedup()
    dedup(sys.stdin.buffer, sys.stdout.buffer)
