import sys

db = sys.argv[1]
in_ = sys.argv[2]
out_ = sys.argv[3]

lengths = {}
id_ = None
db = open(db)
id_ = iter(db)
seqs = iter(db)
for i, s in zip(id_, seqs):
    i = i.strip().split()[0][1:]
    lengths[i] = len(s) - 1  # don't count newline

with open(in_) as infp, open(out_, 'w') as outfp:
    for line in infp:
        gid, start, stop, _ = line.strip().split('\t')
        stop = int(stop)

        if stop > lengths[gid]:
            stop = lengths[gid]

        outfp.write("%s\t%s\t%d\t%s\n" % (gid, start, stop, noclue))
