#!/usr/bin/env python
from optparse import OptionParser

def parse_kg(kgfname):
    f = open(kgfname, 'r')
    locs = []
    for line in f:
        fields = line.rstrip().split('\t')
        chr = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        name = fields[3]
        locs.append((chr, start, end, name))
    return locs

def main():
    parser = OptionParser()
    parser.add_option("-e", "--extend", action="store", dest="extend")
    (options, args) = parser.parse_args()
    locs = parse_kg(args[0])
    for l in locs:
        s = "{chr}\t{start}\t{end}\t{name}"
        print s.format(chr=l[0], start=l[1] - int(options.extend), end=l[1] + int(options.extend), name=l[3])
        
if __name__ == '__main__':
    main()
