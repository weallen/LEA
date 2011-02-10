#!/usr/bin/python
import re

KNOWN_GENE = "/gpfs/home/wallen/data/genome_data/knowngene.bed"
KGXREF = "/gpfs/home/wallen/data/genome_data/kgXref.txt"

def read_kgxref():
    f = open(KGXREF, 'r')
    kgxref = {}
    for line in f:
        fields = line.rstrip().split('\t')
        kgxref[fields[0]] = re.sub(' ', '_', fields[1])
    f.close()
    return kgxref

def read_kg():
    f = open(KNOWN_GENE, 'r')
    known_gene = {}
    for line in f:
        fields = line.rstrip().split('\t')
        chr = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        name = fields[3]
        known_gene[name] = (chr, start, end)
    f.close()
    return known_gene

def main():
    f = open('data/knowngenes.txt','w')
    kg = read_kg()
    kgxref = read_kgxref()
    genes = {}
    for name, pos in kg.iteritems():
        sym = kgxref[name]
        if not sym in genes:
            genes[sym] = pos
        else:
            if genes[sym][0] == pos[0] and \
               genes[sym][1] < pos[1]:
                genes[sym] = pos
    for sym, pos in genes.iteritems():
        f.write("{chr}\t{start}\t{end}\t{name}\n".format(chr=pos[0],\
                                                         start=pos[1],\
                                                         end=pos[2],\
                                                         name=sym))
    f.close()
    
if __name__ == '__main__':
    main()
