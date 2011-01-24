#!/usr/bin/env python
import math
from common import *

def read_knowngene(kgpath):
    kg = {}
    for line in open(kgpath, 'r'):
        fields = line.rstrip().split('\t')
        chr = fields[0]
        start = int(fields[1])
        end = start + 1 
        end = int(fields[2])
        gid = fields[3]
        kg[gid] = (chr, start, end)
    return kg

def read_gene_exp(exppath):
    exp = {}
    for line in open(exppath, 'r'):
        fields = line.rstrip().split('\t')
        gid = fields[0]
        expval = float(fields[1])
        exp[gid] = math.log(expval+1)
    return exp

def main():
    all_exp = []
    kg = read_knowngene("data/"+KNOWNGENE)
    for dset in RNA_DATASETS:
        basename = dset.split('.')[0]
        out = open("data/"+basename+"_exp.txt", 'w')
        dset_exp = read_gene_exp("data/"+dset+'.txt')
        for gid, val in dset_exp.iteritems():
            if gid in kg and gid in dset_exp:
                chr, start, end = kg[gid]
                val = dset_exp[gid]
                out.write("{chr}\t{start}\t{end}\t{val}\n".format(chr=chr, start=start, end=end, val=val))
        out.close()

if __name__=='__main__':
    main()
