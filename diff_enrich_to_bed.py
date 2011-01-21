#!/usr/bin/env python
import sys
from optparse import OptionParser

FDR = 0.05
WINDOW_SIZE = 1000

def num_to_chr(n):
    if n == 21:
        return 'chrY'
    elif n == 20:
        return 'chrX'
    return 'chr' + str(n)

def make_bed(fdr, win_size, args):
    for fname in args:
        diffenrich = open(fname, 'r')
        basename = fname.split(".")[0]
        out = open(basename + '.bed', 'w')
        diffenrich.next()
        for line in diffenrich:
            chr, pos, qval = line.rstrip().split(',')
            if float(qval) < fdr:        
                out.write('{chr}\t{start}\t{end}\n'.format(chr=num_to_chr(int(chr)),\
                                                           start=int(pos), \
                                                           end=int(pos)+win_size))                
def main():
    parser = OptionParser() 
    parser.add_option('-f', '--fdr', action='store', type='float',\
                      dest='fdr', help="FDR rate")
    parser.add_option('-s', '--window-size', action='store', type='int',\
                      dest='win_size', help="window size")
    (options, args) = parser.parse_args()
    if not options.fdr:
        fdr = FDR
    if not options.win_size:
        win_size = WINDOW_SIZE
    make_bed(fdr, win_size, args)

if __name__ == '__main__':    
    main()
