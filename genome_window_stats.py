#!/usr/bin/env python
# Returns a per window count of the number of cpgs across the entire genome
# Adapted from "lookupDnaSequence.py" in the Epigenome Pipeline Package
# by Christoph Bock
import os
import optparse

CHRS = ["chr1", "chr2", "chr3", "chr4", "chr5",\
        "chr6", "chr7", "chr8", "chr9", "chr10",\
        "chr11", "chr12", "chr13", "chr14", "chr15",\
        "chr16", "chr17", "chr18", "chr19", "chrX",\
        "chrY"]
GENOME_PATH = "/gpfs/home/wallen/data/genome_data/mm9_genome"

def num_to_chr(n):
    if n == 21:
        return "chrY"
    if n == 20:
        return "chrX"
    return "chr" + str(n)

def seq_stats(seq):
    length = len(seq)
    non_rep_char = ['A','T','C','G']
    repeat = length
    for letter in non_rep_char:
        repeat -= seq.count(letter)
    seq = seq.upper()
    A = seq.count('A')
    C = seq.count('C')
    G = seq.count('G')
    T = seq.count('T')
    CpG = seq.count('CG')
    CpA = seq.count('CA')
    if G*C > 0:
        O_E = float(CpG*length)/(G*C)
    else:
        O_E = "NA"
    if CpG > 0:
        ratio_ca_vs_cg = float(CpA)/CpG
    else:
        ratio_ca_vs_cg = "NA"
    GC = float(G+C)/length
    return (float(A)/length, float(C)/length, float(G)/length,float(T)/length,\
            float(repeat)/length, GC, float(CpG)/length, float(CpA)/length, O_E,\
            ratio_ca_vs_cg)
            
def lookup_sequence(chrom_file, offset, chrom_start, chrom_end):
    pos = chrom_start / 50 * 51 + chrom_start % 50
    length = (chrom_end - chrom_start + 50) / 50 * 51
    chrom_file.seek(offset + pos, 0)
    seq = chrom_file.read(length)
    seq = seq.replace("\n", "")
    seq = seq[0:(chrom_end - chrom_start)]
    return seq
    
def open_chrom_file(genome_path, chrom):
    chrom_fname = genome_path + os.sep + chrom + '.fa'
    chrom_file = open(chrom_fname, 'r')
    first_byte = chrom_file.read(1)
    offset = 0
    if first_byte == ">":
        offset = len(chrom_file.readline()) + 1
        chrom_file.seek(offset+49, 0)
        seq = chrom_file.read(3)
    valid_letters = ['A','C','G','T','N']
    if not seq[0].upper() in valid_letters \
           or not seq[1] == "\n"\
           or not seq[2].upper() in valid_letters:
        print "Incorrect file format for " + chrom
        raise SystemExit
    return (chrom_file, offset)

def main():
    parser = optparse.OptionParser()
    parser.add_option('--winfile', action='store', type='string', dest='win_file')
    parser.add_option('--winsize', action='store', type='string', dest='win_size')
    parser.add_option('--outfile', action='store', type='string', dest='out_file')
    (options, args) = parser.parse_args()
    chr_files = {}
    offsets = {}
    for chr in CHRS:
        f, o = open_chrom_file(GENOME_PATH, chr)
        chr_files[chr] = f
        offsets[chr] = o
    outfile = open(options.out_file, 'w')
    outfile.write("chr\tstart\tend\tA\tC\tG\tT\trep\tGC\tCpG\tCpA\tO_E\tratioCAvsCG\n")
    outline = "{chr}\t{start}\t{end}\t{Ac}\t{Cc}\t{Gc}\t{Tc}\t{rep}\t{GC}\t{CpG}\t{CpA}\t{O_E}\t{ratio}\n"
    wf = open(options.win_file, 'r')

    count = 0
    for line in wf:
        if count % 1000000 == 0:
            print "Reading window", count
        chr_idx, pos = line.rstrip().split(',') 
        chr = num_to_chr(int(chr_idx))
        start = int(pos) 
        end = int(pos) + int(options.win_size) 
        chr_f = chr_files[chr] 
        offset = offsets[chr] 
#                return (float(A)/length, float(C)/length, float(G)/length,float(T)/length,\ 
#            float(repeat)/length, GC, CpG/length, CpA/length, O_E,\ 
#            ratio_ca_vs_cg) 
        stats  = seq_stats(lookup_sequence(chr_f, offset, start, end))
        outfile.write(outline.format(chr=int(chr_idx), start=start, end=end, Ac=stats[0], Cc=stats[1], Gc=stats[2], Tc=stats[3],\
                                     rep=stats[4], GC=stats[5], CpG=stats[6], CpA=stats[7], O_E = stats[8], ratio = stats[9]))
        count += 1
    wf.close()
    for chr, f in chr_files.iteritems():
        f.close()
    outfile.close()
    
if __name__ == '__main__':
    main()
