#!/usr/bin/env python
import os
import re

ROOT_PATH = "/gpfs/home/wallen/"
BEDTOOLS = ROOT_PATH + "src/bedtools/bin/"
DATA = ROOT_PATH + "experiment/experiment/stavros_data/"
GENOME_DATA = ROOT_PATH + "data/genome_data/"

def load_kgxref():
    kgxref_path = GENOME_DATA + "kgXref.txt"
    kgxref = {}
    f = open(kgxref_path, 'r')
    for line in f:
        gid, sym = line.rstrip().split('\t')
        kgxref[gid] = sym
    return kgxref

def intersect_bed_feature_counts(fname):
    genes = {}
    f = open(fname, 'r')
    for line in f:
        fields = line.rstrip().split('\t')
        gid = fields[6].split("_")[0]
        if gid not in genes:
            genes[gid] = 0
        genes[gid] += 1
    return genes

def make_gene_lists(type):
    kgxref = load_kgxref()
    startdir = os.getcwd()
    if type == "lores":
        path = DATA + "diff_intersect_lores/"
    elif type == "hires":
	path = DATA + "diff_intersect_roi/"
    outpath = DATA + "diff_uniq_genes"
    txpath = DATA + "diff_tx_counts"
    os.chdir(path)
    for dirpath, dirs, files in os.walk(path):
        for name in files:
            if name == ".":
                break
            if name == "..":
                break
            if os.stat(name).st_size == 0:
                break
            feats = intersect_bed_feature_counts(name)
            if type == "lores":
                txfile = open(txpath + "/" + name.split('.')[0] + "_tx_counts.txt", 'w')                            
            else:
                txfile = open(txpath + "/" + name.split('.')[0] + "_tx_counts.txt", 'w')            
            for f, c in feats.iteritems():
                txfile.write("{a}\t{b}\n".format(a=f, b=c))
            genes = set()
            for f, c in feats.iteritems():
                genes.add(kgxref[f])
            if type == "lores":
                outfile = open(outpath + "/" + name.split('.')[0] + "_uniq_genes.txt", 'w')
            else:
                outfile = open(outpath + "/" + name.split('.')[0] + "_uniq_genes.txt", 'w')
            for g in genes:
                outfile.write("{g}\n".format(g=g))
    os.chdir(startdir)
    
def do_diff(file1, file2, outname):
    cmd = BEDTOOLS + "intersectBed" + " -a " + file1 + " -b " + file2 + " -wb > " + outname
    print cmd
    os.system(cmd)

def do_comparison(type):
    path = "."
    if type == "lores":
        path = DATA + "diff_enrich_lores/"
    elif type == "hires":
	path = DATA + "diff_enrich_roi/"
    startdir = os.getcwd()
    os.chdir(path) 
    for dirpath, dirnames, filenames in os.walk(path):
        for name in filenames:
            if name == ".":
                break
            if name == "..":
                break
	
            rootname = name.split('.')
            if type == "lores":
                outfile = DATA + "diff_intersect_lores/" + rootname[0] + "_intersect.bed"
            elif type == "hires":
                outfile = DATA + "diff_intersect_roi/" + rootname[0] + "_intersect.bed"
            index_name = ""
            print name
            if re.search("prom", name):
                print 'promoter'
                index_name = "kg_promoters_10kb.bed"                
            elif re.search("gene", name) or re.search("genes", name):
                index_name = "knowngene.bed"
            elif re.search('intron', name):
                print 'intron'
                index_name = "kg_introns.bed"
            elif re.search('exon', name):
                print 'exon'
                index_name = "kg_exons.bed";
            if index_name != "":
                do_diff(name, GENOME_DATA + "/" + index_name, outfile); 
    os.chdir(startdir) 

def main():
    os.chdir(DATA)
    do_comparison("hires")
    make_gene_lists("hires")
    
    do_comparison("lores")
    make_gene_lists("lores")
if __name__ == '__main__':
    main()
