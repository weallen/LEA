#!/usr/bin/env python
import sys
KGXREF = "/gpfs/home/wallen/data/genome_data/kgXref.txt"
KNOWNGENE = "/gpfs/home/wallen/data/genome_data/knownGene.txt"

class KnownGenes:
    def __init__(self):
        self.tx = {}
        self.cds = {}
        self.exons = {}
        self.kgxref = {}
        
    def add_exon(self, tx_name, chr, estart, eend):
        if tx_name not in self.exons:
            self.exons[tx_name] = []
        self.exons[tx_name].append((estart, eend))

    def add_tx(self, gene_name, tx_name, chr, strand, txstart, txend):
        if gene_name not in self.kgxref:
            self.kgxref[gene_name] = []
        self.kgxref[gene_name].append(tx_name)
        self.tx[tx_name] = (chr, strand, txstart, txend)

    def add_cds(self, tx_name, chr, cdsstart, cdsend):
        self.cds[tx_name] = (chr, cdsstart, cdsend)

    def exons():
        exons = [(name, data) for name, data in self.exons.iteritems()]
        return self.exons

    def tx():
        tx = [(name, data) for name, data in self.cds.iteritems()]
        return tx

    def cds():
        cds = [(name, data) for name, data in self.cds.iteritems()]
        return cds
    
    def promoters(self):
        pass

    def introns(self):
        pass

    def five_utrs(self):
        pass

    def three_utrs(self):
        pass

    def genes(self):
        genes = {}
        for gid, txids in self.kgxref.iteritems():
            chr = ""
            start = 1e10
            end = -1
            for txid in txids:                
                txchr, strand, txstart, txend = self.tx[txid]
                chr = txchr
                if start > txstart:
                    start = txstart
                if end < txend:
                    end = txend
            genes[gid] = (chr, start, end)
        return genes

    def write_gene_file(self, bedname):
        f = open(bedname, 'w')
        for gid, loc in self.genes().iteritems():
            gid_rep = gid.replace(' ','_')
            f.write("{chr}\t{start}\t{end}\t{name}\n".format(chr=loc[0], start=loc[1], end=loc[2], name=gid_rep))
        f.close()
        
def main():
    kg = KnownGenes()    
    kgfile = open(KNOWNGENE, 'r')
    kgxref_file = open(KGXREF, 'r')
    kgxref = {}
    for line in kgxref_file:
        txid, gid = line.rstrip().split('\t')
        kgxref[txid] = gid
        
    for line in kgfile:
        name, chr, strand, txstart, txend, cdsstart, cdsend, exoncount, \
              exonstarts, exonends, _, _ = line.rstrip().split("\t")
        kg.add_tx(kgxref[name], name, chr, strand, int(txstart), int(txend))
        kg.add_cds(name, chr, int(cdsstart), int(cdsend))
        es = exonstarts.split(',')
        ee = exonends.split(',')
        for i in range(int(exoncount)):
            kg.add_exon(name, chr, es[i], ee[i])
    kg.write_gene_file('knowngenes.txt')
    kgxref_file.close()
    kgfile.close()
if __name__ == '__main__':
    main()
    
