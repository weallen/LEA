import numpy as np
from matplotlib import pyplot as plt

DATA_PATH="/gpfs/home/wallen/experiment/experiment/stavros_data/norm"

DATASETS=["icam_hmedip_norm.wig", "moe_hmedip_norm.wig",  "omp_hmedip_norm.wig", "moe_ac3_hmedip_norm.wig", "ngn_hmedip_norm.wig"]

def parse_wig(fname):
    f = open(fname, 'r')
    f.next()
    curr_chr = ''
    #fixedStep chrom=chr1 start=1 step=100 span=100
    for line in f:
        if re.match('^fixedStep', line):            
            curr_chr = line.split()[1].split('=')[1]
    


