#!/usr/bin/env python
import sys
import getopt

import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from PIL import Image

DATA_PATH="/gpfs/home/wallen/experiment/experiment/stavros_data/norm"
SCALE_FACTOR = 100

def chr_to_num(chr):
    if chr == "chrX":
        return 20
    elif chr == "chrY":
        return 21
    elif chr == "chrM":
        return 22
    else:
        return int(chr[3:])

def file_length(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

class NormDipData:
    def __init__(self, num_windows, num_dsets):
        self.chr = np.zeros(num_windows, dtype='short')
        self.norm = np.zeros([num_windows, num_dsets])
        self.num_windows = num_windows
        self.num_dsets = num_dsets

    def load_chr(self, dd_path):
        genome_chr_path = dd_path + "/genome_chr.txt"
        print "Reading", dd_path, "chr..."
        with open(genome_chr_path) as chr:
            for i, l in enumerate(chr):
                self.chr[i] = chr_to_num(l.rstrip('\n'))
                
    def load_norm(self, dset_idx, dd_path):
        genome_norm_path = dd_path + "/genome_norm.txt"        
        print "Reading", dd_path, "norm..."
        with open(genome_norm_path) as norm:
            for i, l in enumerate(norm):
                self.norm[i, dset_idx] = float(l.rstrip('\n'))

def get_chr(ndipdata, chridx):
    return ndipdata.norm[ndipdata.chr == chridx]
    
def load_all_dipdata(dipnames):
    data = {}
    num_windows = file_length(dipnames[0] + "/genome_chr.txt")
    dd = NormDipData(num_windows, len(dipnames))
    dd.load_chr(dipnames[0])
    for dset_idx in range(len(dipnames)):
        print "Loading", dipnames[dset_idx]
        dd.load_norm(dset_idx, dipnames[dset_idx])
    return dd

# hscale must be int > 1
# vscale must be int > 1
def rescale(mtx, hscale, vscale):
    new_len = math.ceil(mtx.shape[0] / hscale)
    rescaled_h = np.zeros([new_len, mtx.shape[1]])
#    rescaled_v = np.zeros([new_len, mtx.shape[1]*vscale])
    for i in xrange(0, rescaled_h.shape[0]-1):
        rescaled_h[i, :] = np.mean(mtx[(i*hscale):((i+1)*hscale), :])
#    for j in xrange(0, rescaled_v.shape[1]):
#        old_idx = math.floor(j / vscale)
#        rescaled_v[:, j] = rescaled_h[:, old_idx]
    return np.repeat(rescaled_h, vscale, axis=1)

def plot_chr(ndd, chridx):
    res = rescale(currchr, 100)
    plt.figure()
    plt.imshow(res)
    
def main():
    dips = sys.argv[1:]
    dd = load_all_dipdata(dips)
    for chridx in range(1, 21):
        save_path = "norm_data_chr_"+str(chridx)+".png"
        currchr = rescale(get_chr(dd, chridx), SCALE_FACTOR, SCALE_FACTOR)
        plt.imsave(arr=currchr, fname=save_path)


if __name__=='__main__':
    main()
