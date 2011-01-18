#!/usr/bin/env python
from PIL import Image
import itertools
from common import *


def sort_chrs_by_size():
    return sorted(CHRS.iteritems(), key=lambda (k, v): (v, k), reverse=True)

def foldr(f, i):
    return lambda s: reduce(f, s, i)

def main():
    nrows = 6
    ncols = 4
    chr_imgs = []
    for chr, _ in CHRS.iteritems():
        img = Image.open('out/' + chr + '.png', 'r').convert("RGBA")
        chr_imgs.append(img)
    max_w, max_h = chr_imgs[0].size
    chr_trans_bg = []
    img_w, img_h = chr_imgs[0].size
    bg = Image.new('RGBA', (img_w*ncols, img_h*nrows), (255,255,255,255))
    bg_w, bg_h = bg.size
    row_idx = 0
    for col_idx in range(len(chr_imgs)):
        x_pos = (col_idx % ncols) * img_w
        y_pos = (row_idx % nrows) * img_h
        bg.paste(chr_imgs[col_idx], (x_pos, y_pos))
        if (col_idx % ncols) == 0:
            row_idx += 1
    bg.save("out/all_chr.png")
    
if __name__ == '__main__':
    main()
