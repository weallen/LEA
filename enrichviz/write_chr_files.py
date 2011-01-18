#!/usr/bin/python

from common import *

HEADER = """<<include ../base.conf>>
<image>\n
dir = ../out
file  = {chrname}.png
24bit = yes
png = yes
#svg = yes
# radius of inscribed circle in image
radius         = {radius}p
background     = white
angle_offset   = -90

auto_alpha_colors = yes
auto_alpha_steps  = 5
</image>

chromosomes = {chrname}

<plots>
type = histogram
"""

PLOT="""
<plot>
# FOR HISTOGRAM
thickness = 1p
color = black
fill_under = yes
fill_color = bs8-2

# FOR HEATMAP
#color =  bd3-2, bd3-3, bd3-4, bd3-5, bd3-6, bd3-7, bd3-8, bd3-9, bd3-10, bd3-11
#stroke_thickness = 0
#min = 0.0
#max = 1.0
r0 = {r0}r
r1 = {r1}r
file = ../data/{dset}
</plot>
"""


def main():
    max_chr_len = 197195432.0 
    print max_chr_len 
    for chr, size in CHRS.iteritems(): 
        f = open('chrs/' + chr + '.conf', 'w') 
        radius = 1500 * (size / max_chr_len) 
        print chr, size/max_chr_len 
        f.write(HEADER.format(chrname=chr, radius=1500))
        r0 = 0.5
        r1 = r0 + .05
        for track in TRACKS:            
            f.write(PLOT.format(dset=track + ".txt", r0=r0, r1=r1))
            r0 = r1 + .01
            r1 = r0 + .05
        f.write("</plots>\n")
        f.close()
        
        
if __name__ == '__main__': 
    main()


