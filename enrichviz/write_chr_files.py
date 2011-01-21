#!/usr/bin/python
import re
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

"""

PLOT="""
<plot>
type = histogram
# FOR HISTOGRAM
thickness = 1p
color = {color}
fill_under = yes
fill_color = {color}

# FOR HEATMAP
#scale_log_base = 2
min = 0.0
max = 1.0
r0 = {r0}r
r1 = {r1}r
file = ../data/{dset}
</plot>
"""

GENE_TRACK="""
<plot>
type = text
color = black
file = ../data/knowngenes.txt
r0 = 0.9r
r1 = 1.0r
label_size = 10p
label_font = condensed
padding = 0p
rpadding = 0p
</plot>
"""


def main():
    max_chr_len = 197195432.0 
    print max_chr_len
    rcols = ""
    for i in range(17):
        rcols += "pg" + str(i) + ","
    rcols += "pg17"
    gcols = ""
    for i in range(17):
        gcols += "rb" + str(i) + ","
    gcols += "rb17"
    for chr, size in CHRS.iteritems(): 
        f = open('chrs/' + chr + '.conf', 'w') 
        radius = 1500 * (size / max_chr_len) 
        print chr, size/max_chr_len 
        f.write(HEADER.format(chrname=chr, radius=2000))
        r0 = 0.5
        r1 = r0 + .05        
        for track in TRACKS:
            if re.search("_hmedip", track): 
                col = 'nyt_blue'
            else: 
                col = 'nyt_red'
            f.write(PLOT.format(dset=track + ".txt", r0=r0, r1=r1, color=col))
            r0 = r1 + .01
            r1 = r0 + .05
        f.write(GENE_TRACK)
        f.write("</plots>\n")
        f.close() 
        
        
if __name__ == '__main__': 
    main()


