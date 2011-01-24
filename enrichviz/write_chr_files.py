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

layers_overflow = hide
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
#scale_log_base = 10
min = 0.0
max = 350.0
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
r0 = 0.92r
r1 = .98r
label_size = 10p
label_font = condensed
padding = 0p
rpadding = 0p
</plot>

<plot>
show = yes
type = tile
file = ../data/knowngene.bed
layers = 12
margin = 0.02u
thickness = 10
orientation = out
color = black
r0 = 0.87r
r1 = 0.90r
background = no
</plot>
"""

EXPR_TRACK="""
<plot>
show = yes
type = scatter
file = ../data/{expr}
glyph = circle
glyph_size = 5
fill_color = {color}
stroke_color = {color}
stroke_thickness = 1
background = no
min = 0.0000
max = 10.00
r0 = 0.81r
r1 = 0.85r
axis = yes
axis_color = lgrey
axis_thickness = 1
axis_spacing = 5
</plot>
"""

MK4_TRACK="""
<plot>
show = yes
type = tile
file = ../data/{mk4}
layers = 15
thickness = 15
margin = 0.02u
color = {color}
background = no
r0 = {r0}r
r1 = {r1}r
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
            r0 = r1 + .001
            r1 = r0 + .05
        f.write(GENE_TRACK)
        colors = ['nyt_orange', 'nyt_yellow']
        for i in range(len(RNA_DATASETS)):
            dset = RNA_DATASETS[i]
            dset_name = dset + "_exp.txt"
            color = colors[i]
            f.write(EXPR_TRACK.format(expr=dset_name, color=color))
        f.write(MK4_TRACK.format(mk4=MK4_DATASETS[0], r0=0.45, r1=0.475, color='nyt_orange'))
        f.write(MK4_TRACK.format(mk4=MK4_DATASETS[1], r0=0.475, r1=0.5, color='nyt_yellow'))
        f.write("</plots>\n")
        f.close() 
        
        
if __name__ == '__main__': 
    main()


