# G6PD pymol (2bh9) colored by CADD max codon score

from pymol import stored
import re, sys, pymol 

infil = open("CADD38_maxcodon_G6PD_2BH9_pos27_515.txt", 'r')
newb = []
for line in infil.readlines(): newb.append(float(line))
infil.close()
alter 2bh9, b = 0.0
alter 2bh9  and n. CA, b=newb.pop(0) # can use alpha2all.py after this
cmd.spectrum("b", "white gray purple", "2bh9", minimum=6.2, maximum=35)
