#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created on 10.09.20

# called by jobarray.sh to map CADD37 phred and raw chr parsed file annotations to dbNSFP C postions files
# run on hoffman

import sys
from OpenTextFile import OpenTextFile

order = sys.argv[1]
chrID = 'chr{}'.format(order)
# path to cadd files on server- updated chr chunks with pos_id phred and raw scores
PATH = '/u/home/m/mfpalafo/project-arboleda/CADD/whole_genome_SNVs_w_Annotations/'

# opentext only read
filename = sys.argv[2] # targetmatch file
dbfile = iter(OpenTextFile(filename, separator=','))
header = next(dbfile)
header = ','.join(header) + ',' + 'CADD37_raw' + ',' + 'CADD37_phred' + '\n'


dbdict = {}
for line in dbfile:
    dbdict[line[0]] = ','.join(line[1:])
    # key = position id, value = all other columns in target match file

outfile = '{}_Heta_targetmatch_CADD37.csv'.format(chrID)

with open(outfile, 'wt') as writefile:
    writefile.write(header)
    print(f"Working on {chrID}")
    caddfilename = '{}_CADD37unique_keyid_raw_phred.tsv.gz'.format(chrID)
    caddfile = iter(OpenTextFile(PATH+caddfilename))
    next(caddfile)
    for line in caddfile:
        if line[0] in dbdict:
            result = line[0] + ',' + dbdict[line[0]] + ',' + line[1] + ',' + line[2] + '\n'
            writefile.write(result)





