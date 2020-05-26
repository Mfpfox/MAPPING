#!/usr/bin/env python
# -*- coding: utf-8 -*-
# COORDINATES : Pmap_dbNSFP_mapped_CADD19.py
# Created on 10/6/19 MFP
import sys
from OpenTextFile import OpenTextFile

# python Pmap_dbNSFP_mapped_CADD19.py $i cktargetmatched_chr${i}.csv
order = sys.argv[1]
chrID = 'chr{}'.format(order)

# path to cadd files
PATH = '/u/home/m/mfpalafo/project-arboleda/CADD/whole_genome_SNVs/'

# opentext only read
filename = sys.argv[2]
dbfile = iter(OpenTextFile(filename, separator=','))
header = next(dbfile)
header = ','.join(header) + ',' + 'CADD_phred_hg19' + '\n'  # header now includes phred19

dbdict = {}     # turn file into key value structure
for line in dbfile:
    dbdict[line[0]] = ','.join(line[1:])
    # key = position, value = all other columns, joinign list with sep=','

outfile = '{}_CADD19_to_dbNSFP.csv'.format(chrID)
with open(outfile, 'wt') as writefile:
    writefile.write(header)
    print(f"Working on {chrID}")
    caddfilename = '{}_CADD.tsv.gz'.format(chrID)
    caddfile = iter(OpenTextFile(PATH+caddfilename))
    next(caddfile)
    for line in caddfile:
        if line[0] in dbdict:
            result = line[0] + ',' + dbdict[line[0]] + ',' + line[1] + '\n'
            writefile.write(result)

"""
%%time
# testing .py before job array
# python Pmap_dbNSFP_mapped_CADD19.py $i cktargetmatched_chr${i}.csv
order = 22
chrID = 'chr{}'.format(order)
# path to cadd files
PATH = '/u/home/m/mfpalafo/project-arboleda/CADD/whole_genome_SNVs/'
# opentext only read
filename = "cktargetmatch_chr22.csv"
dbfile = iter(OpenTextFile(filename, separator=','))
header = next(dbfile)
header = ','.join(header) + ',' + 'CADD_phred_hg19' + '\n'  # header now includes phred19
dbdict = {}     # turn file into key value structure
for line in dbfile:
    dbdict[line[0]] = ','.join(line[1:])
    # key = position, value = all other columns, joinign list with sep=','
outfile = '{}_CADD19_to_dbNSFP.csv'.format(chrID)
with open(outfile, 'wt') as writefile:
    writefile.write(header)
    print(f"Working on {chrID}")
    caddfilename = '{}_CADD.tsv.gz'.format(chrID)
    caddfile = iter(OpenTextFile(PATH+caddfilename))
    next(caddfile)
    for line in caddfile:
        if line[0] in dbdict:
            result = line[0] + ',' + dbdict[line[0]] + ',' + line[1] + '\n'
            writefile.write(result)
"""