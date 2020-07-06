#!/usr/bin/env python
# -*- coding: utf-8 -*-
# COORDINATES : Pmap_pos_overlap_CADD.py
# Created on 12/25/19 MFP

import sys
import csv
import os
from OpenTextFile import OpenTextFile

# argv 1 is chr; 2 is position file that matches chr
order = sys.argv[1] # $i
filename = sys.argv[2] # $f
assembly = int(sys.argv[3]) # 37 or 38
postype = sys.argv[4] # dect or not

chrID = 'chr{}'.format(order)
# input positions file
dbfile = iter(OpenTextFile(filename, separator=','))
header = next(dbfile)

if assembly == 37:
    # path to cadd files
    PATH = '/CADD/whole_genome_SNVs_w_Annotations/'
    # header = pos_hg19,chr,pos_ID
    header = ','.join(header) + ',' + 'Ref' + ',' + 'Alt' + ',' + 'Type' + ',' + 'Length' + ',' + 'AnnoType' + ',' + 'Consequence' + ',' + 'ConsScore' + ',' + 'ConsDetail' + ',' + 'GC' + ',' + 'CpG' + ',' + 'motifECount' + ',' + 'motifEName' + ',' + 'motifEHIPos' + ',' + 'motifEScoreChng' + ',' + 'oAA' + ',' + 'nAA' + ',' + 'GeneID' + ',' + 'FeatureID' + ',' + 'GeneName' + ',' + 'CCDS' + ',' + 'Intron' + ',' + 'Exon' + ',' + 'cDNApos' + ',' + 'relcDNApos' + ',' + 'CDSpos' + ',' + 'relCDSpos' + ',' + 'protPos' + ',' + 'relProtPos' + ',' + 'Domain' + ',' + 'Dst2Splice' + ',' + 'Dst2SplType' + ',' + 'minDistTSS' + ',' + 'minDistTSE' + ',' + 'SIFTcat' + ',' + 'SIFTval' + ',' + 'PolyPhenCat' + ',' + 'PolyPhenVal' + ',' + 'priPhCons' + ',' + 'mamPhCons' + ',' + 'verPhCons' + ',' + 'priPhyloP' + ',' + 'mamPhyloP' + ',' + 'verPhyloP' + ',' + 'GerpRS' + ',' + 'GerpRSpval' + ',' + 'GerpN' + ',' + 'GerpS' + ',' + 'Grantham' + ',' + 'Dist2Mutation' + ',' + 'Freq100bp' + ',' + 'Rare100bp' + ',' + 'Sngl100bp' + ',' + 'Freq1000bp' + ',' + 'Rare1000bp' + ',' + 'Sngl1000bp' + ',' + 'Freq10000bp' + ',' + 'Rare10000bp' + ',' + 'Sngl10000bp' + ',' + 'RawScore' + ',' + 'PHRED' + '\n'
    if postype == 'dect':
        outfile = '{}_CADD_GRCh37_DETECTED_CK.csv'.format(chrID)
    if postype == 'not':
        outfile = '{}_CADD_GRCh37_NOT_DETECTED_CK.csv'.format(chrID)
    # make position file into dictionary
    dbdict = {} # columns-- 0:position 1:chr 2:pos_id19 or pos_id38
    for line in dbfile:
        posi = line[0]
        dbdict[posi] = ','.join(line[1:])
    with open(outfile, 'wt') as writefile:
        writefile.write(header)
        print(f"Working on {chrID}")
    caddfilename = '{}_GRCh37_inclannotationsCADD.tsv.gz'.format(chrID)
    caddfile = iter(OpenTextFile(PATH+caddfilename, separator='\t'))
    next(caddfile) # skip header row
    for line in caddfile:
        # if pos in dict as key
        if line[1] in dbdict:
            result = []
            result.append(line[1])
            pos_id = dbdict[line[1]].split(",")
            result.extend(pos_id)
            result.extend((line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],line[16],line[17],line[18],line[19],line[20],line[21],line[22],line[23],line[24],line[25],line[26],line[27],line[28],line[29],line[30],line[31],line[32],line[33],line[34],line[35],line[36],line[37],line[38],line[39],line[40],line[41],line[42],line[43],line[44],line[65],line[66],line[67],line[68],line[92],line[93],line[94],line[95],line[96],line[97],line[98],line[99],line[100],line[101],line[102],line[105],line[106]))
            with open(outfile, 'a') as out:
                        csvWriter = csv.writer(out)
                        csvWriter.writerow(result)

if assembly == 38:
    PATH = '/CADD/GRCh38/'
    # from h.tsv files on hoffman; header = pos_hg38,chr,pos_ID
    header = ','.join(header) + ',' + 'Ref' + ',' +'Alt' + ',' +'Type' + ',' +'Length' + ',' +'AnnoType' + ',' +'Consequence' + ',' +'ConsScore' + ',' +'ConsDetail' + ',' +'GC' + ',' +'CpG' + ',' +'motifECount' + ',' +'motifEName' + ',' +'motifEHIPos' + ',' +'motifEScoreChng' + ',' +'oAA' + ',' +'nAA' + ',' +'GeneID' + ',' +'FeatureID' + ',' +'GeneName' + ',' +'CCDS' + ',' +'Intron' + ',' +'Exon' + ',' +'cDNApos' + ',' +'relcDNApos' + ',' +'CDSpos' + ',' +'relCDSpos' + ',' +'protPos' + ',' +'relProtPos' + ',' +'Domain' + ',' +'Dst2Splice' + ',' +'Dst2SplType' + ',' +'minDistTSS' + ',' +'minDistTSE' + ',' +'SIFTcat' + ',' +'SIFTval' + ',' +'PolyPhenCat' + ',' +'PolyPhenVal' + ',' +'priPhCons' + ',' +'mamPhCons' + ',' +'verPhCons' + ',' +'priPhyloP' + ',' +'mamPhyloP' + ',' +'verPhyloP' + ',' +'GerpRS' + ',' +'GerpRSpval' + ',' +'GerpN' + ',' +'GerpS' + ',' +'Grantham' + ',' +'Dist2Mutation' + ',' +'Freq100bp' + ',' +'Rare100bp' + ',' +'Sngl100bp' + ',' +'Freq1000bp' + ',' +'Rare1000bp' + ',' +'Sngl1000bp' + ',' +'Freq10000bp' + ',' +'Rare10000bp' + ',' +'Sngl10000bp' + ',' +'RawScore' + ',' +'PHRED' + '\n'
    if postype == 'dect':
        outfile = '{}_CADD_GRCh38_DETECTED_CK.csv'.format(chrID)
    if postype == 'not':
        outfile = '{}_CADD_GRCh38_NOT_DETECTED_CK.csv'.format(chrID)
    dbdict = {}
    for line in dbfile:
        posi = line[0]
        dbdict[posi] = ','.join(line[1:])
    with open(outfile, 'wt') as writefile:
        writefile.write(header)
        print(f"Working on {chrID}")
    caddfilename = '{}h.tsv.gz'.format(chrID)
    caddfile = iter(OpenTextFile(PATH+caddfilename, separator='\t'))
    next(caddfile) # skip header row
    for line in caddfile:
        if line[1] in dbdict:
            result = []
            result.append(line[1])
            pos_id = dbdict[line[1]].split(",")
            result.extend(pos_id)
            result.extend((line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],line[16],line[17],line[18],line[19],line[20],line[21],line[22],line[23],line[24],line[25],line[26],line[27],line[28],line[29],line[30],line[31],line[32],line[33],line[34],line[35],line[36],line[37],line[38],line[39],line[40],line[41],line[42],line[43],line[44],line[75],line[76],line[77],line[78],line[107],line[108],line[109],line[110],line[111],line[112],line[113],line[114],line[115],line[116],line[117],line[123],line[124]))
            with open(outfile, 'a') as out:
                        csvWriter = csv.writer(out)
                        csvWriter.writerow(result)

