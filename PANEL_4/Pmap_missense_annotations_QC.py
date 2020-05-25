# !/usr/bin/env python3
# -*- coding: utf-8 -*-

# Pmap M local ipynb code, Pmap_missense_annotations_QC.py
# markdown M local QC of positions from dbNSFP overlapped with CADD37 or 38 annotations

import os
import sys
import pandas as pd

""" 
dbNSFP SCORE file column names:

columns detected file:
    Amino_acids
    pos_ID 'Q3ZCM7_C354'
    pos_id19 or pos_id38_x

columns not detected file:
    pos_id19 or pos_id38
    Amino_acids,
    pos_ID_falseCKtarget
"""

def create_coordinate_id(df, chrr, pos, ref, alt, assembly):
    # did not format position correctly with leading zeros!
    if assembly == 37:
        df.loc[:,'pos_id19'] = df[chrr].astype(str) + '_' + \
                df[pos].astype(str) + '_' + df[ref].astype(str) + \
                '_' + df[alt].astype(str)
    if assembly == 38:
        df.loc[:,'pos_id38'] = df[chrr].astype(str) + '_' + \
                df[pos].astype(str) + '_' + df[ref].astype(str) + \
                '_' + df[alt].astype(str)
    return df


def format_missense_triple(df, oaacol, naacol):
    #  A|A turns to Ala/Ala
    amino_dict = dict([('A', 'Ala'),('G', 'Gly'), ('I','Ile'), ('L','Leu'), ('P', 'Pro'), ('V','Val'), ('F','Phe'),('W', 'Trp'), ('Y', 'Tyr'), ('D','Asp'),('E','Glu'), ('R','Arg'),('H','His'), ('K','Lys'), ('S','Ser'), ('T', 'Thr'), ('C', 'Cys'), ('M', 'Met'), ('N', 'Asn'), ('Q','Gln')])
    df[oaacol].replace(amino_dict, inplace=True)
    df[naacol].replace(amino_dict, inplace=True)
    ccopy = df[naacol].copy()
    df['Amino_acids'] = df[oaacol].str.cat(ccopy, sep='/')
    return df


def filter_cadd_overlap(df, assembly):
    # [1] filter for missense only
    miss = df[df['Consequence'] == 'NON_SYNONYMOUS'].copy()
    # [2] new pos_id(assembly) to files
    if assembly == 37:
        miss = create_coordinate_id(miss, 'chr', 'pos_hg19', 'Ref', 'Alt', assembly)
    if assembly == 38:
        miss = create_coordinate_id(miss, 'chr', 'pos_hg38', 'Ref', 'Alt', assembly)
    # [3] new missense type column in 3 letter format with '/' sep {oAA, nAA}
    miss = format_missense_triple(miss, 'oAA', 'nAA')

    return miss

    # [4] concat all files from 37: DECT or NOT ... 38: DECT or NOT


def main():
    os.chdir('/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/CADDmapped/RESULT_pos_overlap_dbNSFPcoordinates')

    chrls = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']

    # GRCh37
    for order in chrls:
        chrID = 'chr{}'.format(order)
        file1 = '{}_CADD_GRCh37_DETECTED_CK.csv'.format(chrID)
        out1 = 'MISSENSE_{}_CADD_GRCh37_DETECTED_CK.csv'.format(chrID)
        file2 = '{}_CADD_GRCh37_NOT_DETECTED_CK.csv'.format(chrID)
        out2 = 'MISSENSE_{}_CADD_GRCh37_NOT_DETECTED_CK.csv'.format(chrID)
        df1 = pd.read_csv(file1, low_memory=False)
        df2 = pd.read_csv(file2, low_memory=False)
        df1out = filter_cadd_overlap(df1, 37)
        df2out = filter_cadd_overlap(df2, 37)
        print("saving detected and not detected GRCh37 ", chrID)
        print()
        df1out.to_csv(out1, index=False)
        df2out.to_csv(out2, index=False)

    # GRCh38
    for order in chrls:
        chrID = 'chr{}'.format(order)
        file1 = '{}_CADD_GRCh38_DETECTED_CK.csv'.format(chrID)
        out1 = 'MISSENSE_{}_CADD_GRCh38_DETECTED_CK.csv'.format(chrID)
        file2 = '{}_CADD_GRCh38_NOT_DETECTED_CK.csv'.format(chrID)
        out2 = 'MISSENSE_{}_CADD_GRCh38_NOT_DETECTED_CK.csv'.format(chrID)
        df1 = pd.read_csv(file1, low_memory=False)
        df2 = pd.read_csv(file2, low_memory=False)
        df1out = filter_cadd_overlap(df1, 38)
        df2out = filter_cadd_overlap(df2, 38)
        print("saving detected and not detected GRCh38 ", chrID)
        print()
        df1out.to_csv(out1, index=False)
        df2out.to_csv(out2, index=False)

main()
