# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd

# Pmap_dbNSFP_CADD_merge_DECT19.py

# see markdown M local for all code, this script is tailored for detected hg19 cadd positions merging with posid19 from dbNSFP col

'''
filter files: 

col_det_posid19.csv --- ALLCON_{}_CADD_GRCh37_DETECTED_CK.csv

col_notdet_posid19.csv --- ALLCON_{}_CADD_GRCh37_NOT_DETECTED_CK.csv

col_det_posid38.csv --- ALLCON_{}_CADD_GRCh38_DETECTED_CK.csv

col_notdet_posid38.csv --- ALLCON_{}_CADD_GRCh38_NOT_DETECTED_CK.csv

'''

def main():
    # CHANGE PATH NAME HERE
    os.chdir('/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/CADDmapped/RESULT_CADDv14_pos_overlap_dbNSFPcoordinates/ALL_CONSEQUENCES/37det/')

    chrls = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']

    for order in chrls:
        chrID = 'chr{}'.format(order)

        # CHANGE FILE NAMES HERE
        file1 = 'ALLCON_{}_CADD_GRCh37_DETECTED_CK.csv'.format(chrID)
        out1 = 'FILTERONposid_{}_CADD_GRCh37_DETECTED_CK.csv'.format(chrID)
        df1 = pd.read_csv(file1, low_memory=False)

        # CHANGE FILTER FILE NAME HERE
        fildf = pd.read_csv('col_det_posid19.csv')

        print("filter file shape: ")
        print(fildf.shape)
        print("chr cadd file shape: ")
        print(df1.shape)

        dropme = ['Type', 'Length', 'AnnoType','ConsScore', 'ConsDetail','oAA', 'nAA', 'GeneID','FeatureID', 'GeneName', 'CCDS', 'Intron', 'Exon', 'cDNApos','relcDNApos', 'CDSpos', 'relCDSpos', 'relProtPos', 'Domain','Dst2Splice', 'Dst2SplType', 'minDistTSS', 'minDistTSE']

        df1.drop(dropme, axis=1, inplace=True)
        df1.drop_duplicates(keep='first', inplace = True)

        # CHANGE MERGE ON HERE
        mer = pd.merge(df1, fildf, how='inner', on=['pos_id19'])
        print("merge shape : ")
        print(mer.shape)
        mer.to_csv(out1, index=False)
        print("done with ", order)
main()

