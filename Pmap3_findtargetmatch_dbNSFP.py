#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Pmap_findtargetmatch_dbNSFP.py

import pandas as pd
import os
from ast import literal_eval

def findtargetmatch(df):
    newdf = []
    for index, row in df.iterrows():
        newrow = []
        key19 = row['pos_id19']
        key38 = row['pos_id38']
        oaa = row['aaref']
        aaa = row['aaalt']
        aapos = row['matched_aapos']
        vep = row['HGVSp_VEP']
        cadd = row['CADD_phred_hg38']
        ukbid = row['matched_UKBID']
        posdic = row['pos_dict']
        newrow.append(key19)  # 0
        newrow.append(key38)  # 1
        newrow.append(oaa)   # 2
        newrow.append(aaa)   # 3
        newrow.append(aapos)  # 4
        newrow.append(vep)  # 5
        newrow.append(cadd)  # 6
        newrow.append(ukbid)  # 7
        newrow.append(posdic)  # 8

        python_dic = literal_eval(posdic)
        matched = []

        iv = int(aapos)
        if iv in python_dic.keys():
            dicval = python_dic[iv]
            if dicval == oaa:
                matched.append(oaa + str(iv))

        # if matched is empty
        if not matched:
            newrow.append("False")
        else:
            newrow.append(matched)  # 9
        newdf.append(newrow)

    header = ['pos_id19', 'pos_id38', 'aaref', 'aaalt', 'matched_aapos', 'HGVSp_VEP', 'CADD_phred_hg38', 'matched_UKBID', 'pos_dict', 'matched_target']
    final = pd.DataFrame(newdf)
    final.columns = header
    return final


def main():
    ref = "CYSLYS_everlabeled_groupby_entry_dictOfPositions_TODO0_3991.csv"
    refdf = pd.read_csv(ref)
    refdf.columns = ['matched_UKBID', 'labeled_pos_count', 'pos_dict']
    ref = refdf[['matched_UKBID', 'pos_dict']].copy()

    # read in each dbNSFP ckfiltered file and add key ids for 19 and 38
    infiles = ['ckfiltered_chr1.csv',
               'ckfiltered_chr2.csv',
               'ckfiltered_chr3.csv',
               'ckfiltered_chr4.csv',
               'ckfiltered_chr5.csv',
               'ckfiltered_chr6.csv',
               'ckfiltered_chr7.csv',
               'ckfiltered_chr8.csv',
               'ckfiltered_chr9.csv',
               'ckfiltered_chr10.csv',
               'ckfiltered_chr11.csv',
               'ckfiltered_chr12.csv',
               'ckfiltered_chr13.csv',
               'ckfiltered_chr14.csv',
               'ckfiltered_chr15.csv',
               'ckfiltered_chr16.csv',
               'ckfiltered_chr17.csv',
               'ckfiltered_chr18.csv',
               'ckfiltered_chr19.csv',
               'ckfiltered_chr20.csv',
               'ckfiltered_chr21.csv',
               'ckfiltered_chr22.csv',
               'ckfiltered_chrX.csv',
               'ckfiltered_chrY.csv']

    outfiles = ['cktargetmatch_chr1.csv',
                'cktargetmatch_chr2.csv',
                'cktargetmatch_chr3.csv',
                'cktargetmatch_chr4.csv',
                'cktargetmatch_chr5.csv',
                'cktargetmatch_chr6.csv',
                'cktargetmatch_chr7.csv',
                'cktargetmatch_chr8.csv',
                'cktargetmatch_chr9.csv',
                'cktargetmatch_chr10.csv',
                'cktargetmatch_chr11.csv',
                'cktargetmatch_chr12.csv',
                'cktargetmatch_chr13.csv',
                'cktargetmatch_chr14.csv',
                'cktargetmatch_chr15.csv',
                'cktargetmatch_chr16.csv',
                'cktargetmatch_chr17.csv',
                'cktargetmatch_chr18.csv',
                'cktargetmatch_chr19.csv',
                'cktargetmatch_chr20.csv',
                'cktargetmatch_chr21.csv',
                'cktargetmatch_chr22.csv',
                'cktargetmatch_chrX.csv',
                'cktargetmatch_chrY.csv']

    for ii, oo in zip(infiles, outfiles):
        filename = ii
        out = oo
        df = pd.read_csv(filename, low_memory=False,converters={'pos(1-based)':
                         '{:0>9}'.format, 'hg19_pos(1-based)': '{:0>9}'.format})
        # create the keys on the large complete df
        df.loc[:, 'pos_id38'] = df['#chr'].astype(str) + '_' + \
            df['pos(1-based)'].astype(str) + '_' + df["ref"].astype(str) + \
            '_' + df["alt"].astype(str)  # keyid 38
        df.loc[:, 'pos_id19'] = df['hg19_chr'].astype(str) + '_' + \
            df['hg19_pos(1-based)'].astype(str) + '_' + df["ref"].astype(str) + \
            '_' + df["alt"].astype(str)  # keyid 19

        # simplify df
        simdf = df[['aaref', 'aaalt', 'HGVSp_VEP', 'CADD_phred', 'matched_UKBID', 'matched_aapos', 'pos_id38', 'pos_id19']].copy()
        simdf.columns = ['aaref', 'aaalt', 'HGVSp_VEP', 'CADD_phred_hg38', 'matched_UKBID', 'matched_aapos','pos_id38','pos_id19']

        # merging simplified dataframe with uniprot position dict
        mer = pd.merge(simdf, ref, how='inner', on=['matched_UKBID'])
        print('shape simplified dbNSFP file : ', simdf.shape)
        print('shape of merged file: ', mer.shape)
        # search for labeled position
        founddf = findtargetmatch(mer)
        founddf.to_csv(out, index=False)
        df.to_csv("addedkey_"+ii, index=False)
        print("done with : ", out)
        print()
        print()
main()