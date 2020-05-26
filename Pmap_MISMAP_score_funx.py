#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# adapted code from markdown F

# 2 versions!

# 1 for all ENSP, 1 prosequence per row
# df85 = mismap_score_allprosequence(v85, UKB)

# 1 for unique ENSP, list of prosequence per row (see markdown J)
# df85 = mismap_score_unique_prosequences(v85, UKB)


'''
VERSION 1
shape ENSEMBL df:  (10183, 21)
shape ukb df:  (3953, 6)
shape post merge:  (10183, 26)
final df shape:  (10183, 34)
---
shape ENSEMBL df:  (10395, 21)
shape ukb df:  (3953, 6)
shape post merge:  (10395, 26)
final df shape:  (10395, 34)
---
shape ENSEMBL df:  (10612, 21)
shape ukb df:  (3953, 6)
shape post merge:  (10612, 26)
final df shape:  (10612, 34)
---
shape ENSEMBL df:  (10663, 21)
shape ukb df:  (3953, 6)
shape post merge:  (10663, 26)
final df shape:  (10663, 34)
---
shape ENSEMBL df:  (10564, 21)
shape ukb df:  (3953, 6)
shape post merge:  (10564, 26)
final df shape:  (10564, 34)
'''

'''
# VERSION 2
shape ENSEMBL df:  (3953, 21)
shape ukb df:  (3953, 6)
shape post merge:  (3953, 26)
final df shape:  (3953, 34)
---
shape ENSEMBL df:  (3953, 21)
shape ukb df:  (3953, 6)
shape post merge:  (3953, 26)
final df shape:  (3953, 34)
---
shape ENSEMBL df:  (3953, 21)
shape ukb df:  (3953, 6)
shape post merge:  (3953, 26)
final df shape:  (3953, 34)
---
shape ENSEMBL df:  (3953, 21)
shape ukb df:  (3953, 6)
shape post merge:  (3953, 26)
final df shape:  (3953, 34)
---
shape ENSEMBL df:  (3953, 21)
shape ukb df:  (3953, 6)
shape post merge:  (3953, 26)
final df shape:  (3953, 34)
'''

import os
import sys
import numpy as np
import pandas as pd
from ast import literal_eval

def list2col(s, colname):
    fseries = pd.DataFrame(np.array(s).reshape(-1, 1))
    fseries.columns = [colname]
    return fseries


def mismap_score_allprosequence(enspdf, uniprot):
    # new columns to add
    count_col = []
    count_cys_col = []
    count_lys_col = []
    miss_col = []
    miss_count_c = []
    miss_count_k = []
    corrFrac_col = []
    missFrac_col = []
    # merge two df
    mer = pd.merge(enspdf, uniprot, how='inner', on=['ID'])
    # QC merge
    print('shape ENSEMBL df: ', enspdf.shape)
    print('shape ukb df: ', uniprot.shape)
    print('shape post merge: ', mer.shape)
    for index, row in mer.iterrows():
        ensp_len = row['Length']
        pep = row['proSequence']
        ukb_id = row['ID']
        tot_tar = row['labeled_pos_count']
        posCK = row['pos_dict']
        tot_c = int(row['count_C_targets'])
        tot_k = int(row['count_K_targets'])
        # positions found counter
        count = 0
        count_cys = 0
        count_lys = 0
        # evaluate as dictionary
        python_dict = literal_eval(posCK)
        # loop thru position dict
        for key in python_dict:
            k = int(key)
            # convert pos # to index 0-based
            i = k - 1
            # if index within ensembl protein string
            if i < int(ensp_len):
                # get AA from peptide string
                AA = pep[i]
                # get AA from pos dic
                checker = python_dict[key]
                # match, increase count
                if AA == checker:
                    count += 1
                    if AA == 'C':
                        count_cys += 1
                    if AA == 'K':
                        count_lys += 1
        # fraction correct pos
        cor_per = count/tot_tar
        # fraction missed
        miss = tot_tar - count
        missC = tot_c - count_cys
        missK = tot_k - count_lys
        miss_per = miss/tot_tar
        # add columns
        count_col.append(count)
        count_cys_col.append(count_cys)
        count_lys_col.append(count_lys)
        miss_col.append(miss)
        miss_count_c.append(missC)
        miss_count_k.append(missK)
        corrFrac_col.append(round(cor_per, 3))
        missFrac_col.append(round(miss_per, 3))
    # add columns
    count_col2 = list2col(count_col, "found_count")
    count_cys_col2 = list2col(count_cys_col, "found_count_C")
    count_lys_col2 = list2col(count_lys_col, "found_count_K")
    miss_col2 = list2col(miss_col, "missed_count")
    miss_count_c2 = list2col(miss_count_c, "missed_count_C")
    miss_count_k2 = list2col(miss_count_k, "missed_count_K")
    corrFrac_col2 = list2col(corrFrac_col, "correct_frac")
    missFrac_col2 = list2col(missFrac_col, "missed_frac")
    df = pd.concat([mer, count_col2,count_cys_col2, count_lys_col2, miss_col2,miss_count_c2, miss_count_k2, corrFrac_col2, missFrac_col2], axis=1)
    print("final df shape: ", df.shape)
    print()
    return df



def mismap_score_unique_prosequences(enspdf, uniprot):
    # new columns to add
    count_col = []
    count_cys_col = []
    count_lys_col = []
    miss_col = []
    miss_count_c = []
    miss_count_k = []
    corrFrac_col = []
    missFrac_col = []
    # merge two df
    mer = pd.merge(enspdf, uniprot, how='inner', on=['ID'])
    # QC merge
    print('shape ENSEMBL df: ', enspdf.shape)
    print('shape ukb df: ', uniprot.shape)
    print('shape post merge: ', mer.shape)
    for index, row in mer.iterrows():
        # pep is a list of seqeunces
        pep = row['proSequence']
        total_proseq = int(row['count_proSequence']) # number of unique ensp prosequences
        ukb_id = row['ID']
        posCK = row['pos_dict']
        # multiplying # per ukb ID by number of ENSP sequences mapping to UKB ID
        tot_tar = int(row['labeled_pos_count']) * total_proseq
        tot_c = int(row['count_C_targets']) * total_proseq
        tot_k = int(row['count_K_targets']) * total_proseq
        # positions found counter, reset with new row/UKBID
        count = 0
        count_cys = 0
        count_lys = 0
        # evaluate as dictionary
        python_dict = literal_eval(posCK)
        # looping thru prosequence list
        for proseq in pep:
            ensp_len = len(proseq)
            # loop thru position dict
            for key in python_dict:
                k = int(key)
                # convert pos # to index 0-based
                i = k - 1
                # if index within ensembl protein string
                if i < int(ensp_len):
                    # get AA from peptide string
                    AA = proseq[i]
                    # get AA from pos dic
                    checker = python_dict[key]
                    # match, increase count
                    if AA == checker:
                        count += 1
                        if AA == 'C':
                            count_cys += 1
                        if AA == 'K':
                            count_lys += 1
        # fraction correct pos
        cor_per = count/tot_tar  # frac corr is total found / # pos in ukb IDs * # ensp seq checked
        # fraction missed
        miss = tot_tar - count
        missC = tot_c - count_cys
        missK = tot_k - count_lys
        miss_per = miss/tot_tar
        # add columns
        count_col.append(count)
        count_cys_col.append(count_cys)
        count_lys_col.append(count_lys)
        miss_col.append(miss)
        miss_count_c.append(missC)
        miss_count_k.append(missK)
        corrFrac_col.append(round(cor_per, 3))
        missFrac_col.append(round(miss_per, 3))
    # add columns
    count_col2 = list2col(count_col, "found_count")
    count_cys_col2 = list2col(count_cys_col, "found_count_C")
    count_lys_col2 = list2col(count_lys_col, "found_count_K")
    miss_col2 = list2col(miss_col, "missed_count")
    miss_count_c2 = list2col(miss_count_c, "missed_count_C")
    miss_count_k2 = list2col(miss_count_k, "missed_count_K")
    corrFrac_col2 = list2col(corrFrac_col, "correct_frac")
    missFrac_col2 = list2col(missFrac_col, "missed_frac")
    df = pd.concat([mer, count_col2,count_cys_col2, count_lys_col2, miss_col2,miss_count_c2, miss_count_k2, corrFrac_col2, missFrac_col2], axis=1)
    print("final df shape: ", df.shape)
    print()
    return df

