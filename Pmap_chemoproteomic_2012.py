#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on
for Project Pmap, task A
@author: mfpfox
data: 11.6.19

original 4119 chemoproteomic UKB IDs mapped to 11_2012 provided fasta file, lost 26 IDs (not found), making new ID count 4093. The following code takes the 4093 IDs found in fasta file 11_2012 and compared detected CYS and LYS positions to the original data reference.

expected all positions to be matched in 4093 protein sequences. following this task, data was mapped to 2018 canonical reference proteins.

code creates 2 types of files for each amino acid subset (pantarget or reactivity)- positions of detected CK found OR positions not all found per protein

"""
import sys
import os
import pandas as pd

os.chdir("/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/CYS_LYS_2012_to2018UKBpositions")

def split_ID(df1, col):
    new1 = df1[col].str.split("_", n=1, expand = True)
    df1["ID"] = new1[0]
    return df1

def count_unique_IDs(df, col):
    idset = set(df[col])
    print("unique IDs in df: ", len(idset))
    return idset

# aalist is from fasta sequences, chemdf is chemoproteomic data subset
def position_overlap(aalist, chemdf, found_name, notfound_name):
    # create list of reactivity and pan target positions
    subsetls = list(chemdf.pos_ID)
    print("chemo pos_ID length: ", len(subsetls))
    # overlap with positions from fasta sequences
    overlapAA = set(aalist).intersection(subsetls)
    print("overlapping AA's in chemoproteomic data subset: ", len(overlapAA))
    # make dataframe with overlap results
    newls = list(overlapAA)
    newdf = pd.DataFrame(newls)
    newdf.columns = ['pos_ID']
    # count unique ID
    splitdf = split_ID(newdf, 'pos_ID')
    idset = count_unique_IDs(splitdf, 'ID')

    # label the subset data frame with new column whether its in fasta sequence or not
    done_overlapped = addcolumnconditional(overlapAA, chemdf, "pos_ID", "found_pos_ID")

    # seperating True and False then saving False
    TRUEonly = done_overlapped[done_overlapped["found_pos_ID"] == "True"]
    FALSEonly = done_overlapped[done_overlapped["found_pos_ID"] == "False"]
    print("Found posID TRUE df shape: ", TRUEonly.shape)
    print("Found posID FALSE df shape: ", FALSEonly.shape)

    # saving dfs
    FALSEonly.reset_index(drop=True,inplace=True)
    FALSEonly.to_csv(notfound_name,index=False)
    TRUEonly.reset_index(drop=True,inplace=True)
    TRUEonly.to_csv(found_name,index=False)
    return idset

def main(): 
    # subsets filtered for 4093 detected IDs overlapping 11_2012 fasta
    cpan = pd.read_csv("CYS_pantarget_11_2012IDs_6157.csv")
    crea = pd.read_csv("CYS_reactivity_11_2012IDs_1483.csv")
    kpan = pd.read_csv("LYS_pantarget_11_2012IDs_8350.csv")
    krea = pd.read_csv("LYS_reactivity_11_2012IDs_4576.csv")

    # 11_2012 release with 4093 detected IDs
    #ukbfasta = pd.read_csv("sprot_11_2012_detectedIDs_4093.csv")
    fastaC = pd.read_csv("sprot_11_2012_all_Cpos.csv")
    fastaK = pd.read_csv("sprot_11_2012_all_Kpos.csv")
    # make into list with only pos_ID
    allC = list(fastaC.key_id)
    allK = list(fastaK.key_id)

    print()
    CRidset = position_overlap(allC, crea, "pos_found_C_react.csv", "pos_not_found_C_react.csv")
    print()
    CPidset = position_overlap(allC, cpan, "pos_found_C_pan.csv", "pos_not_found_C_pan.csv")
    print()
    KRidset = position_overlap(allK, krea, "pos_found_K_react.csv", "pos_not_found_K_react.csv")
    print()
    KPidset = position_overlap(allK, kpan, "pos_found_K_pan.csv", "pos_not_found_K_pan.csv")

    totCid = set(list(CRidset) + list(CPidset))
    print("total Cys group IDs: ", len(totCid))
    totKid = set(list(KRidset) + list(KPidset))
    print("total Lye group IDs: ", len(totKid))

    allIDs = set(list(totCid) + list(totKid))
    print(len(allIDs))
