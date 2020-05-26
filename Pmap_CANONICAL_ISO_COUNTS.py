#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Pmap_CANONICAL_ISO_COUNTS.py
# 9/23/19

import os
import sys
import numpy as np
import pandas as pd
import csv
import argparse
import Bio
from Bio import SeqIO
from ast import literal_eval 
import difflib 
from statistics import mean
sys.path.append('/Users/mariapalafox/Desktop/Toolbox')
from all_funx import *
os.chdir("/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/CANONICAL_number/")

# transform uniprot canonical + isofrom KBSP human proteome fasta file into csv
def fastaToDF(filename,dbtype): 
    identifiers = [seq_record.id for seq_record in SeqIO.parse(filename, "fasta")]
    lengths = [len(seq_record.seq) for seq_record in SeqIO.parse(filename, "fasta")]
    proSeq = [seq_record.seq for seq_record in SeqIO.parse(filename, "fasta")]
    if dbtype == "uniprot":
        splitAcc = []
        splitEntry = []
        for id in identifiers:
            splitID = id.split('|')
            acc = splitID[1]
            splitAcc.append(acc)
            entryName = splitID[2]
            splitEntry.append(entryName)
        s1 = pd.Series(splitAcc, name='ID')
        s2 = pd.Series(splitEntry, name= 'entryName')
        s3 = pd.Series(lengths, name='Length')
        s4 = pd.Series(proSeq, name='proSequence')
        series = [s1,s2,s3,s4]
        df = pd.concat(series, axis=1)
        return(df)

# split ID with isoform details by "-" and group by stable ID
def splitUniprotID(df, colname): 
    proIDsplit = df[colname].str.split("-", n=1, expand=True)
    proIDsplit.columns = ['stable_entry', 'isoformNumber']
    dfmerge = pd.concat([df,proIDsplit], axis=1)
    # converting canonical seq with no isofrom into 0 value
    dfmerge.isoformNumber.fillna(value=0, inplace=True)
    # converting isoform column into categorical
    dfmerge.isoformNumber = dfmerge.isoformNumber.astype('object')
    # idea, make the isoformnum and len in tuple, so i can compare if canonical 0 is -1 in ID mapping file
    stable_group = dfmerge.groupby(['stable_entry'])['isoformNumber'].apply(list).reset_index()
    return stable_group

# create list from stable_entry df for each column
def grabTrue_df(df, colname):
    dfTrue = df[df[colname] == 'True']
    dfTrue.reset_index(drop=True, inplace=True)
    return dfTrue

def main():
    #dffasta = "UniprotKB_isoforms_canonical_filtered_homosapiens_9606_2018_06release.fasta"
    #dfuni = fastaToDF(dffasta, "uniprot")
    #dfuni.to_csv("fasta2csv_iso_canon_ukb_2018.csv",index=False)
    dfuni=pd.read_csv("fasta2csv_iso_canon_ukb_2018.csv", index=False)
    # total lines = 195341
    # extract column with ID and Sequence only
    dfukb = dfuni[['ID','Length']].copy()
    # parse ID column in proteome file
    stable_fasta = splitUniprotID(dfukb, "ID")
    describeMe(stable_fasta)
    # read in IDdat file
    idat = pd.read_table("HUMAN_9606_idmapping_column1_unique.txt", header=None)
    idat.columns = ['ID']
    # parse ID column in proteome file
    stable_idat = splitUniprotID(idat, "ID")
    describeMe(stable_idat)
    # add additional column if uniprot ID are in ccds set
    stable_idat.sort_values(by=['stable_entry'], inplace = True)
    stable_fasta.sort_values(by=['stable_entry'], inplace = True)

    # changing column names preparing for merge of all rows and subset ccds rows
    fasta = stable_fasta[['stable_entry','isoformNumber']].copy()
    fasta.columns = ['stable_entry', 'fasta_isoformNumber']
    idat = stable_idat[['stable_entry','isoformNumber']].copy()
    idat.columns = ['stable_entry', 'idat_isoformNumber']

    # creating merge on stable entry
    merge170k = pd.merge(idat,fasta,how='inner',on=['stable_entry'])
    # add length column to this merged file before parsing out TRUE ccds subset rows
    merge170k['idat_len'] = merge170k['idat_isoformNumber'].str.len()
    merge170k['fasta_len'] = merge170k['fasta_isoformNumber'].str.len()

    # SUB SET LABELING
    # reading in 18432 simple labels
    ccds = pd.read_csv("simple18432IDs.txt",header=None, names = ['ID'])
    labeled = pd.read_csv("all3991_everLabeled_entries_simple.csv")
    sharedlabeled = pd.read_csv("uniprotIDs3979.csv", header=None, names = ['ID'])
    # making lists
    lsccds = ccds.ID.tolist() # in 18432 from which come the labeled subsets
    lsccds.remove('ID') # len now 18432
    labeled.columns = ['ID']
    lslabeled = labeled.ID.tolist() # labeled before
    lssharedlabeled = sharedlabeled.ID.tolist() # shared and labeled before
    # labeling this df with the list of made to represent eachs subset
    merge170k = addcolumnconditional(lslabeled, merge170k, 'stable_entry', 'labeledSet_3991')
    merge170k = addcolumnconditional(lssharedlabeled, merge170k, 'stable_entry', 'labeledShared_3979')
    merge170k = addcolumnconditional(lsccds, merge170k, 'stable_entry', 'UKBccds_9606_set_18432')
    checkColumnValues(merge170k, 'labeledSet_3991') # all IDs present
    checkColumnValues(merge170k, 'labeledShared_3979') # all IDs present
    checkColumnValues(merge170k, 'UKBccds_9606_set_18432') # all IDs present
    merge170k.to_csv("merged_fasta_idat_can_iso_columnsFORnumber_isoforms_per_stableID_173324.csv", index=False)
    
    # saving single isoform entries
    singleIsoform = merge170k[(merge170k['idat_len'] == 1) & (merge170k['fasta_len'] ==1)]
    checkColumnValues(singleIsoform,"idat_isoformNumber" )
    checkColumnValues(singleIsoform,"fasta_isoformNumber" )
    checkColumnValues(singleIsoform,"labeledSet_3991" )
    checkColumnValues(singleIsoform,"labeledShared_3979" )
    checkColumnValues(singleIsoform,"UKBccds_9606_set_18432" )
    describeMe(singleIsoform)
    singleIsoform.to_csv("merged_fasta_idat_singleIsoform_entries_162688.csv", index=False)

    # pull out all non single Isoform rows
    multiIsoform = merge170k[merge170k['idat_len'] != 1]
    describeMe(multiIsoform) # 10635 rows
    checkColumnValues(multiIsoform,'idat_isoformNumber')
    checkColumnValues(multiIsoform,'fasta_isoformNumber')
    checkColumnValues(multiIsoform,"idat_isoformNumber")
    checkColumnValues(multiIsoform,"fasta_isoformNumber")
    checkColumnValues(multiIsoform,"labeledSet_3991")
    checkColumnValues(multiIsoform,"labeledShared_3979")
    checkColumnValues(multiIsoform,"UKBccds_9606_set_18432")
    # getting difference in Isoform numbers from idat - fasta
    multiIsoform = multiIsoform.assign(idatLen_minus_fastaLen = multiIsoform["idat_len"] - multiIsoform["fasta_len"])
    checkColumnValues(multiIsoform, 'idatLen_minus_fastaLen')
    # 10628 had 1 difference, only keep these rows
    multiIsoform = multiIsoform[multiIsoform['idatLen_minus_fastaLen']==1]
    describeMe(multiIsoform)
    checkColumnValues(multiIsoform,"labeledSet_3991")
    checkColumnValues(multiIsoform,"labeledShared_3979")
    checkColumnValues(multiIsoform,"UKBccds_9606_set_18432")
    multiIsoform.to_csv("merged_fasta_idat_multiIsoform_1moreIDinIDATfiltered_10628.csv", index=False)
    multiIsoform.reset_index(drop=True, inplace=True)

"""
QC
# FLAG CHECK FOR NUMBER OF ISOFORMS DIFF BTW IDAT - FASTA = 1
# RESULT: ALL ROWS HAD MATCHED UKB ID -1 DIFF BTW ID SOURCES

    # how many rows in idat have 0 in the list, which means i can drop it
    multiIsoform['flag0_idat'] = multiIsoform.apply(lambda x: 0 in x['idat_isoformNumber'], axis=1).astype(int)

    multiIsoform['flag0_fasta'] = multiIsoform.apply(lambda x: 0 in x['fasta_isoformNumber'], axis=1).astype(int)

    checkColumnValues(multiIsoform, 'flag0_idat')
    checkColumnValues(multiIsoform, 'flag0_fasta')

# all rows have both 0 in idat isoform numbers column and 0 in fasta column

 flag0_idat  Count
0           1  10628
   flag0_fasta  Count
0            1  10628
(DATA 10628, 11)

"""
    # FIGURE CODE for isoform count
    KBisoforms = merge170k[['stable_entry','fasta_len','labeledSet_3991','labeledShared_3979','UKBccds_9606_set_18433']].copy()
    KBisoforms.shape
    KBisoforms.to_csv("isoformCounts_UKBstable_entry_fastawIsos_173324.csv",index=False)

    sidat = multiIsoform.idat_isoformNumber
    sfas = multiIsoform.fasta_isoformNumber
    # getting canonical isoform number
    DAT = []
    FAS = []
    CAN = []
    for i in range(len(sidat)):  
        lsIdat = sidat[i]
        lsFas = sfas[i]
        DAT.append(lsIdat)
        FAS.append(lsFas)
        dif = list(set(lsIdat) - set(lsFas))
        CAN.append(dif)

    canon_col = pd.DataFrame(CAN, columns = ['canonical_iso_number'])
    final = pd.concat([multiIsoform, canon_col], axis=1)
    # shape 10628, 10
    canonicalcount = final.canonical_iso_number.value_counts()
    canonicalcount.to_csv("hist_canonical_num_counts_fasta2018_idmap2018.csv", index=False)

    failed = final[final['canonical_iso_number'] != '1']
    failed.reset_index(inplace=True, drop=True)
    failed.columns = ['stable_entry', 'idat_isoformNumber', 'fasta_isoformNumber', 'idat_len',
       'fasta_len', 'labeledSet_3991', 'labeledShared_3979',
       'UKBccds_9606_set_18432', 'idatLen_minus_fastaLen',
       'canonical_iso_number']
    failed.to_csv("canonical_NOT_1_IDs_fastavsidmap_288.csv", index=False)

    labTRUE = grabTrue_df(final, 'labeledSet_3991')
    labenspTRUE = grabTrue_df(final,'labeledShared_3979')
    # should be 18432 shape
    print(labTRUE.shape)
    print(labenspTRUE.shape)
    in3991 = labTRUE.canonical_iso_number.value_counts()
    in3979 = labenspTRUE.canonical_iso_number.value_counts()
    in91 = pd.DataFrame(in3991)
    in79 = pd.DataFrame(in3979)
    in91.reset_index(inplace=True)
    in79.reset_index(inplace=True)
    in91.columns = ['Isoform_number','Count']
    in79.columns = ['Isoform_number','Count']
    in91.to_csv("hist_in3991_2507.csv",index=False)
    in79.to_csv("hist_in3979_2504.csv", index=False)
    not1_3991 = labTRUE[labTRUE['canonical_iso_number'] != '1']
    not1_3979 = labenspTRUE[labenspTRUE['canonical_iso_number'] != '1']
    print(not1_3979.shape)
    print(not1_3991.shape)
    not1_3991.to_csv("canonical_NOT_1_IDs_in3991_56.csv",index=False)
    not1_3979.to_csv("canonical_NOT_1_IDs_in3979_56.csv",index=False)

main()

