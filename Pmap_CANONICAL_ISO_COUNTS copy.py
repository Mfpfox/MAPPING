#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Pmap_CANONICAL_ISO_COUNTS.py
# 9/23/19

# packages
import os
import sys
import numpy as np
import pandas as pd
import csv
import argparse
import Bio
from Bio import SeqIO
from ast import literal_eval # for mismap_score func
import difflib # diff btw 2 strings
from statistics import mean

sys.path.append('/Users/mariapalafox/Desktop/Toolbox')
from all_funx import *

# setting path as CANONICAL project folder
os.chdir("/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/CANONICAL_number/")

# transform uniprot canonical + isofrom KBSP human proteome fasta file into csv
def fastaToDF(filename,dbtype): 
    # parse sequence fasta file
    identifiers = [seq_record.id for seq_record in SeqIO.parse(filename, "fasta")]
    lengths = [len(seq_record.seq) for seq_record in SeqIO.parse(filename, "fasta")]
    proSeq = [seq_record.seq for seq_record in SeqIO.parse(filename, "fasta")]
    # splitting identifiers into 2 seperate list
    if dbtype == "uniprot":
        splitAcc = []
        splitEntry = []
        for id in identifiers:
            splitID = id.split('|')
            acc = splitID[1]
            splitAcc.append(acc)
            entryName = splitID[2]
            splitEntry.append(entryName)
        #converting lists to pandas Series
        s1 = pd.Series(splitAcc, name='ID')
        s2 = pd.Series(splitEntry, name= 'entryName')
        s3 = pd.Series(lengths, name='Length')
        s4 = pd.Series(proSeq, name='proSequence')
        series = [s1,s2,s3,s4]
        df = pd.concat(series, axis=1)
        return(df)

# split by -, group by stable ID
# splitUniprotID(dfukb, "ID")
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


# does this merger df solve the mystery of how many canonical sequences are -1 version ?
def main():
    dffasta = "UniprotKB_isoforms_canonical_filtered_homosapiens_9606_2018_06release.fasta"
    dfuni = fastaToDF(dffasta, "uniprot")
    #dfuni.to_csv("fasta2csv_iso_canon_ukb_2018.csv",index=False)
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

    # reading in 18432 simple labels
    ccds = pd.read_csv("simple18432IDs.txt",header=None, names = ['ID'])
    labeled = pd.read_csv("all3991_everLabeled_entries_simple.csv")
    sharedlabeled = pd.read_csv("uniprotIDs3979.csv", header=None, names = ['ID'])
    # making lists
    lsccds = ccds.ID.tolist() # in 18433 from which come the labeled subsets
    lsccds.remove('ID') # len now 18432
    labeled.columns = ['ID']
    lslabeled = labeled.ID.tolist() # labeled before
    lssharedlabeled = sharedlabeled.ID.tolist() # shared and labeled before
    # labeling this df with the list of made to represent eachs subset
    merge170k = addcolumnconditional(lslabeled, merge170k, 'stable_entry', 'labeledSet_3991')
    merge170k = addcolumnconditional(lssharedlabeled, merge170k, 'stable_entry', 'labeledShared_3979')
    merge170k = addcolumnconditional(lsccds, merge170k, 'stable_entry', 'UKBccds_9606_set_18433')
    checkColumnValues(merge170k, 'labeledSet_3991') # all IDs present
    checkColumnValues(merge170k, 'labeledShared_3979') # all IDs present
    checkColumnValues(merge170k, 'UKBccds_9606_set_18432') # all IDs present
    merge170k.to_csv("merged_fasta_idat_can_iso_columnsFORnumber_isoforms_per_stableID_173324.csv", index=False)
    # ----------------------------------------------------------------------
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
"""
    # how many rows in idat have 0 in the list, which means i can drop it
    multiIsoform['flag0_idat'] = multiIsoform.apply(lambda x: 0 in x['idat_isoformNumber'], axis=1).astype(int)
    multiIsoform['flag0_fasta'] = multiIsoform.apply(lambda x: 0 in x['fasta_isoformNumber'], axis=1).astype(int)
    checkColumnValues(multiIsoform, 'flag0_idat')
    checkColumnValues(multiIsoform, 'flag0_fasta')
    # all rows have both 0 in idat isoform numbers column and 0 in fasta column
"""
    idatTRUE = grabTrue_df(stable_idat, 'UKBccds_9606_set_18432')
    fastaTRUE = grabTrue_df(stable_fasta, 'UKBccds_9606_set_18432')
    fastaTRUE = fastaTRUE[['stable_entry','isoformNumber']].copy()
    fastaTRUE.columns = ['stable_entry', 'fasta_isoformNumber']
    idatTRUE = idatTRUE[['stable_entry','isoformNumber']].copy()
    idatTRUE.columns = ['stable_entry', 'idat_isoformNumber']
    mergeIDmapping = pd.merge(idatTRUE,fastaTRUE,how='inner',on=['stable_entry'])
    # all labeled uniprot IDs
    all3991 = "all3991_everLabeled_entries_simple.csv"
    shared3979 = "uniprotIDs3979.csv"
    labeled = pd.read_csv(all3991)
    sharedlabeled = pd.read_csv(shared3979, header=None, names = ['ID'])
    labeled.columns = ['ID']
    lslabeled = labeled.ID.tolist() # labeled before
    lssharedlabeled = sharedlabeled.ID.tolist() # shared and labeled before
    # labeling this df with the list of made to represent eachs subset
    addcolumnconditional(lslabeled, mergeIDmapping, 'stable_entry', 
        'labeledSet_3991')
    addcolumnconditional(lssharedlabeled, mergeIDmapping, 'stable_entry', 
        'labeledShared_3979')
    checkColumnValues(mergeIDmapping, 'labeledSet_3991')
    checkColumnValues(mergeIDmapping, 'labeledShared_3979')
    # saved file on jupyter

    # read merge file in
    # add len of list columns 
    # cp idat file and remove 0 for rows greater than 1 len
    # compare missing number ...the idat is 1 more in len but after remving 0 they should be same len
    # if diff in list is 0 and 1 than 1 is canonical and so forth

    # add len of list column
    # find rows with 1 in fasta file for rows with more than 1 isonumber
    # greater than len 1 rows look for canonical counts



main()

