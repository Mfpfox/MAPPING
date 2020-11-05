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

#os.chdir("/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/CYS_LYS_2012_to2018UKBpositions")

"""
import sys
import os
import pandas as pd
import argparse
from maplib import *

parser = argparse.ArgumentParser()
parser.add_argument("fasta", help="FASTA to search list of CpDAA keys against, accepts UniProtKB, Gencode, and RefSeq files.")
parser.add_argument("aalist", help=".csv file with column of CpDAA keys to search for in FASTA file sequences (e.g. 'P11313_K170').")
parser.add_argument("outfile", help="name of excel file output with tabs for Found and NotFound CpDAA keys in reference FASTA (e.g. 'Cys_isoTOP2012')")

args = parser.parse_args()
fasta  = args.fasta
aalist = args.aalist
outfile = args.outfile

# aalist is from fasta sequences, chemdf is chemoproteomic data subset
def CpDAA_position_check(): 

    # import ukb fasta file


    # convert to .csv file
    

    # make df of all AA pos keys (e.g. P11413_K170) input AA
    


    fastaC = pd.read_csv("sprot_11_2012_all_Cpos.csv")

    #  all AA keys as list
    aalist = list(fastaC.key_id)

    ## file with CpDAA keys, column name assumed as "pos_ID"
    chemdf = pd.read_csv("CYS_reactivity_11_2012IDs_1483.csv")

# save as multi tab excel sheet instead
    found_name = "pos_found_C_react.csv"
    notfound_name = "pos_not_found_C_react.csv" 

    # create list of detected keys
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

    # labels chemoproteomic results df w/ col for whether key found in reference fasta file sequences
    done_overlapped = addcolumnconditional(overlapAA, chemdf, "pos_ID", "found_pos_ID")

    # seperating True and False then saving as tabs of excel file
    TRUEonly = done_overlapped[done_overlapped["found_pos_ID"] == "True"]
    FALSEonly = done_overlapped[done_overlapped["found_pos_ID"] == "False"]
    print("Found pos_ID keys in fasta reference seq: ", TRUEonly.shape)
    print("Not found pos_ID keys in fasta reference seq: ", FALSEonly.shape)

    # saving dfs
    FALSEonly.reset_index(drop=True,inplace=True)
    FALSEonly.to_csv(notfound_name,index=False)
    TRUEonly.reset_index(drop=True,inplace=True)
    TRUEonly.to_csv(found_name,index=False)

    totCid = set(list(idset))
    print("total Cys group IDs: ", len(totCid))


if __name__ == '__main__':
    CpDAA_position_check()