#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 11/14/19 Pmap_protein_sequence_comparison.py

"""
- input files must have  ID and proSequence column (fasta -> CSV output files)
- compares the prosequence in df against ID:seqeunce dict created, adds column True for identical False for not identical
- seqeunceDifference funx takes only False identity rows and identifies WHAT is different about them
    - added code for positions (1 based) that are different
- made from MARKDOWN 'B'

# COLUMN NAME NOTES: 
- $difference_18vs12 = the amino acid in 2018 sequence that is not in 2012 sequence
- $posdiff_ref_shorter_pro = the position (1 based not index) of difference between 2 prosequences
"""

import os
import sys
import pandas as pd
from ast import literal_eval
import difflib

sys.path.append("/Users/mariapalafox/Desktop/Toolbox")
from all_funx import *



def identicalSequenceCheck(dfref, dfalt):
    ref_dic = dict(zip(dfref.ID, dfref.proSequence))
    ref_len = dict(zip(dfref.ID, dfref.Length))
    newcol = []
    len2012 = []
    for index, row in dfalt.iterrows():
        pep = row['proSequence']
        ukb_id = row['ID']
        # retrieve uniprot fasta seq using reference dictionary
        mypep = ref_dic[ukb_id]
        len12 = ref_len[ukb_id]
        str(mypep)
        str(pep)
        if mypep == pep:
            newcol.append('True')
            len2012.append(len12)
        if mypep != pep: 
            newcol.append('False')
            len2012.append(len12)
    # add identity score list as new column
    dfalt.loc[:,'identical_seq'] = newcol
    dfalt.loc[:,'length_2012_fasta'] = len2012
    dfalt.loc[:,'len12_minus_len18'] = dfalt['length_2012_fasta'] - dfalt['Length']
    print(dfalt.shape)
    checkColumnValues(dfalt, "identical_seq")
    return dfalt

# assumes ID and proSequence column
# $$$
# dfref is all 2012 seqeucne dfalt is only non identical sequences 2018
def sequenceDifference(dfref, dfalt):
    ref_dic = dict(zip(dfref.ID, dfref.proSequence))
    diffcol = []
    posdiffcol = []
    for index, row in dfalt.iterrows():
        altpep = row['proSequence']
        ukb_id = row['ID']
        # retrieve 2012seq from dict
        refpep = ref_dic[ukb_id]
        str(refpep)
        str(altpep)
        lenref = len(refpep)
        lenalt = len(altpep)

        # len check is for posout_list2, shorter pep must be in range(len)
        if lenref == lenalt:
            output_list = [altpep[i] for i in range(len(altpep)) if altpep[i] != refpep[i]]
            posoutput_list2 = [i+1 for i in range(len(altpep)) if altpep[i] != refpep[i]]
            out = ', '.join(output_list)
            #posout = ', '.join(posoutput_list2)
            diffcol.append(out)
            posdiffcol.append(posoutput_list2)

        # 2012 protein longer than 2018
        if lenref > lenalt:
            output_list = [altpep[i] for i in range(len(altpep)) if altpep[i] != refpep[i]]
            posoutput_list2 = [i+1 for i in range(len(altpep)) if altpep[i] != refpep[i]]
            out = ', '.join(output_list)
            #posout = ', '.join(posoutput_list2)
            diffcol.append(out)
            posdiffcol.append(posoutput_list2)

        # 2018 protein longer than 2012
        if lenref < lenalt:
            output_list = [altpep[i] for i in range(len(refpep)) if refpep[i] != altpep[i]]
            posoutput_list2 = [i+1 for i in range(len(refpep)) if refpep[i] != altpep[i]]
            out = ', '.join(output_list)
            #posout = ', '.join(posoutput_list2)
            diffcol.append(out)
            posdiffcol.append(posoutput_list2)

    # add diff as new column
    dfalt.loc[:, 'difference_18vs12'] = diffcol
    # add pos diff as new column
    dfalt.loc[:, 'posdiff_ref_shorter_pro'] = posdiffcol
    print("new dfalt shape with added columns: ", dfalt.shape)
    return dfalt


def main():
    # file with ID col for IDs you want to only compare sequences of
    ids = pd.read_csv("UKBIDs_2018_3964.csv")
    fasta12 = pd.read_csv("cravattlab_11_2012.csv") # what paper says is ref
    fasta18 = pd.read_csv("UKB_CCDS_2018_06_fasta.csv")
    print("ID file shape: ", ids.shape)
    print("fasta 2012 shape: ", fasta12.shape)
    print("fasta 2018 shape: ", fasta18.shape)
    print()

    # only keep IDs/sequences from both fastas that match an ID in IDlist
    idlist = list(ids.ID)
    seq2012 = addcolumnconditionalDrop(idlist, fasta12, 'ID', 'in_IDlist')
    seq2018 = addcolumnconditionalDrop(idlist, fasta18, 'ID', 'in_IDlist')

    # call funx that compares sequences
    checked2018 = identicalSequenceCheck(seq2012, seq2018)
    # saving all 2018 sequences in ID list with TRUE FALSE seq comparison col
    checked2018.to_csv("UKBref2.0_2018_T3931_F33.csv", index=False)

    # isolate FALSE sequence identity rows
    falseonly = checked2018[checked2018['identical_seq'] == 'False']

    # $$$ outputs what exactly is different about seqeunces
    diff2018 = sequenceDifference(seq2012, falseonly)
    diff2018.to_csv("UKBseq2.0_diff_33nonidentical_proSeq_2012vs2018seq.csv", index=False)
