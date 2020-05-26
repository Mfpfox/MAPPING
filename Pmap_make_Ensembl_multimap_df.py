#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
creates merge of xref files (containing stable IDs and uniprot accessions for each release) with ensembl fasta files processed to contain stable and version ensembl IDs.

from merged file i can calculate stable ID multimapping to set of shared uniprot IDs. 

final task was using stable ID key added for merge (ENSP_ENST_ENSG) to create version of merge ensembl files that contain shared stable IDs, dfs will be used to calculate version ID differences controlling for set of stable IDs.

'''

import pandas as pd
import os
import sys
sys.path.append("/Users/mariapalafox/Desktop/Toolbox/")
from all_funx import *


# drops rows where ID does not match list IDs 'Lever'
def dropNotLabeled(df, sharedls, colname):
    df['Shared_IDs'] = np.where(df[colname].isin(sharedls), "True", "False")
    df2 = df[df['Shared_IDs'] == "True"].copy()
    print("sorting ID column and reseting index: ")
    df2.sort_values(by=[colname], inplace=True)
    df2.reset_index(drop=True, inplace=True)
    print("shape final cleaned df: ", df2.shape)
    return df2

def stableID_key(df, option):
    if option == 'xref':
        df['stableID_key'] = df['gene_stable_id'].astype(str) + \
'_' + df['transcript_stable_id'].astype(str) + \
'_' + df['protein_stable_id'].astype(str)
    if option == 'fasta':
        df['stableID_key'] = df['ENSG'].astype(str) + \
'_' + df['ENST'].astype(str) + \
'_' + df['ENSP'].astype(str)


def sort_col(df2):
    df2.sort_values(by=['xref', 'stableID_key'],inplace=True)
    df2.reset_index(drop=True, inplace=True)
    return df2

def shared_IDs(df1, df2, df3, df4, df5, colname):
    # SHARED SET OF IDS
    # unique set of IDs in each release
    uni85 = set(df1[colname].tolist())
    uni92 = set(df2[colname].tolist())
    uni94 = set(df3[colname].tolist())
    uni96 = set(df4[colname].tolist())
    uni97 = set(df5[colname].tolist())
    # IDs in all releases
    shared = uni85 & uni92 & uni94 & uni96 & uni97
    print("set of shared IDs: ")
    print(len(set(shared)))
    shared = list(shared)
    # filtering release tsv files for only shared uniprot IDs
    shared85 = dropNotLabeled(df1, shared, colname)
    shared92 = dropNotLabeled(df2, shared, colname)
    shared94 = dropNotLabeled(df3, shared, colname)
    shared96 = dropNotLabeled(df4, shared, colname)
    shared97 = dropNotLabeled(df5, shared, colname)
    uniqueCount(shared85, colname)
    uniqueCount(shared92, colname)
    uniqueCount(shared94, colname)
    uniqueCount(shared96, colname)
    uniqueCount(shared97, colname)
    return shared85, shared92, shared94, shared96, shared97


def shared_uniprot_IDs():
    # reference 2018 uniprot IDs
    ID = pd.read_csv("UKBIDs_2018_3963.csv")
    # making ID list to file
    idls = list(ID.ID)
    print(len(idls))
    # ensembl cross-ref files for UKB stable --- ENSP ENST ENSG for releases
    t85 = pd.read_csv("Homo_sapiens.GRCh37.85.uniprot.tsv", delimiter=r"\s+")
    t92 = pd.read_csv("Homo_sapiens.GRCh38.92.uniprot.tsv", delimiter=r"\s+")
    t94 = pd.read_csv("Homo_sapiens.GRCh38.94.uniprot.tsv", delimiter=r"\s+")
    t96 = pd.read_csv("Homo_sapiens.GRCh38.96.uniprot.tsv", delimiter=r"\s+")
    t97 = pd.read_csv("Homo_sapiens.GRCh38.97.uniprot.tsv", delimiter=r"\s+")

    found85 = addcolumnconditionalDrop(idls, t85, 'xref', 'found_3963set')
    found92 = addcolumnconditionalDrop(idls, t92, 'xref', 'found_3963set')
    found94 = addcolumnconditionalDrop(idls, t94, 'xref', 'found_3963set')
    found96 = addcolumnconditionalDrop(idls, t96, 'xref', 'found_3963set')
    found97 = addcolumnconditionalDrop(idls, t97, 'xref', 'found_3963set')
    print()
    print("How many unique IDs from 3963 set in xref files?")
    uniqueCount(found85, 'xref')
    uniqueCount(found92, 'xref')
    uniqueCount(found94, 'xref')
    uniqueCount(found96, 'xref')
    uniqueCount(found97, 'xref')

# ------------------------------------------------
    shared85, shared92, shared94, shared96, shared97 = shared_IDs(found85, found92, found94, found96, found97, 'xref')

# ------------------------------------------
    stableID_key(shared85, 'xref')
    stableID_key(shared92, 'xref')
    stableID_key(shared94, 'xref')
    stableID_key(shared96, 'xref')
    stableID_key(shared97, 'xref')
    # all xref files have unique ENSG_ENST_ENSP values, ready to map to fasta files
    # shared85.to_csv("xref_v85_shared3955_10188.csv", index=False)
    # shared92.to_csv("xref_v92_shared3955_10399.csv", index=False)
    # shared94.to_csv("xref_v94_shared3955_10616.csv", index=False)
    # shared96.to_csv("xref_v96_shared3955_10667.csv", index=False)
    # shared97.to_csv("xref_v97_shared3955_10568.csv", index=False)

    # simplify xref files
    shared85 = shared85[['xref', 'stableID_key']].copy()
    shared92 = shared92[['xref', 'stableID_key']].copy()
    shared94 = shared94[['xref', 'stableID_key']].copy()
    shared96 = shared96[['xref', 'stableID_key']].copy()
    shared97 = shared97[['xref', 'stableID_key']].copy()

    # -------------------FASTA MERGE----------------------------------
    # fasta to csv files with stable ID added by script
    v85 = pd.read_csv("v85Homo_sapiens.GRCh37.pep.all.csv")
    v92 = pd.read_csv("v92Homo_sapiens.GRCh38.pep.all.csv")
    v94 = pd.read_csv("v94Homo_sapiens.GRCh38.pep.all.csv")
    v96 = pd.read_csv("v96Homo_sapiens.GRCh38.pep.all.csv")
    v97 = pd.read_csv("v97Homo_sapiens.GRCh38.pep.all.csv")
    # add stableID_key
    stableID_key(v85, 'fasta')
    stableID_key(v92, 'fasta')
    stableID_key(v94, 'fasta')
    stableID_key(v96, 'fasta')
    stableID_key(v97, 'fasta')
    # merge on stable id to append uniprot ID with ensembl stable + version IDs
    mer85 = pd.merge(v85,shared85,how='inner',on=['stableID_key'])
    mer92 = pd.merge(v92,shared92,how='inner',on=['stableID_key'])
    mer94 = pd.merge(v94,shared94,how='inner',on=['stableID_key'])
    mer96 = pd.merge(v96,shared96,how='inner',on=['stableID_key'])
    mer97 = pd.merge(v97,shared97,how='inner',on=['stableID_key'])
    # confirm all 3955 IDs mapped to fasta files
    uniqueCount(mer85, 'xref')
    uniqueCount(mer92, 'xref')
    uniqueCount(mer94, 'xref')
    uniqueCount(mer96, 'xref')
    uniqueCount(mer97, 'xref')

    # sort on xref and keyID then saving files
    mer85 = sort_col(mer85)
    mer92 = sort_col(mer92)
    mer94 = sort_col(mer94)
    mer96 = sort_col(mer96)
    mer97 = sort_col(mer97)
    # mer85.to_csv("v85_fasta_merge_xref_3955IDs_10188.csv", index=False)
    # mer92.to_csv("v92_fasta_merge_xref_3955IDs_10399.csv", index=False)
    # mer94.to_csv("v94_fasta_merge_xref_3955IDs_10616.csv", index=False)
    # mer96.to_csv("v96_fasta_merge_xref_3955IDs_10667.csv", index=False)
    # mer97.to_csv("v97_fasta_merge_xref_3955IDs_10568.csv", index=False)

    # ---------- shared stable ID keys in all ensembl releases ------
    m85, m92, m94, m96, m97 = shared_IDs(mer85, mer92, mer94, mer96, mer97, 'stableID_key')
    uniqueCount(m85, 'xref')
    uniqueCount(m92, 'xref')
    uniqueCount(m94, 'xref')
    uniqueCount(m96, 'xref')
    uniqueCount(m97, 'xref')

    # saving set for version ID difference in panel 2
    # m85.to_csv("v85_fasta_merge_xref_8865IDs_3889xref.csv", index=False)
    # m92.to_csv("v92_fasta_merge_xref_8865IDs_3889xref.csv", index=False)
    # m94.to_csv("v94_fasta_merge_xref_8865IDs_3889xref.csv", index=False)
    # m96.to_csv("v96_fasta_merge_xref_8865IDs_3889xref.csv", index=False)
    # m97.to_csv("v97_fasta_merge_xref_8865IDs_3889xref.csv", index=False)



