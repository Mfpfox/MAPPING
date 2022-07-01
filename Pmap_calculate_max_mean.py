# !/usr/bin/env python3
# -*- coding: utf-8 -*-
# feb 5 2021 finding max and mean scores for all amino acid positions
import os
import sys
import numpy as np
import pandas as pd
import csv
from statistics import mean 

# only 1/4526 proteins from dbnsfp were missing AAs

# script updates
## 1. changed max() to np.nanmax() so max and mean DANN Fathmm counts match. previous version of results show more missing values for max than for mean, but max and mean missing are equal for DANN and FATHMM
## 2. on merge, prev. was "inner" now doing "left" to prevent 188 aa lost (P42166 = 188 missing amino acids compared to  ukb protein length)

def addcolumnconditional(mapList, df, dfcol, newcol):
    mendel = []
    for g in df[dfcol]:
        if g in mapList:
            mendel.append("True")
        else:
            mendel.append("False")
    df.loc[:, newcol] = mendel
    return df

def add_maxmean_col(df, dfanno, s1, s2, s3, s4, outputfile):
    max1 = s1 + ".max"
    mean1 = s1 + ".mean"
    max2 = s2 + ".max"
    mean2 = s2 + ".mean"
    max3 = s3 + ".max"
    mean3 = s3 + ".mean"
    max4 = s4 + ".max"
    mean4 = s4 + ".mean"
    df[s1] = df[s1].astype(float)
    df[s2] = df[s2].astype(float)
    df[s3] = df[s3].astype(float)
    df[s4] = df[s4].astype(float)
    dfg1 = df.groupby('posID',sort=False)[s1].apply(list)
    dfg1 = pd.DataFrame(dfg1)
    dfg2 = df.groupby('posID',sort=False)[s2].apply(list)
    dfg2 = pd.DataFrame(dfg2)
    dfg3 = df.groupby('posID',sort=False)[s3].apply(list)
    dfg3 = pd.DataFrame(dfg3)
    dfg4 = df.groupby('posID',sort=False)[s4].apply(list)
    dfg4 = pd.DataFrame(dfg4)
    #######CADD raw###############
    dfg1[mean1] = dfg1[s1].apply(lambda x: np.nanmean(x))
    dfg1[max1] = dfg1[s1].apply(lambda x: np.nanmax(x))
    dfg1.reset_index(inplace=True)
    dfg1.drop(s1, axis=1, inplace=True)
    dfanno = pd.merge(dfanno, dfg1, on=['posID'], how='left')
    ##########CADD phred############
    dfg2[mean2] = dfg2[s2].apply(lambda x: np.nanmean(x))
    dfg2[max2] = dfg2[s2].apply(lambda x: np.nanmax(x))
    dfg2.reset_index(inplace=True)
    dfg2.drop(s2, axis=1, inplace=True)
    dfanno = pd.merge(dfanno, dfg2, on=['posID'], how='left')
    ##########DANN############
    dfg3[mean3] = dfg3[s3].apply(lambda x: np.nanmean(x))
    dfg3[max3] = dfg3[s3].apply(lambda x: np.nanmax(x))
    dfg3.reset_index(inplace=True)
    dfg3.drop(s3, axis=1, inplace=True)
    dfanno = pd.merge(dfanno, dfg3, on=['posID'], how='left')
    ###########fathmmMKL###########
    dfg4[mean4] = dfg4[s4].apply(lambda x: np.nanmean(x))
    dfg4[max4] = dfg4[s4].apply(lambda x: np.nanmax(x))
    dfg4.reset_index(inplace=True)
    dfg4.drop(s4, axis=1, inplace=True)
    dfanno = pd.merge(dfanno, dfg4, on=['posID'], how='left')
    ######################
    dfanno = dfanno.round(3)
    print(dfanno.head(5))
    print("saving merge")
    dfanno.to_csv(outputfile, index=False)


file = sys.argv[1] # selectCols_appV1_chr${i}.csv
outputfile = sys.argv[2] # codonScores_allLetters_chr${i}.csv
df = pd.read_csv(file, low_memory=False)
dfanno = df[["posID", "aaref", "matched_UKBID","matched_aapos"]].copy()
dfanno.drop_duplicates(inplace=True)
df = df[["posID","CADD38raw","CADD38phred","DANN","fathmmMKL"]].copy()
# convert "." to na value
df.replace(".", np.NaN, inplace=True)
# accounting for aa pos with missing scores
missing = df[df.isnull().any(axis=1)]
pos = list(set(missing['posID']))
dfanno = addcolumnconditional(pos, dfanno, 'posID', 'missingScores')
add_maxmean_col(df, dfanno, "CADD38raw","CADD38phred","DANN","fathmmMKL",outputfile)
