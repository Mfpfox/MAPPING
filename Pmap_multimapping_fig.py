#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# multimapping 'fold change of stable IDs' figure

import pandas as pd
import os
import sys
sys.path.append("/Users/mariapalafox/Desktop/Toolbox/")
from all_funx import *

class mydict(dict):
    def __init__(self):
        self = dict()
    def add(self, key, value):
        self[key] = value

def countUniqueValues(df, colname):
    dict1 = mydict()
    dict2 = mydict()
    for i in colname:
        colum = df[i].tolist()
        setcol = set(colum)
        dict1.add(i, len(setcol))
        dict2.add(i, setcol)
    print(dict1)
    return(dict2)

def dictValueOverlap(d1, d2, d3, d4, d5):
    # loop over each level save value
    newdic = mydict()
    newdiclen = mydict()
    levs = ['gene_stable_id', 'transcript_stable_id', 'protein_stable_id']
    for i in levs:
        v1 = d1[i]
        v2 = d2[i]
        v3 = d3[i]
        v4 = d4[i]
        v5 = d5[i]
        overla = v1 & v2 & v3 & v4 & v5
        newdic.add(i, overla)
        newdiclen.add(i, len(overla))
    print(newdiclen) # {'gene_stable_id': 4170, 'transcript_stable_id': 8939, 'protein_stable_id': 8939}
    return newdic 

def main():
    # stable ID columns
    collist = ["gene_stable_id", "transcript_stable_id", "protein_stable_id", "xref"]
    df85 = pd.read_csv("sharedv85_10272_uniqueENSP.csv")
    df92 = pd.read_csv("sharedv92_10479_uniqueENSP.csv")
    df94 = pd.read_csv("sharedv94_10699_uniqueENSP.csv")
    df96 = pd.read_csv("sharedv96_10750_uniqueENSP.csv")
    df97 = pd.read_csv("sharedv97_10650_uniqueENSP.csv")
    dic85 = countUniqueValues(df85, collist)
    dic92 = countUniqueValues(df92, collist)
    dic94 = countUniqueValues(df94, collist)
    dic96 = countUniqueValues(df96, collist)
    dic97 = countUniqueValues(df97, collist)
    dfoverlap = dictValueOverlap(dic85, dic92, dic94, dic96, dic97)
    print(dfoverlap)
main()