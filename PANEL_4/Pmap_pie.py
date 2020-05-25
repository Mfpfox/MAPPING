#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11.23.19
for Project MAP
@author: mfpfox

goal: pie chart of total CK and detected CK.

OUTPUT:
checking labeled column values:
  ever_C_ID  Count
0     False  15517
1      True   2915

checking labeled column values:
  ever_K_ID  Count
0     False  15809
1      True   2623

shape ukbCK labeled =  (18432, 30)
TOTAL protiens with only Cys detected:
(1340, 30)
TOTAL Cys in the proteins:  19812

TOTAL protiens with only Lysine detected:
(1048, 30)
TOTAL Lysine in the proteins:  34910

TOTAL protiens with both CK detected:
(1575, 30)
TOTAL Lysine in the proteins:  61917
TOTAL Cysteine in the proteins:  14944

TOTAL protiens with all C detected:
(2915, 30)
TOTAL C in the proteins:  34756

TOTAL protiens with all K detected:
(2623, 30)
TOTAL K in the proteins:  96827
"""
import os
import sys
import pandas as pd

sys.path.append("/Users/mariapalafox/Desktop/Toolbox")
from all_funx import *

os.chdir("/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/CYS_LYS_2012_to2018UKBpositions/REFERENCE_2018_mappedCK/SUBSETS/")


# files
# ever labeled Cys IDs
# ever labeled Lys IDs
# Both labeled

def main():
    # ukb file with be one annotated
    everC = "CYS2018_everlabeled_gene_valueCounts_2915.csv"
    everK = "LYS2018_everlabeled_gene_valueCounts_2623.csv"
    uni = "countsfromfastaCCDS_filter_counts_all_AA.csv"

    cys = read_pd(everC)
    lys = read_pd(everK)
    ukb = read_pd(uni)
    dropUnnamed(ukb)

    print("total of all lysines detected: ", lys.counts.sum())
    print("total of all cysteine detected: ", cys.counts.sum())

    Cls = cys.ID.tolist()
    Kls = lys.ID.tolist()

    ukbC = addcolumnconditional(Cls, ukb, 'ID', 'ever_C_ID')
    ukbCK = addcolumnconditional(Kls, ukbC, 'ID', 'ever_K_ID')
    print("shape ukbCK labeled = ", ukbCK.shape)
    print()

    # only C detected proteins
    Conly = ukbCK[(ukbCK['ever_C_ID'] == "True") & (ukbCK['ever_K_ID'] == "False")]
    print("TOTAL protiens with only Cys detected: ")
    print(Conly.shape)
    print("TOTAL Cys in the proteins: ", Conly.C.sum())
    print()

    # only K detected proteins
    Konly = ukbCK[(ukbCK['ever_C_ID'] == "False") & (ukbCK['ever_K_ID'] == "True")]
    print("TOTAL protiens with only Lysine detected: ")
    print(Konly.shape)
    print("TOTAL Lysine in the proteins: ", Konly.K.sum())
    print()

    # both C K detecte proteins
    bothck = ukbCK[(ukbCK['ever_C_ID'] == "True") & (ukbCK['ever_K_ID'] == "True")]
    print("TOTAL protiens with both CK detected: ")
    print(bothck.shape)
    print("TOTAL Lysine in the proteins: ", bothck.K.sum())
    print("TOTAL Cysteine in the proteins: ", bothck.C.sum())
    print()

    # all C IDs
    Call = ukbCK[ukbCK['ever_C_ID'] == "True"]
    print("TOTAL protiens with all C detected: ")
    print(Call.shape)
    print("TOTAL C in the proteins: ", Call.C.sum())
    print()

    Kall = ukbCK[ukbCK['ever_K_ID'] == "True"]
    print("TOTAL protiens with all K detected: ")
    print(Kall.shape)
    print("TOTAL K in the proteins: ", Kall.K.sum())

main()



