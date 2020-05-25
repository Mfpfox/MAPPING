#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Pmap_ultimate_UKB.py

import pandas as pd
import os 
import sys
os.chdir("/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/TSV_UNIPROT_xref/")
    sys.path.append("/u/home/m/mfpalafo/project-arboleda/Toolbox")
from all_funx import *

def main():
    # abundance file
    myukb = pd.read_csv("CYSLYS_failure2map/countsfromfastaCCDS_filter_normalize_all_AA.csv")
    udf = myukb[['ID', 'Length', 'proSequence','C', 'K']].copy()

    # labeled accessions
    hitset = "CYS_LYS_2012_to2018UKBpositions/all3991_everLabeled_entries_simple.csv"
    hset = pd.read_csv(hitset)
    hls = hset['Entry'].tolist()

    dflabeled = addcolumnconditional(hls, udf, 'ID','chemoproteomicSet')
    dflabeled.to_csv("UKBCCDS_abundanceCK_detected_18432.csv",index=False)

main()
