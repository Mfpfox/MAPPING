# !/usr/bin/env python3
# -*- coding: utf-8 -*-

# filtering using updated count of uniprot IDs as 3,840 (number shared post dbNSFP mapping)

import os
import sys
import pandas as pd

def filter_by_col(infil, outname, fil, filcol):
    inn = read_pd(infil)
    fil_list = fil[filcol].tolist()
    filtered = addcolumnconditionalDrop(fil_list, inn, filcol, 'in_set')
    filtered.to_csv(outname, index=False)

def main():
    os.chdir("/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/CADDmapped/RESULT_CADDv14_pos_overlap_dbNSFPcoordinates/")
    # search files with pos_ID (ukbid_aapos)
    detfilter = "PULL_annotations/SEARCH_hg19_DETECTED_104475.csv"
    notfilter = "PULL_annotations/SEARCH_hg19_NOT_DETECTED_1222911.csv"
    detfil = pd.read_csv(detfilter, usecols=['pos_ID'])
    notfil = pd.read_csv(notfilter, usecols=['pos_ID'])
    detfil.drop_duplicates(keep='first', inplace=True)
    notfil.drop_duplicates(keep='first', inplace=True)
    print("detected post drop duplicate shape of df: ", detfil.shape)
    print("not detected post drop duplicate shape: ", notfil.shape)
    print(detfil.head(2))
    print()
    print(notfil.head(2))
     # saving non redundant on pos_ID search files
    os.chdir("/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/CADDmapped/BUG_FIX/NOT_3840_FILTERED")
    detfil.to_csv("FILTER_detected_pos_ID_14925.csv", index=False)
    notfil.to_csv("FILTER_notdetected_pos_ID_174702.csv", index=False)
    # K not detected
    filter_by_col("CADDspecific_K_notdetected_1015295x2rows.csv", "FILTERED_CADDspecific_K_notdetected_174702posID.csv", notfil, "pos_ID")
    # C not detected 
    filter_by_col("CADDspecific_C_notdetected_244506x2rows.csv", "FILTERED_CADDspecific_C_notdetected_174702posID.csv", notfil, "pos_ID")
    # K detected 
    filter_by_col("CADDspecific_K_detected_63014x2rows.csv", "FILTERED_CADDspecific_K_detected_14925.csv", detfil, "pos_ID")
    # C detected 
    filter_by_col("CADDspecific_C_detected_43050x2rows.csv", "FILTERED_CADDspecific_C_detected_14925posID.csv", detfil, "pos_ID")
main()



















