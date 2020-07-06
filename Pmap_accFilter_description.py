#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# getting all files unique ID count, shape, and name in new description file

import os
import pandas as pd
import csv

def main():
    # create description file
    res = "accFiltering_description.csv"
    os.system("touch %s" % (res))
    with open(res, 'w') as out:
        csvWriter = csv.writer(out)

    # input files
    infiles = ['accFiltered_chr1.csv','accFiltered_chr2.csv','accFiltered_chr3.csv','accFiltered_chr4.csv','accFiltered_chr5.csv','accFiltered_chr6.csv','accFiltered_chr7.csv','accFiltered_chr8.csv','accFiltered_chr9.csv','accFiltered_chr10.csv','accFiltered_chr11.csv','accFiltered_chr12.csv','accFiltered_chr13.csv','accFiltered_chr14.csv','accFiltered_chr15.csv','accFiltered_chr16.csv','accFiltered_chr17.csv','accFiltered_chr18.csv','accFiltered_chr19.csv','accFiltered_chr20.csv','accFiltered_chr21.csv','accFiltered_chrX.csv','accFiltered_chrY.csv']
    # loop thru input
    for v in infiles:
        row = []
        df = pd.read_csv(v,low_memory=False)
        unqlen = len(df['matched_UKBID'].unique())
        row.append(v)
        row.append(df.shape)
        row.append(unqlen)
        with open(res, 'a') as out:
            csvWriter = csv.writer(out)
            csvWriter.writerow(row)
        print("done with ", v)
main()