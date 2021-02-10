#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 10/2/19 Pmap_dbNSFP_CysLysfilter.py

# packages
import numpy as np
import pandas as pd
import csv


def main():
    # read in header file
    h = pd.read_csv("h.txt",sep='\t')
    hlist = h.columns
    print(len(hlist))
    print(hlist)
    # 376 IDs
    header = hlist

    infiles = ['accFiltered_chr1.csv',
    'accFiltered_chr2.csv',
    'accFiltered_chr3.csv',
    'accFiltered_chr4.csv',
    'accFiltered_chr5.csv',
    'accFiltered_chr6.csv',
    'accFiltered_chr7.csv',
    'accFiltered_chr8.csv',
    'accFiltered_chr9.csv',
    'accFiltered_chr10.csv',
    'accFiltered_chr11.csv',
    'accFiltered_chr12.csv',
    'accFiltered_chr13.csv',
    'accFiltered_chr14.csv',
    'accFiltered_chr15.csv',
    'accFiltered_chr16.csv',
    'accFiltered_chr17.csv',
    'accFiltered_chr18.csv',
    'accFiltered_chr19.csv',
    'accFiltered_chr20.csv',
    'accFiltered_chr21.csv',
    'accFiltered_chrX.csv',
    'accFiltered_chrY.csv']

    outtfiles = ['ckfiltered_chr1.csv',
    'ckfiltered_chr2.csv',
    'ckfiltered_chr3.csv',
    'ckfiltered_chr4.csv',
    'ckfiltered_chr5.csv',
    'ckfiltered_chr6.csv',
    'ckfiltered_chr7.csv',
    'ckfiltered_chr8.csv',
    'ckfiltered_chr9.csv',
    'ckfiltered_chr10.csv',
    'ckfiltered_chr11.csv',
    'ckfiltered_chr12.csv',
    'ckfiltered_chr13.csv',
    'ckfiltered_chr14.csv',
    'ckfiltered_chr15.csv',
    'ckfiltered_chr16.csv',
    'ckfiltered_chr17.csv',
    'ckfiltered_chr18.csv',
    'ckfiltered_chr19.csv',
    'ckfiltered_chr20.csv',
    'ckfiltered_chr21.csv',
    'ckfiltered_chrX.csv',
    'ckfiltered_chrY.csv']

    for ii, oo in zip(infiles,outtfiles):
        filename = ii
        outputfile = oo
        with open(filename, newline='') as file:
            csvReader = csv.reader(file)
            header = next(csvReader)
            # create output file and write to file immediately
            os.system("touch %s" % (outputfile))
            with open(outputfile, 'w') as out:
                csvWriter = csv.writer(out)
                # header
                csvWriter.writerow(header)
            for row in csvReader:
                aa = row[4]
                if aa == 'K' or aa == 'C':
                    with open(outputfile, 'a') as out:
                        csvWriter = csv.writer(out)
                        csvWriter.writerow(row)
            print("done with: ", ii)

main()
