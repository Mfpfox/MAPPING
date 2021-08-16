#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import csv
import sys
import os

def grabIDaapos(lsID, lspos, lsdet):
    result = []
    for ii, val in enumerate(lsID):
        if val != '.':
            if '-' not in val:
                if val in lsdet:
                    aa = lspos[ii]
                    # ID, aapos, index is ls returned
                    result.append(val)
                    result.append(aa)
                    result.append(ii)
                    return result
    return('0')

def main():
    filename = sys.argv[1]
    outfile = sys.argv[2]
    accession = sys.argv[3]
    # read chr specific and labeled before uniprot accession
    ukbid = pd.read_csv(accession, header=None)
    labels = ukbid[0].tolist()
    labelset = set(labels)
    with open(filename, newline='') as file:
        # read in file, save header
        csvReader = csv.reader(file, delimiter='\t')
        header = next(csvReader)
        # add new col
        header.append('matched_UKBID')
        header.append('matched_aapos')
        header.append('matched_index')
        # create and write to outfile
        os.system("touch %s" % (outfile))
        with open(outfile, 'w') as out:
            csvWriter = csv.writer(out)
            csvWriter.writerow(header)
        # loop over rows
        for row in csvReader:
            # create list
            uniprotIDs = row[16]
            aapos = row[11]
            lsID = uniprotIDs.split(";")
            lsAA = aapos.split(";")
            addme = grabIDaapos(lsID, lsAA, labelset)
            # if return is 0 go to next row
            # else add row to outfile
            if len(addme) == 3:
                # turn list into string
                row.append(addme[0])  # ukb id
                row.append(addme[1])  # aapos
                row.append(addme[2])  # index
                with open(outfile, 'a') as out:
                        csvWriter = csv.writer(out)
                        csvWriter.writerow(row)
    print("done with : ", outfile)
main()
