#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 10/2/19 Pmap_dbNSFP_parselabeledID.py

# packages
import numpy as np
import pandas as pd
import csv


def main():

    # dbNSFP list
    infiles = ['dbNSFP4.0a_variant.chr1',
    'dbNSFP4.0a_variant.chr2',
    'dbNSFP4.0a_variant.chr3',
    'dbNSFP4.0a_variant.chr4',
    'dbNSFP4.0a_variant.chr5',
    'dbNSFP4.0a_variant.chr6',
    'dbNSFP4.0a_variant.chr7',
    'dbNSFP4.0a_variant.chr8',
    'dbNSFP4.0a_variant.chr9',
    'dbNSFP4.0a_variant.chr10',
    'dbNSFP4.0a_variant.chr11',
    'dbNSFP4.0a_variant.chr12',
    'dbNSFP4.0a_variant.chr13',
    'dbNSFP4.0a_variant.chr14',
    'dbNSFP4.0a_variant.chr15',
    'dbNSFP4.0a_variant.chr16',
    'dbNSFP4.0a_variant.chr17',
    'dbNSFP4.0a_variant.chr18',
    'dbNSFP4.0a_variant.chr19',
    'dbNSFP4.0a_variant.chr20',
    'dbNSFP4.0a_variant.chr21',
    'dbNSFP4.0a_variant.chrX',
    'dbNSFP4.0a_variant.chrY']
    # output files
    outputfile = ['accFiltered_chr1.csv',
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

    # chr specific accessions that have been labeled and in UP5640
    chrspecific = ['chr1_only_accessions_from3901_proteome5640bed.csv',
    'chr2_only_accessions_from3901_proteome5640bed.csv',
    'chr3_only_accessions_from3901_proteome5640bed.csv',
    'chr4_only_accessions_from3901_proteome5640bed.csv',
    'chr5_only_accessions_from3901_proteome5640bed.csv',
    'chr6_only_accessions_from3901_proteome5640bed.csv',
    'chr7_only_accessions_from3901_proteome5640bed.csv',
    'chr8_only_accessions_from3901_proteome5640bed.csv',
    'chr9_only_accessions_from3901_proteome5640bed.csv',
    'chr10_only_accessions_from3901_proteome5640bed.csv',
    'chr11_only_accessions_from3901_proteome5640bed.csv',
    'chr12_only_accessions_from3901_proteome5640bed.csv',
    'chr13_only_accessions_from3901_proteome5640bed.csv',
    'chr14_only_accessions_from3901_proteome5640bed.csv',
    'chr15_only_accessions_from3901_proteome5640bed.csv',
    'chr16_only_accessions_from3901_proteome5640bed.csv',
    'chr17_only_accessions_from3901_proteome5640bed.csv',
    'chr18_only_accessions_from3901_proteome5640bed.csv',
    'chr19_only_accessions_from3901_proteome5640bed.csv',
    'chr20_only_accessions_from3901_proteome5640bed.csv',
    'chr21_only_accessions_from3901_proteome5640bed.csv',
    'chrX_only_accessions_from3901_proteome5640bed.csv',
    'chrY_only_accessions_from3901_proteome5640bed.csv']


    for ii, oo, a in zip(infiles, outputfile, chrspecific):
        filename = ii
        outfile = oo
        accession = a

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

            # create and write to outfile
            os.system("touch %s" % (outfile))
            with open(outfile, 'w') as out:
                csvWriter = csv.writer(out)
                csvWriter.writerow(header)

            # loop over rows
            for row in csvReader:
                # clean ID column
                uniprotIDs = row[16]
                nodot = uniprotIDs.replace("."," ")
                cleanrow = nodot.replace(";"," ")

                # create list then set
                splitrow = cleanrow.split()
                dbcolset = set(splitrow)
                # find matching accession between two sets
                addme = list(labelset & dbcolset)

                # if return is empty go to next row
                # else add row to outfile
                if len(addme) != 0:
                    # turn list into string
                    addstring = "".join(addme)
                    # add to row
                    row.append(addstring)
                    with open(outfile, 'a') as out:
                            csvWriter = csv.writer(out)
                            csvWriter.writerow(row)

        print("done with : ", outfile)

main()
