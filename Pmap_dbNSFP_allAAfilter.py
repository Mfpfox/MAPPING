import pandas as pd
import csv
import sys
import os
# Feb. 5, 2021 for all letters on CpDAA app v1.0

# 1. makes header out of selected column names, posID in index0
# 2. writes header to outfile
# 3. loops over rows, if aaref and aaalt column dont have . or X, 
    # 4. make posID for uniprot id and grch38 and 37 keyIDs
    # 5. apend row positions to newrow[] in order of header
    # 6. write newrow to outfile

## output files will next be checked for missense per codon, compared to uniprot protein lengths, and scored by max and mean missense for each posID

filename = sys.argv[1] # 'allAA_IDFiltered_chr1.csv'
outputfile = sys.argv[2] # 'selectCols_appV1_chr1.csv'

notref = ['X', '.'] # to match uniprot length, remove X aa from ref/alt
notalt = ['X', '.']

with open(filename, newline='') as file:
    csvReader = csv.reader(file)
    oldheader = next(csvReader) # dont use
    header = ["posID" , "aaref", "aaalt", 
        "matched_UKBID" , "matched_aapos", "keyID38", "keyID19", 
        "chr38" , "pos38" , "ref" , "alt" , "refcodon" , "codonpos" , 
        "CADD38raw" , "CADD38phred" , "DANN" , "fathmmMKL"]
    # create output and write immediately
    os.system("touch %s" % (outputfile))
    with open(outputfile, 'w') as out:
        csvWriter = csv.writer(out)
        csvWriter.writerow(header) # adds custom header
    for row in csvReader:
        newrow = []
        aaref = row[4] ## removed stop gain and other strand '.' rows
        aalt = row[5]
        if (aaref not in notref) and (aalt not in notalt):
            posID = row[376] + '_' + aaref + row[377] #UKBID_AApos
            pos38 = '{:0>9}'.format(row[1])
            pos19 = '{:0>9}'.format(row[8])
            keyID38 = row[0] + '_' + pos38 + '_' + row[2] + '_' + row[3]
            keyID19 = row[7] + '_' + pos19 + '_' + row[2] + '_' + row[3]
            newrow.append(posID)
            newrow.append(aaref)
            newrow.append(aalt)
            newrow.append(row[376])
            newrow.append(row[377])
            newrow.append(keyID38)
            newrow.append(keyID19)
            newrow.append(row[0])
            newrow.append(row[1])
            newrow.append(row[2])
            newrow.append(row[3])
            newrow.append(row[29])
            newrow.append(row[30])
            newrow.append(row[101])
            newrow.append(row[103])
            newrow.append(row[104])
            newrow.append(row[106])
            with open(outputfile, 'a') as out:
                csvWriter = csv.writer(out)
                csvWriter.writerow(newrow)
    print("done filtering ", filename)



