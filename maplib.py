# !/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Pmap functions

import os
import sys
import numpy as np
import pandas as pd
import csv
from Bio import SeqIO
from ast import literal_eval # for mismap_score func
import difflib # diff btw 2 strings
from statistics import mean
import jellyfish
import textdistance

### Basics ###

def create_coordinate_id(df, chrr, pos, ref, alt):
    # variables are colnames of df
    df.loc[:,'pos_coordinate'] = df[chrr].astype(str) + '_' + \
    df[pos].astype(str) + '_' + df[ref].astype(str) + \
    '_' + df[alt].astype(str)
    return df

def header_index(df):
    header = df.columns
    d = {header[i] : i for i in range(0, len(header))}
    print(d)

def stableID_key(df):
    df['stableID_key'] = df['gene_stable_id'].astype(str) + \
'_' + df['transcript_stable_id'].astype(str) + \
'_' + df['protein_stable_id'].astype(str)

def replace_col_value(df, colname, val_og, val_new):
    # replace col values 
    df[colname].replace(val_og, val_new, inplace=True)
    return df

def drop_duplicates(df, col):
    df.drop_duplicates(subset = col, keep = 'first', inplace = True) 
    return df

def change_col_type(df, col):
    # changing column types in pandas
    df[col] = pd.to_numeric(df[col])
    dffinal = df.sort_values(by=[col])
    return dffinal

def add_thresholds(df):
    # adding new threshold by old numbers from nature paper
    df.loc[df.reactivity <= 2, 'Threshold'] = 'High' 
    df.loc[(df.reactivity > 2) & (df.reactivity <= 5), 'Threshold'] = 'Medium' 
    df.loc[df.reactivity > 5, 'Threshold'] = 'Low'
    return df

def subtract_cols(df, col1, col2, newcol):
    # new column from subtracting other 2
    df = df.assign(newcol = df[col1] - df[col2])
    return df

def map_column_values(refdf, df, dickey, dicvalue):
    # adds gene symbol col
    # make sure refdf has dickey and dicvalue column
    # df must have dickey column 'ID'
    newcol = dicvalue # 'HGNC_name'
    refdic = dict(zip(refdf[dickey], refdf[dicvalue])) # makes dict
    df[newcol] = df[dickey] # new col 'gene' name with df key values
    df[dicvalue] = df[dicvalue].map(refdic)
    return df

def format_missense_triple(df, oaacol, naacol):
    #  A|A turns to Ala/Ala
    amino_dict = dict([('A', 'Ala'),('G', 'Gly'), ('I','Ile'), ('L','Leu'), ('P', 'Pro'), ('V','Val'), ('F','Phe'),('W', 'Trp'), ('Y', 'Tyr'), ('D','Asp'),('E','Glu'), ('R','Arg'),('H','His'), ('K','Lys'), ('S','Ser'), ('T', 'Thr'), ('C', 'Cys'), ('M', 'Met'), ('N', 'Asn'), ('Q','Gln')])
    df[oaacol].replace(amino_dict, inplace=True)
    df[naacol].replace(amino_dict, inplace=True)
    ccopy = df[naacol].copy()
    df['Amino_acids'] = df[oaacol].str.cat(ccopy, sep='/')
    return df

def setme(s):
    ls = list(s)
    finalset = set(ls)
    return finalset

def avgme(s):
    return mean(s)

def rangefinder(s):
    maxx = max(s)
    minn = min(s)
    diff = maxx - minn
    return diff

def sortbycols(df, colnamelist):
    # sorting df based on 2 columns, colname in list
    df.sort_values(by=colnamelist, inplace=True)
    df.reset_index(inplace=True, drop=True)

def sortbycol(df, colname):
    # sorting df based on 1 col
    df.sort_values(by=[colname], inplace=True)
    df.reset_index(inplace=True, drop=True)

def countChar(ls):
    # using count() to get C and K #
    Cys = 0
    Lys = 0
    for st in COUNTck:
        C = st.count('C')
        K = st.count('K')
        Cys += C
        Lys += K
    print(Cys)
    print(Lys)

def dropNotLabeled(df, Lever):
    # drops rows where ID does not match list IDs 'Lever'
    df['SharedUKB_IDs'] = np.where(df['xref'].isin(Lever), "True", "False")
    df2 = df[df['SharedUKB_IDs'] == "True"].copy()
    print("sorting ID column and resetting index: ")
    df2.sort_values(by=['xref'], inplace=True)
    df2.reset_index(drop=True, inplace=True)
    print("shape final cleaned df: ", df2.shape)
    return df2

########### more complex functions ##################

def feature_id_col(filename, outfile, subset):
    # from M19_12_27 markdown
    if subset == 'dect':
        # adding feature_ID column from parsing enst list col from dbNSFP with matched_index col
        # modeled after function from Pmap_parseID_correction.py
        with open(filename, newline='') as file:
            # read in file, save header
            csvReader = csv.reader(file)
            header = next(csvReader)
            # create and write to outfile
            os.system("touch %s" % (outfile))
            with open(outfile, 'w') as out:
                csvWriter = csv.writer(out)
                csvWriter.writerow(header)
            # loop over rows
            for row in csvReader:
                matchI = int(row[6])
                changeindex = [21,26,28,29,31,32,34,35,37,42,44,45,46,47,49,50,52,53,55,56,74,75]
                for coli in changeindex:
                    rowvalue = row[coli]
                    splitval = rowvalue.split(";")
                    if matchI < len(splitval):
                        newval = splitval[matchI]
                        row[coli] = newval
                with open(outfile, 'a') as out:
                        csvWriter = csv.writer(out)
                        csvWriter.writerow(row)
        print("done with : ", outfile)
    if subset == 'notdect':
        # adding feature_ID column from parsing enst list col from dbNSFP with matched_index col
        # modeled after function from Pmap_parseID_correction.py
        with open(filename, newline='') as file:
            # read in file, save header
            csvReader = csv.reader(file)
            header = next(csvReader)
            # create and write to outfile
            os.system("touch %s" % (outfile))
            with open(outfile, 'w') as out:
                csvWriter = csv.writer(out)
                csvWriter.writerow(header)
            # loop over rows
            for row in csvReader:
                matchI = int(row[5])
                changeindex = [6,19,21,22,24,25,27,28,30,35,37,38,39,40,42,43,45,46,48,49,67,68]
                for coli in changeindex:
                    rowvalue = row[coli]
                    splitval = rowvalue.split(";")
                    if matchI < len(splitval):
                        newval = splitval[matchI]
                        row[coli] = newval
                with open(outfile, 'a') as out:
                        csvWriter = csv.writer(out)
                        csvWriter.writerow(row)
        print("done with : ", outfile)

def avg_reactivity(df):
    # turning series into list
    listCreact = list(df.reactivity)
    # turning list into list of lists
    llcreact = [x.split() for x in listCreact]
    # removed last index of sublist which contained the / value
    for sublist in llcreact:
        del sublist[-1]
    # getting average of lists in list
    avgg = []
    for sublist in llcreact:
        test_list = [float(i) for i in sublist]
        s = sum(test_list)
        l = len(sublist)
        a = s/l
        avgg.append(a)
    # formating decimal place
    avgf = ["{0:.2f}".format(x) for x in avgg]
    print("shape of df: ", df.shape)
    print("shape of averaged column: ", len(avgf))
    # convert this into a column for CreactTRUE df
    avgR = pd.DataFrame(np.array(avgf).reshape(-1, 1))
    avgR.columns = ['avgReactivity']
    df2 = pd.concat([df, avgR], axis=1)
    df3 = df2[['pos_ID', 'avgReactivity']].copy()
    return df3

def mismaplines_dynamic(df):
# expanding this out to include dynamic levels - 1 2 3 
    diffline = []
    for index, row in df.iterrows():
        ukbid = row['ID']
        ls = row['frac_missed']
        #python_ls = literal_eval(ls) 
        python_ls = ls
        lenLS = len(set(python_ls))
        if lenLS == 5:
            diffline.append("5")
        if lenLS == 4:
            diffline.append("4")
        if lenLS == 3:
            diffline.append("3")
        if lenLS == 2:
            diffline.append("2")
        if lenLS == 1:
            diffline.append("1")
    df.loc[:,'dynamic_slope_scores'] = diffline
    print(df.shape)


### protein seq comparison ###

def sequenceDifference(dfref, dfalt):
    # assumes ID and proSequence column
    ref_dic = dict(zip(dfref.ID, dfref.proSequence))
    diffcol = []
    for index, row in dfalt.iterrows():
        pep = row['proSequence']
        ukb_id = row['ID']
        # retrieve uniprot fasta seq using reference dictionary
        mypep = ref_dic[ukb_id]
        str(mypep)
        if mypep == pep:
            diffcol.append('nodiff')
        if mypep != pep:
            output_list = [li for li in difflib.ndiff(mypep, pep) if li[0] != ' ']
            out = ', '.join(output_list)
            diffcol.append(out)
    # add identity score list as new column
    dfalt.loc[:, 'difference_12vs18'] = diffcol
    print(dfalt.shape)
    return dfalt

def get_pos_dictionary(pdcol):
    # creates list of Cys and Lys position keys
    list_of_dicts = []
    for row in pdcol.tolist():
        listofkey = []
        listofvalues = []
        for i in row:
            v = i[0]
            listofvalues.append(v)
            k = i[1:]
            listofkey.append(k)
        di = dict(zip(listofkey, listofvalues))
        intdict = {int(oldkey): val for oldkey, val in di.items()}
        sdi = dict(sorted(intdict.items(), reverse=False))
        list_of_dicts.append(sdi)
    return list_of_dicts

def identicalSequenceCheck(dfref, dfalt):
    # used for ukb to ukb 2018 2012 seq comparison dfref is 2018
    # assumes ID and proSequence column
    # compares the sequence in df and dict value, adds column based on comparison results
    # True is they are identical, False is they are different
    ref_dic = dict(zip(dfref.ID, dfref.proSequence))
    ref_len = dict(zip(dfref.ID, dfref.Length))
    newcol = []
    len2012 = []
    for index, row in dfalt.iterrows():
        pep = row['proSequence']
        ukb_id = row['ID']
        # retrieve uniprot fasta seq using reference dictionary
        mypep = ref_dic[ukb_id]
        len12 = ref_len[ukb_id]
        str(mypep)
        if mypep == pep:
            newcol.append('True')
            len2012.append(len12)
        if mypep != pep: 
            newcol.append('False')
            len2012.append(len12)
    # add identity score list as new column
    dfalt.loc[:, 'identical_seq'] = newcol
    dfalt.loc[:, 'length_2012_fasta'] = len2012
    dfalt.loc[:, 'len12_minus_len18'] = dfalt['length_2012_fasta'] - dfalt['Length']
    print(dfalt.shape)
    checkColumnValues(dfalt, "identical_seq")
    return dfalt

def identicalSequences(dfEnsp, ref_dic, newcol):
    # for comparison of ukb seq to ensp  seq
    # input df with proSequence col, a dict with ID and prosequence values
    # outputs a dataframe with new column based on match or not
    # outputs list of IDs that have different protein sequences
    # outputs counts of True False in new col
    res = []
    for index, row in dfEnsp.iterrows():
        pep = row['proSequence']
        idd = row['ENSP']
        # check pep to dict pep sequence
        mypep = ref_dic[idd]
        str(mypep)
        if mypep == pep:
            res.append('True')
        if mypep != pep:
            res.append('False')
    # add new column
    dfEnsp.loc[:, newcol] = res
    # pull out rows with different peptides same ENSP save as list
    dfFALSE = dfEnsp[dfEnsp[newcol] == 'False']
    dfFALSE_list = dfFALSE['ENSP'].tolist()
    # count # of same and different peptide for same ENSP
    checkColumnValues(dfEnsp, newcol)
    # print IDs that have different pep
    print(dfFALSE_list)
    return dfEnsp

def sequenceDistance(dfEnsp, ref_dic, newcolresult, hamming, hammingNorm, levenshtein, levenshteinNorm):
    res = []
    ham = []
    hamnorm = []
    lev = []
    levnorm = []
    serSeq = dfEnsp['proSequence'].copy()
    serID = dfEnsp['stableID_key'].copy()
    for inx, val in serSeq.items():
        pep = str(val)
        p = pep.strip()
        idd = str(serID[inx])
        # check pep to dict pep sequence
        mypep = ref_dic[idd]
        str(mypep)
        # identical
        if mypep == p:
            res.append('True')
            ham.append('identical')
            hamnorm.append('identical')
            lev.append('identical')
            levnorm.append('identical')
        # not identical to canonical
        if mypep != p: 
            res.append('False')
            # calculates hamming distance, penalizes positional differences, edit based distance
            ham.append(textdistance.hamming(mypep, p))
            # normalized hamming = # mismatched positions/ len of longer sequence
            hamnorm.append(textdistance.hamming.normalized_distance(mypep, p))
            # levenshtein score is edit based but not not penalized position, insertion at pos 1 is jsut 1 diff
            lev.append(textdistance.levenshtein(mypep,p))
            levnorm.append(textdistance.levenshtein.normalized_distance(mypep, p))
    dfEnsp.loc[:,newcolresult] = res
    dfEnsp.loc[:,hamming] = ham
    dfEnsp.loc[:,hammingNorm] = hamnorm
    dfEnsp.loc[:,levenshtein] = lev
    dfEnsp.loc[:,levenshteinNorm] = levnorm
    return dfEnsp
