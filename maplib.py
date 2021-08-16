# !/usr/bin/env python3
# -*- coding: utf-8 -*-
## author: maria f. palafox

import os
import sys
import numpy as np
import pandas as pd
import csv
from Bio import SeqIO
from ast import literal_eval 
import difflib 
from statistics import mean
import jellyfish
import textdistance
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# [1]
def UKBposID(df, IDcol, AAposCol):
    # make uniprotID + AApos key column
    df['pos_ID'] = df[IDcol].astype(str) + \
        '_' + df[AAposCol].astype(str)
    return(df)

# [2]
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

# [3]
def add_thresholds(df):
    # adding reactivity threshold labels, assumes colname 'reactivity'
    # High no lower bound ... r <= 2
    # Med 2 < r <= 5
    # Low r > 5
    # used for 2012 experiments
    # for 2019 experiments first filtered out rows w/ r < 0.5
    df.loc[df.reactivity <= 2, 'Threshold'] = 'High' 
    df.loc[(df.reactivity > 2) & (df.reactivity <= 5), 'Threshold'] = 'Medium' 
    df.loc[df.reactivity > 5, 'Threshold'] = 'Low'
    return df

# [4]
def checkColumnValues(df, col):
    print(df[col].value_counts().reset_index().rename(columns={'index':col, col:'Count'}))


# [5]
def uniprotkbFasta2csv(filename):
    identifiers = [seq_record.id for seq_record in SeqIO.parse(filename, "fasta")]
    lengths = [len(seq_record.seq) for seq_record in SeqIO.parse(filename, "fasta")]
    proSeq = [seq_record.seq for seq_record in SeqIO.parse(filename, "fasta")]
    # splitting identifiers into 2 seperate list
    splitAcc = []
    splitEntry = []
    proseq = []
    for id in identifiers:
        splitID = id.split('|')
        acc = splitID[1]
        splitAcc.append(acc)
        entryName = splitID[2]
        splitEntry.append(entryName)
    for ps in proSeq:
        s = str(ps)
        proseq.append(s)
    s1 = pd.Series(splitAcc, name='ID')
    s2 = pd.Series(splitEntry, name= 'entryName')
    s3 = pd.Series(lengths, name='Length')
    s4 = pd.Series(proseq, name='proSequence')
    series = [s1,s2,s3,s4]
    df = pd.concat(series, axis=1)
    return(df)

# [6]
def make_aapos_key(df, aa):
    # assumes colnames of 'proSequence' and 'ID' to make
    # residue-level key id
    annotation = []
    for index, row in df.iterrows():
        raw_seq = str(row['proSequence'])
        entry = str(row["ID"])
        for i, j in enumerate(raw_seq):
            if j == aa:
                pos = str(i+1) #1 index correction
                annotation.append(entry + '_' + aa + pos)
    final = pd.DataFrame(annotation)
    final.columns = ['pos_ID']
    print("Number of", aa, " amino acids in fasta", final.shape)
    return(final)

# [7]
def addcolumnconditional(mapList, df, dfcol, newcol):
    ls = []
    for g in df[dfcol]:
        if g in mapList:
            ls.append("True")
        else:
            ls.append("False")
    df.loc[:, newcol] = ls
    print("results of addColumnConditional: ")
    checkColumnValues(df, newcol)
    print()
    return df

# [8]
def addcolumnconditionalDrop(mapList, df, dfcol, newcol):
    ls = []
    for g in df[dfcol]:
        if g in mapList:
            ls.append("True")
        else:
            ls.append("False")
    df.loc[:, newcol] = ls
    print("pre drop df shape: ", df.shape)
    df = df[df[newcol] == "True"].copy()
    df.drop(newcol, axis=1, inplace=True)
    df.reset_index(inplace=True, drop=True)
    print("post drop df shape: ", df.shape)
    return df


# [10]
def makeCHRfiles(df, savename):
    # accepts df w/ uniprot IDs and col for 'chr'
    # savename is extenstion for chr chunk files, e.g. chr1'_savename.csv'
    chrlist = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7',
    'chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15',
    'chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX',
    'chrY']
    for i in chrlist:
        ch = str(i)
        dff = df[df['chr'] == ch].copy()
        dff.drop(['chr'], inplace=True, axis=1) # drop col with chr
        print(i, dff.shape)
        dff.to_csv(ch+savename,index=False)

# [11]
def split_ID(df1, col):
    new1 = df1[col].str.split("_", n=1, expand = True)
    df1["ID"] = new1[0]
    return df1

# [12]
def uniqueCount(df, colname):
    print("total values in column: ", len(df[colname]))
    print("total unique values: ", len(df[colname].unique()))
    print()

# [13]
def ENSPfasta2DF(filename, release):
    # Pmap_ensembl_fasta.py - parse sequence fasta file
    identifiers = [seq_record.id for seq_record in
                   SeqIO.parse(filename, "fasta")]
    descr = [seq_record.description for seq_record in
             SeqIO.parse(filename, "fasta")]
    lengths = [len(seq_record.seq) for seq_record in
               SeqIO.parse(filename, "fasta")]
    proSeq = [seq_record.seq for seq_record in SeqIO.parse(filename, "fasta")]
    ensp = []
    enspv = []
    proseq = []
    ensg = []
    ensgv = []
    enst = []
    enstv = []
    assembly = []
    chrom = []
    DNAstart = []
    DNAstop = []
    if release == "v85":
        for id in identifiers:
            splitID = id.split('.')
            stable = splitID[0]
            ensp.append(stable)
            enspv.append(id)
        for ps in proSeq:
            s = str(ps)
            proseq.append(s)
        for row in descr:
            splitrow = row.split(" ")
            for i, val in enumerate(splitrow):
                if 'chromosome' in val:
                    locationSplit = val.split(":")
                    leng = len(locationSplit)
                    if leng == 6:
                        grch = locationSplit[1]
                        assembly.append(grch)
                        chrr = locationSplit[2]
                        chrom.append(chrr)
                        start = locationSplit[3]
                        DNAstart.append(start)
                        stopl = locationSplit[4]
                        DNAstop.append(stopl)
                    else:
                        assembly.append(None)
                        chrom.append(None)
                        DNAstart.append(None)
                        DNAstop.append(None)
                if 'ENSG' in val:
                    gene = val.split(":")[1]
                    stablegene = gene.split(".")[0]
                    ensgv.append(gene)
                    ensg.append(stablegene)
                if 'ENST' in val:
                    tx = val.split(":")[1]
                    stabletx = tx.split(".")[0]
                    enstv.append(tx)
                    enst.append(stabletx)
        nA = pd.Series(assembly, name='Assembly')
        nC = pd.Series(chrom, name='chromosome')
        nSta = pd.Series(DNAstart, name='start')
        nSto = pd.Series(DNAstop, name='stop')
        s2 = pd.Series(enspv, name='ENSPv')
        s1 = pd.Series(ensp, name='ENSP')
        s7 = pd.Series(lengths, name='ENSPlength', dtype=int)
        s8 = pd.Series(proseq, name='proSequence')
        s4 = pd.Series(enstv, name='ENSTv')
        s3 = pd.Series(enst, name='ENST')
        s6 = pd.Series(ensgv, name='ENSGv')
        s5 = pd.Series(ensg, name='ENSG')
        series = [s1, s2, s7, s8, s3, s4, s5, s6, nA, nC, nSta, nSto]
        df = pd.concat(series, axis=1)
        df.dropna(inplace=True)
        df.reset_index(drop=True, inplace=True)
        return (df)
    # release is not v85
    else:
        HGNCsymbol = []
        HGNCdescription = []
        for id in identifiers:
            splitID = id.split('.')
            stable = splitID[0]
            ensp.append(stable)
            enspv.append(id)
        for ps in proSeq:
            s = str(ps)
            proseq.append(s)
        for row in descr:
            splitrow = row.split(" ")
            for i, val in enumerate(splitrow):
                if 'chromosome' in val:
                    locationSplit = val.split(":")
                    leng = len(locationSplit)
                    if leng == 6:
                        grch = locationSplit[1]
                        assembly.append(grch)
                        chrr = locationSplit[2]
                        chrom.append(chrr)
                        start = locationSplit[3]
                        DNAstart.append(start)
                        stopl = locationSplit[4]
                        DNAstop.append(stopl)
                    else:
                        assembly.append(None)
                        chrom.append(None)
                        DNAstart.append(None)
                        DNAstop.append(None)
                if 'ENSG' in val:
                    gene = val.split(":")[1]
                    stablegene = gene.split(".")[0]
                    ensgv.append(gene)
                    ensg.append(stablegene)
                if 'ENST' in val:
                    tx = val.split(":")[1]
                    stabletx = tx.split(".")[0]
                    enstv.append(tx)
                    enst.append(stabletx)
                if 'gene_symbol' in val:
                    sym = val.split(":")[1]
                    HGNCsymbol.append(sym)
                if 'description' in val:
                    # get index of val in row
                    pos = i
                    restOfRow = splitrow[pos:]
                    restOfRow = ' '.join(restOfRow)
                    HGNCdescription.append(restOfRow)
        nS = pd.Series(HGNCsymbol, name='HGNCsymbol')
        nD = pd.Series(HGNCdescription, name='HGNCdescription')
        nA = pd.Series(assembly, name='Assembly')
        nC = pd.Series(chrom, name='chromosome')
        nSta = pd.Series(DNAstart, name='start')
        nSto = pd.Series(DNAstop, name='stop')
        s2 = pd.Series(enspv, name='ENSPv')
        s1 = pd.Series(ensp, name='ENSP')
        s7 = pd.Series(lengths, name='ENSPlength', dtype=int)
        s8 = pd.Series(proseq, name='proSequence')
        s4 = pd.Series(enstv, name='ENSTv')
        s3 = pd.Series(enst, name='ENST')
        s6 = pd.Series(ensgv, name='ENSGv')
        s5 = pd.Series(ensg, name='ENSG')
        series = [s1, s2, s7, s8, s3, s4, s5, s6, nA, nC, nSta, nSto, nS, nD]
        df = pd.concat(series, axis=1)
        df.dropna(inplace=True)
        df.reset_index(drop=True, inplace=True)
        return (df)

# [14]
def create_coordinate_id(df, chrr, pos, ref, alt):
    # variables are colnames of df
    df.loc[:,'pos_coordinate'] = df[chrr].astype(str) + '_' + \
    df[pos].astype(str) + '_' + df[ref].astype(str) + \
    '_' + df[alt].astype(str)
    return df

# [15]
def header_index(df):
    header = df.columns
    d = {header[i] : i for i in range(0, len(header))}
    print(d)

# [16]
def stableID_key(df):
    df['stableID_key'] = df['gene_stable_id'].astype(str) + \
'_' + df['transcript_stable_id'].astype(str) + \
'_' + df['protein_stable_id'].astype(str)

# [17]
def replace_col_value(df, colname, val_og, val_new):
    # replace col values 
    df[colname].replace(val_og, val_new, inplace=True)
    return df

# [18]
def drop_duplicates(df, col):
    df.drop_duplicates(subset = col, keep = 'first', inplace = True) 
    return df

# [19]
def change_col_type(df, col):
    # changing column types in pandas
    df[col] = pd.to_numeric(df[col])
    dffinal = df.sort_values(by=[col])
    return dffinal

# [20]
def subtract_cols(df, col1, col2, newcol):
    # new column from subtracting other 2
    df = df.assign(newcol = df[col1] - df[col2])
    return df

# [21]
def map_column_values(refdf, df, dickey, dicvalue):
    # adds gene symbol col
    # make sure refdf has dickey and dicvalue column
    # df must have dickey column 'ID'
    newcol = dicvalue # 'HGNC_name'
    refdic = dict(zip(refdf[dickey], refdf[dicvalue])) # makes dict
    df[newcol] = df[dickey] # new col 'gene' name with df key values
    df[dicvalue] = df[dicvalue].map(refdic)
    return df

# [22]
def format_missense_triple(df, oaacol, naacol):
    #  A|A turns to Ala/Ala
    amino_dict = dict([('A', 'Ala'),('G', 'Gly'), ('I','Ile'), ('L','Leu'), ('P', 'Pro'), ('V','Val'), ('F','Phe'),('W', 'Trp'), ('Y', 'Tyr'), ('D','Asp'),('E','Glu'), ('R','Arg'),('H','His'), ('K','Lys'), ('S','Ser'), ('T', 'Thr'), ('C', 'Cys'), ('M', 'Met'), ('N', 'Asn'), ('Q','Gln')])
    df[oaacol].replace(amino_dict, inplace=True)
    df[naacol].replace(amino_dict, inplace=True)
    ccopy = df[naacol].copy()
    df['Amino_acids'] = df[oaacol].str.cat(ccopy, sep='/')
    return df

# [23]
def setme(s):
    ls = list(s)
    finalset = set(ls)
    return finalset

# [24]
def avgme(s):
    return mean(s)

# [25]
def rangefinder(s):
    maxx = max(s)
    minn = min(s)
    diff = maxx - minn
    return diff

# [26]
def sortbycols(df, colnamelist):
    # sorting df based on 2 columns, colname in list
    df.sort_values(by=colnamelist, inplace=True)
    df.reset_index(inplace=True, drop=True)

# [27]
def sortbycol(df, colname):
    # sorting df based on 1 col
    df.sort_values(by=[colname], inplace=True)
    df.reset_index(inplace=True, drop=True)

# 
def addMaxMean(df, scorecol):
    maxscore = scorecol + "_max"
    meanscore = scorecol + "_mean"
    # grouping on pos_ID
    df[scorecol] = df[scorecol].astype(float)
    df = df.groupby('pos_ID',sort=False)[scorecol].apply(list)
    # convert back to pd dataframe
    df = pd.DataFrame(df)
    df[meanscore] = df[scorecol].apply(lambda x: mean(x))
    df[maxscore] = df[scorecol].apply(lambda x: max(x))
    df.reset_index(inplace=True)
    df.drop(scorecol, axis=1, inplace=True)
    return df

# [28]
def dropNotLabeled(df, Lever):
    # drops rows where ID does not match list IDs 'Lever'
    df['SharedUKB_IDs'] = np.where(df['xref'].isin(Lever), "True", "False")
    df2 = df[df['SharedUKB_IDs'] == "True"].copy()
    print("sorting ID column and resetting index: ")
    df2.sort_values(by=['xref'], inplace=True)
    df2.reset_index(drop=True, inplace=True)
    print("shape final cleaned df: ", df2.shape)
    return df2

# [29]
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

# [30]
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

# [31]
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
