# !/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Pmap functions - mapping paper version of all_funx.py

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

"""
def str_identity_check(str1, str2):
    # position of difference between two strings
    if len(str1) >= len(str2)
    output_list = [li for li in difflib.ndiff(str1, str2) if li[0] != ' ']
#print(posoutput_list)
output_list2 = [li for li in difflib.ndiff(altpep, refpep) if li[0] != ' ']
posoutput_list2 = [i+1 for i in range(len(altpep)) if altpep[i] != refpep[i]]
print()
print(output_list2)
print(posoutput_list2)



# new col for R label True is HIGH value and False if LOWMEDIUM
# new col for panreactive True if Target False if panreactive
def relabel_reactiveC(row):
    if row['Cys_react_threshold'] == "High":
        return 'True'
    if row['Cys_react_threshold'] == "Low":
        return 'False'
    if row['Cys_react_threshold'] == "Medium":
        return 'False'
    return 'Other'
def relabel_reactiveK(row):    
    if row['Lys_react_threshold'] == "High":
        return 'True'
    if row['Lys_react_threshold'] == "Low":
        return 'False'
    if row['Lys_react_threshold'] == "Medium":
        return 'False'
    return 'Other'
def relabel_targetC(row):
    if row['Cys_target_label'] == "Target":
        return 'True'
    if row['Cys_target_label'] == "panReactive":
        return 'False'
    return 'Other'
def relabel_targetK(row):    
    if row['Lys_target_label'] == "Target":
        return 'True'
    if row['Lys_target_label'] == "panReactive":
        return 'False'
    return 'Other'
# calling relabel functions
react_cys['HIGH_subgroup'] = react_cys.apply(lambda row: relabel_reactiveC(row), axis=1)
react_lys['HIGH_subgroup'] = react_lys.apply(lambda row: relabel_reactiveK(row), axis=1)
pan_cys['TARGET_subgroup'] = pan_cys.apply(lambda row: relabel_targetC(row), axis=1)
pan_lys['TARGET_subgroup'] = pan_lys.apply(lambda row: relabel_targetK(row), axis=1)
"""

########### Basics #################

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


#-----adding genes and uniprot ID ----------------

def map_column_values(refdf, df, dickey, dicvalue):
    # make sure refdf has dickey and dicvalue column
    # df must have dickey column 'ID'
    newcol = dicvalue # 'HGNC_name'
    refdic = dict(zip(refdf[dickey], refdf[dicvalue])) # makes dict
    df[newcol] = df[dickey] # new col 'gene' name with df key values
    df[dicvalue] = df[dicvalue].map(refdic)
    return df

def add_genename(df, colkey):
    genefile = "/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/TSV_UNIPROT_xref/IDs/GENE_ID_KEY_NOEXCEL_OPENING.csv"
    genekey=pd.read_csv(genefile, usecols=[0,3]) # 0=ID 3=HGNC_name
    print(genekey.shape)
    df2 = map_column_values(genekey, df, 'sUKBID', 'HGNC_name')
    return df2
#---------------------


def format_missense_triple(df, oaacol, naacol):
    #  A|A turns to Ala/Ala
    amino_dict = dict([('A', 'Ala'),('G', 'Gly'), ('I','Ile'), ('L','Leu'), ('P', 'Pro'), ('V','Val'), ('F','Phe'),('W', 'Trp'), ('Y', 'Tyr'), ('D','Asp'),('E','Glu'), ('R','Arg'),('H','His'), ('K','Lys'), ('S','Ser'), ('T', 'Thr'), ('C', 'Cys'), ('M', 'Met'), ('N', 'Asn'), ('Q','Gln')])
    df[oaacol].replace(amino_dict, inplace=True)
    df[naacol].replace(amino_dict, inplace=True)
    ccopy = df[naacol].copy()
    df['Amino_acids'] = df[oaacol].str.cat(ccopy, sep='/')
    return df


### functions that work on list #### used with .apply()

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

######

def sortbycols(df, colnamelist):
    # sorting df based on 2 columns, colname in list
    df.sort_values(by=colnamelist, inplace=True)
    df.reset_index(inplace=True, drop=True)
    #sortbycol(fasta85, colnamelist)

def sortbycol(df, colname):
    # sorting df based on 1 col
    df.sort_values(by=[colname], inplace=True)
    df.reset_index(inplace=True, drop=True)
    #dflist = [TF85,TF92,TF94,TF96,TF97, key85, key92, key94, key96, key97]
    #for i in dflist:
        #sortbycol(i, 'ENSP')

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


def group_scores(df, colname):
    # hamming_normalized_dist
    # levenshtein_normalized_dist
    lsdf = df.groupby('stableID_key')[colname].apply(list)
    lsdf = pd.DataFrame(lsdf)
        # concat all frames together
        #set_columns = pd.concat([lsdf, avgdf], axis=1)
        # problem with duplicated colnames from concat files
    lsdf['avg_distscr'] = lsdf[colname].apply(avgme)
    lsdf['range_distscr'] = lsdf[colname].apply(rangefinder)
    lsdf['max_dist_scr'] = lsdf[colname].apply(lambda x: max(x))
    print("group df shape: ", lsdf.shape)
    print(lsdf.head(1))
    return lsdf


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
    # mismaplines_dynamic(gall) # checkColumnValues(gall, 'dynamic_slope_scores')
    # MAP dynamic slope score col back to merged file
    #all_ref_slope = dict(zip(gall.ID, gall.dynamic_slope_scores))
    #mall['dynamic_slope_scores'] = mall['ID']
    #mall.dynamic_slope_scores = mall.dynamic_slope_scores.map(all_ref_slope)



########### Large file processing - K markdown #######

def pandachunks_append(filename, chunksiz):
    # read the large csv file with specified chunksize 
    df_chunk = pd.read_csv(filename, chunksize=chunksiz)
    chunk_list = []  # Each chunk is in df format
    for chunk in df_chunk:  
        chunk_filter = chunk_preprocessing(chunk)
        chunk_list.append(chunk_filter)
    df_concat = pd.concat(chunk_list)

def pandachunks(in_f, out_f, colnames, chunksize):
    # saves each chunk, does not append
    df_chunk = pd.read_csv(in_f, low_memory=False, chunksize=chunksize)
    for chunk in df_chunk:  
        chunk.columns = colnames 
        chunk_filter = pd.merge(chunk,non19,how='inner',on=['pos_id19'])
        chunk_filter.to_csv(out_f, index=False, header=False, mode='a')
    # pandachunks('not_detected_merge_appended.csv', 'FILTERchunks.csv', colnames, 10000)

############### BED file functions #############

def formatHGVSvep(element):
    # input = "c.960T>G;c.951T>G;c.1062T>G"
    # output = 'T/G'
    bun = re.findall(r"c.\d+(\D>\D)",element)
    str1 = ''.join(bun[0])
    fin = str1.replace('>','/')
    return fin

def make_bedfile(): 
    detected = "SCORE_annotation_3840_CYS_LYS_detected.csv"
    notdet = "SCORE_annotation_3840_CYS_LYS_NOT_allcol_detected.csv"
    dect = pd.read_csv(detected, usecols=["pos_id19", "hg19_chr", "hg19_pos(1-based)", "Amino_acids", "HGVSc_VEP", "cds_strand", 'pos_ID'], low_memory=False)
    notdect = pd.read_csv(notdet, usecols=["pos_id19", "hg19_chr", "hg19_pos(1-based)", "Amino_acids", "HGVSc_VEP", "cds_strand", 'pos_ID_falseCKtarget'], low_memory=False)
    # change to int type
    notdect['hg19_pos(1-based)'] = notdect['hg19_pos(1-based)'].astype(int)
    dect['hg19_pos(1-based)'] = dect['hg19_pos(1-based)'].astype(int)
    # new unique ID for bed file and VEP QC
    dect['bedID'] = dect['pos_id19'] + ";" + dect['pos_ID'] + ";" + dect['Amino_acids']
    notdect['bedID'] = notdect['pos_id19'] + ";" + notdect['pos_ID_falseCKtarget'] + ";" + notdect['Amino_acids']
    # MISSENSE BEDFILE 
    # simplify df
    keepdet = ['hg19_chr', 'hg19_pos(1-based)','HGVSc_VEP', 'cds_strand', 'bedID']
    bed_dect = dect[keepdet].copy()
    bed_not = notdect[keepdet].copy()
    # # change order of columns
    cols = ['hg19_chr', 'hg19_pos(1-based)', 'HGVSc_VEP', 'cds_strand', 'bedID']
    bed_dect = bed_dect[cols]
    bed_not = bed_not[cols]
    bed_dect['end'] = bed_dect['hg19_pos(1-based)'].copy()
    bed_not['end'] = bed_not['hg19_pos(1-based)'].copy()
    # rename
    cols = ['hg19_chr', 'start', 'cds_change', 'cds_strand', 'bedID', 'end']
    bed_dect.columns = cols
    bed_not.columns = cols
    # change order again
    cols = ['hg19_chr', 'start', 'end', 'cds_change', 'cds_strand', 'bedID']
    bed_dect = bed_dect[cols]
    bed_not = bed_not[cols]
    bed_dect['cds_change'] = bed_dect['cds_change'].apply(formatHGVSvep)
    bed_not['cds_change'] = bed_not['cds_change'].apply(formatHGVSvep)
    bed_not.to_csv("BED_NOT_DETECTED_HGVSchange_1222911.tsv", sep='\t', index=False)
    bed_dect.to_csv("BED_DETECTED_HGVSchange_104475.tsv", sep='\t', index=False)

########## protein seq comparison ###################

def sequenceDifference(dfref, dfalt):
    # assumes ID and proSequence column
    # $$$
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


################### MISMAP SCORE PROJECT

def get_pos_dictionary(pdcol):
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
    #list_of_dicts = get_pos_dictionary(groupedCYSLYS2_sort['pos'])
    #{2556: 'C', 1280: 'C', 991: 'C', 838: 'K', 849: 'K', 1339: 'K', 1425: 'K', 1221: 'K', 1044: 'K', 300: 'K', 1232: 'K', 239: 'K'}


def identicalSequenceCheck(dfref, dfalt):
    # Pmap assumes ID and proSequence column
    # used for ukb to ukb 2018 2012 seq comparison dfref is 2018
    # assumes ID and proSequence column
    # compares the sequence in df and dict value, adds column based on comparison results
    # True is they are identical, False is they are different, no info on HOW diff they are
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


def sequenceDistance(dfEnsp, ref_dic, newcolresult, hamming, hammingNorm, levenshtein):
    # identicalSequence including distance scores
    # version that includes scores for TRUE sequences as well as FALSE (matching canonical uniprot seq)
    res = []
    ham = []
    hamnorm = []
    lev = []
    serSeq = dfEnsp['proSequence'].copy()
    serID = dfEnsp['UniProtSP_xref'].copy()
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
#             ham.append('identical')
#             hamnorm.append('identical')
#             lev.append('identical')
            ham.append(textdistance.hamming(mypep, p))
            # normalized hamming = # mismatched positions/ len of longer sequence
            hamnorm.append(textdistance.hamming.normalized_distance(mypep, p))
            # levenshtein score is edit based but not not penalized position, insertion at pos 1 is jsut 1 diff
            lev.append(textdistance.levenshtein(mypep,p))
        # not identical to canonical
        if mypep != p:
            res.append('False')
            # calculates hamming distance, penalizes positional differences, edit based distance
            ham.append(textdistance.hamming(mypep, p))
            # normalized hamming = # mismatched positions/ len of longer sequence
            hamnorm.append(textdistance.hamming.normalized_distance(mypep, p))
            # levenshtein score is edit based but not not penalized position, insertion at pos 1 is jsut 1 diff
            lev.append(textdistance.levenshtein(mypep,p))
    # add new column
    dfEnsp.loc[:, newcolresult] = res
    dfEnsp.loc[:, hamming] = ham
    dfEnsp.loc[:, hammingNorm] = hamnorm
    dfEnsp.loc[:, levenshtein] = lev
    return dfEnsp


####### COMBINE DF #####################

def combineDFS(optio):
    # combined for input into R
    # both true and false identity scores
    if optio == "all":
        df85 = pd.read_csv("groupedMISMAP_score_85_3979.csv")
        df92 = pd.read_csv("groupedMISMAP_score_92_3979.csv")
        df94 = pd.read_csv("groupedMISMAP_score_94_3979.csv")
        df96 = pd.read_csv("groupedMISMAP_score_96_3979.csv")
        df97 = pd.read_csv("groupedMISMAP_score_97_3979.csv")
        df85['release'] = 85
        df92['release'] = 92
        df94['release'] = 94
        df96['release'] = 96
        df97['release'] = 97
        # combining dfs on axis 
        df_merge = pd.concat([df85, df92, df94, df96, df97])
        print(df_merge.head(10))
        print(df_merge.tail(10))
        print(df_merge.shape)
        print("all identity scores len: ", 3979*5)
        df_merge.to_csv("Rmerge_MISMAP_AllIdentityOnly_3979.csv",index=False)
    if optio == "false":
        #only false
        df85 = pd.read_csv("groupedMISMAP_score_85_FALSEidentity_1805.csv")
        df92 = pd.read_csv("groupedMISMAP_score_92_FALSEidentity_1805.csv")
        df94 = pd.read_csv("groupedMISMAP_score_94_FALSEidentity_1805.csv")
        df96 = pd.read_csv("groupedMISMAP_score_96_FALSEidentity_1805.csv")
        df97 = pd.read_csv("groupedMISMAP_score_97_FALSEidentity_1805.csv")
        df85['release'] = 85
        df92['release'] = 92
        df94['release'] = 94
        df96['release'] = 96
        df97['release'] = 97
        # combining dfs on axis
        df_merge = pd.concat([df85, df92, df94, df96, df97])
        print(df_merge.head(10))
        print(df_merge.tail(10))
        print(df_merge.shape)
        print("only false identity scores len: ", 1805*5)
        df_merge.to_csv("Rmerge_MISMAP_FalseIdentityOnly_1805.csv",index=False)
    #     combineDFS("all")
    #     combineDFS("false")




def Pmap_merge():
    infalse = ["CADDmapped_falsetarget_chr10.csv","CADDmapped_falsetarget_chr11.csv","CADDmapped_falsetarget_chr12.csv","CADDmapped_falsetarget_chr13.csv","CADDmapped_falsetarget_chr14.csv","CADDmapped_falsetarget_chr15.csv","CADDmapped_falsetarget_chr16.csv","CADDmapped_falsetarget_chr17.csv","CADDmapped_falsetarget_chr18.csv","CADDmapped_falsetarget_chr19.csv","CADDmapped_falsetarget_chr1.csv","CADDmapped_falsetarget_chr20.csv","CADDmapped_falsetarget_chr21.csv","CADDmapped_falsetarget_chr22.csv","CADDmapped_falsetarget_chr2.csv","CADDmapped_falsetarget_chr3.csv","CADDmapped_falsetarget_chr4.csv","CADDmapped_falsetarget_chr5.csv","CADDmapped_falsetarget_chr6.csv","CADDmapped_falsetarget_chr7.csv","CADDmapped_falsetarget_chr8.csv","CADDmapped_falsetarget_chr9.csv","CADDmapped_falsetarget_chrX.csv","CADDmapped_falsetarget_chrY.csv"]
    mergeme = ["addedkey_ckfiltered_chr10.csv","addedkey_ckfiltered_chr11.csv","addedkey_ckfiltered_chr12.csv","addedkey_ckfiltered_chr13.csv","addedkey_ckfiltered_chr14.csv","addedkey_ckfiltered_chr15.csv","addedkey_ckfiltered_chr16.csv","addedkey_ckfiltered_chr17.csv","addedkey_ckfiltered_chr18.csv","addedkey_ckfiltered_chr19.csv","addedkey_ckfiltered_chr1.csv","addedkey_ckfiltered_chr20.csv","addedkey_ckfiltered_chr21.csv","addedkey_ckfiltered_chr22.csv","addedkey_ckfiltered_chr2.csv","addedkey_ckfiltered_chr3.csv","addedkey_ckfiltered_chr4.csv","addedkey_ckfiltered_chr5.csv","addedkey_ckfiltered_chr6.csv","addedkey_ckfiltered_chr7.csv","addedkey_ckfiltered_chr8.csv","addedkey_ckfiltered_chr9.csv","addedkey_ckfiltered_chrX.csv","addedkey_ckfiltered_chrY.csv"]
    chrorder = ["chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr1","chr20","chr21","chr22","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"]
    for ii,ll,mm in zip(infalse, chrorder, mergeme):
        df = pd.read_csv(ii)
        merg = pd.read_csv(mm, low_memory=False)
        chrlabel = ll


##### from R prep 10.18.19 markdown ##########
def assemblyGrouper(df, CK):
    # cleaning CK specific dbNSFP filtered files with cadd scores from hg19 and hg38
    # spliting df to add assembly column (doubles in size) for R figure
    # group labels are helpfl for parsing out later
    
    if CK == 'C':       
        test2 = df[['pos_id19', 'pos_id38', 'xref', 'matched_target', 'lost_amino', 'gained_amino',
       'Amino_acids', 'pos_ID', 'Cys_reactivity', 'Cys_react_threshold',
       'Cys_target_label', 'geneNamePrimary']].copy()
        test3 = test2.copy()
        test2['CADD_score'] = df['CADD_phred_hg19'].copy()
        test2['assembly'] = 'hg19'
        test3['CADD_score'] = df['CADD_phred_hg38'].copy()
        test3['assembly'] = 'hg38'
        test4 = pd.concat([test2,test3])
        test4.reset_index(inplace=True,drop=True)
        print("combined shape : ", test4.shape)
        print("copy shapes : ", test2.shape, test3.shape)
        return test4

    if CK == 'K':
        test2 = df[['pos_id19', 'pos_id38', 'xref', 'matched_target', 'lost_amino', 'gained_amino',
       'Amino_acids', 'pos_ID', 'Lys_reactivity', 'Lys_react_threshold','Lys_target_label','geneNamePrimary']].copy()
        test3 = test2.copy()
        test2['CADD_score'] = df['CADD_phred_hg19'].copy()
        test2['assembly'] = 'hg19'
        test3['CADD_score'] = df['CADD_phred_hg38'].copy()
        test3['assembly'] = 'hg38'   
        test4 = pd.concat([test2,test3])
        test4.reset_index(inplace=True,drop=True)
        print("combined shape : ", test4.shape)
        print("copy shapes : ", test2.shape, test3.shape)
        return test4

    if CK == 'N':
        test2 = df[['pos_id19', 'pos_id38', 'xref', 'matched_target','lost_amino', 'gained_amino',
       'Amino_acids', 'pos_ID_falseCKtarget', 'geneNamePrimary']].copy()
        test3 = test2.copy()
        test2['CADD_score'] = df['CADD_phred_hg19'].copy()
        test2['assembly'] = 'hg19'
        test3['CADD_score'] = df['CADD_phred_hg38'].copy()
        test3['assembly'] = 'hg38'   
        test4 = pd.concat([test2,test3])
        test4.reset_index(inplace=True,drop=True)
        print("combined shape : ", test4.shape)
        print("copy shapes : ", test2.shape, test3.shape)
        return test4

        simdf = df[["pos_id19","pos_id38","aaref","aaalt","matched_aapos","CADD_phred_hg38","matched_UKBID","CADD_phred_hg19","CADDdiff_38minus19","Amino_acids","pos_ID_falseCKtarget"]].copy()

        # merging simplified dataframe with uniprot position dict
        mer = pd.merge(simdf, merg, how='inner', on=['pos_id19'])
        print('shape simplified df pre merge : ', simdf.shape)
        print('shape of merged file: ', mer.shape)

        Foutname = "CADDmapped_falsetargets_merged_" + chrlabel + ".csv"
        mer.to_csv(Foutname, index=False)
        print("done with: ", Foutname)

