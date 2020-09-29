#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# identicalProSequenceChecker.py
# demonstration of script at stableID_proSequence_checker.ipynb

import os
import sys
import argparse
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from ast import literal_eval
import difflib

# CHANGE PATH
os.chdir("/Users/mariapalafox/Desktop/MAPPING/BEST_PRACTICE")
#-------------------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument("ensemblFasta", help="This is the Ensembl release-specific input fasta, example: v94Homo_sapiens.GRCh38.pep.all.fa")

parser.add_argument("ukbfasta", help="This is the UniProtKB input fasta")

parser.add_argument("IDxref", help="This is the Ensembl cross-reference file specific for mapping to UniProt IDs, example: Homo_sapiens.GRCh38.94.uniprot.tsv")

parser.add_argument("outfileName", help="This is the name of saved .csv")

args = parser.parse_args()
ensemblFasta = args.ensemblFasta
ukbfa = args.ukbfasta
xrefile = args.IDxref
outfileName = args.outfileName

#------------------ funx -------------------------------
def checkColumnValues(df, col):
    print(df[col].value_counts().reset_index().rename(columns={'index':col, col:'Count'}))

def addColumnLabelDrop(mapList, df, dfcol, newcol):
    # same as above but also drops False value rows
    print("mapList set length: ", len(set(mapList)))
    mendel = []
    for g in df[dfcol]:
        if g in mapList:
            mendel.append("True")
        else:
            mendel.append("False")
    df.loc[:, newcol] = mendel
    checkColumnValues(df, newcol)
    print("dropping False rows")
    df.drop(df[df[newcol] == "False"].index, inplace=True)
    df.reset_index(inplace=True, drop=True)
    print("final shape: ", df.shape)
    return df

def uniqueCount(df, colname):
    print("len of col: ", len(df[colname]))
    print("len of col set: ", len(df[colname].unique()))
    print()

def ENSPfasta2DF(filename):
    # Pmap_ensembl_fasta.py - parse sequence fasta file
    identifiers = [seq_record.id for seq_record in SeqIO.parse(filename, "fasta")]
    descr = [seq_record.description for seq_record in SeqIO.parse(filename, "fasta")]
    lengths = [len(seq_record.seq) for seq_record in SeqIO.parse(filename, "fasta")]
    proSeq = [seq_record.seq for seq_record in SeqIO.parse(filename, "fasta")]
    # create ensembl id and seq list objects
    ensp = []
    enspv = []
    proseq = []
    ensg = []
    ensgv = []
    enst = []
    enstv = []
    # create description info cols
    assembly = []
    chrom = []
    DNAstart = []
    DNAstop = []
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
    # loop over other row info
    for row in descr:
        splitrow = row.split(" ")
        for i,val in enumerate(splitrow):
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
                HGNCdescription.append(restOfRow) # adding all info
    # add description cols         
    nS = pd.Series(HGNCsymbol, name='HGNCsymbol')
    nD = pd.Series(HGNCdescription, name='HGNCdescription')
    nA = pd.Series(assembly, name='Assembly')
    nC = pd.Series(chrom, name= 'chromosome')
    nSta = pd.Series(DNAstart, name='start')
    nSto = pd.Series(DNAstop, name='stop')
    # ensembl id type cols and seq
    s2 = pd.Series(enspv, name='ENSPv')
    s1 = pd.Series(ensp, name= 'ENSP')
    s7 = pd.Series(lengths, name='ENSPlength', dtype=int)
    s8 = pd.Series(proseq, name='proSequence')
    s4 = pd.Series(enstv, name='ENSTv')
    s3 = pd.Series(enst, name= 'ENST')
    s6 = pd.Series(ensgv, name='ENSGv')
    s5 = pd.Series(ensg, name= 'ENSG')
    # create df and return
    series = [s1, s2, s7, s8, s3, s4, s5, s6, nA, nC, nSta, nSto, nS, nD]
    df = pd.concat(series, axis=1)
    print("dropping None value rows in df")
    df.dropna(inplace=True)
    df.reset_index(drop=True, inplace=True)
    return(df)

def UKBfasta2DF(filename):
    identifiers = [seq_record.id for seq_record in SeqIO.parse(filename, "fasta")]
    lengths = [len(seq_record.seq) for seq_record in SeqIO.parse(filename, "fasta")]
    proSeq = [seq_record.seq for seq_record in SeqIO.parse(filename, "fasta")]
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
    s3 = pd.Series(lengths, name='UKBlength')
    s4 = pd.Series(proseq, name='proSequence')
    series = [s1,s2,s3,s4]
    df = pd.concat(series, axis=1)
    return(df)

def identicalSequenceCheck(dfref, dfalt):
    # Pmap_protein_sequence_comparison.py
    # dfref must have 'ID' and 'proSequence' colnames
    # dfalt must have 'xref' and 'proSequence' colnames
    ref_dic = dict(zip(dfref.ID, dfref.proSequence)) # from uniprot fasta
    newcol = []
    for index, row in dfalt.iterrows(): # ensembl df
        pep = row['proSequence']
        ukb_id = row['xref'] # colname
        # retrieve uniprot fasta seq using dictionary
        pepRef = ref_dic[ukb_id]
        str(pepRef)
        str(pep)
        if pepRef == pep:
            newcol.append('True')
        if pepRef != pep: 
            newcol.append('False')
    # add identity score list as new column
    dfalt.loc[:,'identical2UKBfasta'] = newcol
    print(dfalt.shape)
    checkColumnValues(dfalt, "identical2UKBfasta")
    print("dropping non identical protein seq")
    dfalt = dfalt[dfalt['identical2UKBfasta'] == "True"].copy()
    print(dfalt.shape)
    dfalt.reset_index(inplace=True, drop=True)
    return dfalt
#-------------------------------------------------------------

def main():
    ## import ensembl fasta
    enspDF = ENSPfasta2DF(ensemblFasta)

    ## import ukb fasta 
    ukbDF = UKBfasta2DF(ukbfa)

    ## list of ukbIDs to filter xref
    listUKB = list(ukbDF['ID'])

    ## import xref 
    xref = pd.read_table(xrefile)

    ## filter ensembl xref file by ukbIDs list
    xref = addColumnLabelDrop(listUKB, xref, 'xref', 'xrefInUKBfasta')

    ## filter ensembl fasta for only ENSP contained in xref file 
        # xref was filtered for UKBIDs in ref proteome
    xref = xref[['protein_stable_id', 'xref']].copy()
    xref.columns = ['ENSP', 'xref']
    listENSP = list(set(xref['ENSP']))
    enspDF = addColumnLabelDrop(listENSP, enspDF, 'ENSP', 'enspInXref')

    ## merge ensembl files based on ENSP- 
        # xref was filtered by list of UniProt IDs from UKB fasta
        # enspDF is the ensembl fasta df filtered for ENSP contained in filtered xref file
    mergeENSP = pd.merge(xref, enspDF, on=['ENSP'])

    ## check ensembl pro sequence identity against canonical ukb seq from fasta
    finaldf = identicalSequenceCheck(ukbDF, mergeENSP)

    # change colnames 
    finaldf.columns = ['ENSP', 'UKBID', 'ENSPv', 'Length', 'proSequence', 
    'ENST', 'ENSTv','ENSG', 'ENSGv', 'Assembly', 'chromosome', 
    'start', 'stop','HGNCsymbol', 'HGNCdescription', 'enspInXref', 
    'identical2UKBfasta']
    # simplify
    keepme = ['ENSP', 'UKBID', 'ENSPv', 'Length', 'proSequence', 
    'ENST', 'ENSTv','ENSG', 'ENSGv', 'Assembly', 'chromosome', 
    'start', 'stop','HGNCsymbol', 'HGNCdescription']
    finaldf = finaldf[keepme].copy()

    # removing rows that dont have std chrom coordinate from Ensembl fasta file
    chrlist = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    finaldf= finaldf[finaldf['chromosome'].isin(chrlist)]
    print("------Summary of final df------------")
    checkColumnValues(finaldf, 'chromosome')
    uniqueCount(finaldf, 'ENSP')
    uniqueCount(finaldf, 'ENST')
    uniqueCount(finaldf, 'ENSG')
    uniqueCount(finaldf, 'UKBID')
    # save df
    finaldf.to_csv(outfileName, index=False)


if __name__ == '__main__':
    main()

