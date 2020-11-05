#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# protocol_identicalProSequenceChecker.py
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



parser = argparse.ArgumentParser()

parser.add_argument("EnsemblFasta", help="This is the Ensembl release-specific fasta filename, example: Homo_sapiens.GRCh38.pep.all.fa")

parser.add_argument("StableIDkey", help="This is the Ensembl-UniProt stable ID cross-reference filename from the same release as the Ensembl fasta file provided, example: Homo_sapiens.GRCh38.94.uniprot.tsv")

parser.add_argument("EnsemblRelease", help="This is the name of specific Ensembl release (ex. v94) that match the fasta and xref file")

parser.add_argument("UniprotFasta", help="This is the UniProtKB fasta filename")

parser.add_argument("UniprotDatabase", help="This is the UniProt database name, options are 'SP' for SWISSPROT or 'SPTREMBL' for SWISSPROT & TREMBL")

args = parser.parse_args()
ensemblFasta = args.EnsemblFasta
xrefile = args.StableIDkey
release = args.EnsemblRelease
ukbfa = args.UniprotFasta
UKBdb = args.UniprotDatabase 
outfileName = release + "CheckedProSeq_ENSPtoUKB.csv"



def checkColumnValues(df, col):
    print(df[col].value_counts().reset_index().rename(columns={'index':col, col:'Count'}))

def addColumnLabelDrop(mapList, df, dfcol, newcol):
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
    print("total values in column: ", len(df[colname]))
    print("total unique values: ", len(df[colname].unique()))
    print()

def ENSPfasta2DF(filename, release):
    # Pmap_ensembl_fasta.py - parse sequence fasta file
    identifiers = [seq_record.id for seq_record in SeqIO.parse(filename, "fasta")]
    descr = [seq_record.description for seq_record in SeqIO.parse(filename, "fasta")]
    lengths = [len(seq_record.seq) for seq_record in SeqIO.parse(filename, "fasta")]
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
        nA = pd.Series(assembly, name='Assembly')
        nC = pd.Series(chrom, name= 'chromosome')
        nSta = pd.Series(DNAstart, name='start')
        nSto = pd.Series(DNAstop, name='stop')
        s2 = pd.Series(enspv, name='ENSPv')
        s1 = pd.Series(ensp, name= 'ENSP')
        s7 = pd.Series(lengths, name='ENSPlength', dtype=int)
        s8 = pd.Series(proseq, name='proSequence')
        s4 = pd.Series(enstv, name='ENSTv')
        s3 = pd.Series(enst, name= 'ENST')
        s6 = pd.Series(ensgv, name='ENSGv')
        s5 = pd.Series(ensg, name= 'ENSG')
        series = [s1, s2, s7, s8, s3, s4, s5, s6, nA, nC, nSta, nSto]
        df = pd.concat(series, axis=1)
        df.dropna(inplace=True)
        df.reset_index(drop=True, inplace=True)
        return(df)
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
                    HGNCdescription.append(restOfRow)
        nS = pd.Series(HGNCsymbol, name='HGNCsymbol')
        nD = pd.Series(HGNCdescription, name='HGNCdescription')
        nA = pd.Series(assembly, name='Assembly')
        nC = pd.Series(chrom, name= 'chromosome')
        nSta = pd.Series(DNAstart, name='start')
        nSto = pd.Series(DNAstop, name='stop')
        s2 = pd.Series(enspv, name='ENSPv')
        s1 = pd.Series(ensp, name= 'ENSP')
        s7 = pd.Series(lengths, name='ENSPlength', dtype=int)
        s8 = pd.Series(proseq, name='proSequence')
        s4 = pd.Series(enstv, name='ENSTv')
        s3 = pd.Series(enst, name= 'ENST')
        s6 = pd.Series(ensgv, name='ENSGv')
        s5 = pd.Series(ensg, name= 'ENSG')
        series = [s1, s2, s7, s8, s3, s4, s5, s6, nA, nC, nSta, nSto, nS, nD]
        df = pd.concat(series, axis=1)
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




def main():
    print("Importing Ensembl fasta")
    enspDF = ENSPfasta2DF(ensemblFasta, release)
    print("Shape of Ensembl fasta dataframe: ",enspDF.shape)

    print("Importing UniProt fasta")
    ukbDF = UKBfasta2DF(ukbfa)
    print("Shape of UniProt fasta dataframe: ",ukbDF.shape)
    listUKB = list(ukbDF['ID'])

    xref = pd.read_table(xrefile)
    if UKBdb == "SP":
        xref = xref[xref['db_name']== 'Uniprot/SWISSPROT'].copy() 
        print("Xref SP UniProt IDs kept")
    if UKBdb == "SPTREMBL":
        print("Xref SP and TREMBL UniProt IDs kept")

    ## filter ensembl xref file by ukbIDs list
    print("filtering Ensembl xref file by UniProt fasta stable IDs")
    xref = addColumnLabelDrop(listUKB, xref, 'xref', 'xrefInUKBfasta')

    ## filter ensembl fasta for ENSP in xref file (filtered by uniprot fasta)
    xref = xref[['protein_stable_id', 'xref']].copy()
    xref.columns = ['ENSP', 'xref']
    listENSP = list(set(xref['ENSP']))
    print("filtering  Ensembl peptide df by stable ENSP IDs from xref file")
    enspDF = addColumnLabelDrop(listENSP, enspDF, 'ENSP', 'enspInXref')       

    ## merge xref file (with ENSP & UKB ID) with Ensembl fasta (with ENSP, ENST, ENSG, proSeq...) 
    print("combining Ensembl xref (ENSP - UKBID cross-ref) with Ensembl fasta DF (contains Ensembl biotype IDs both stable and versioned)")
    mergeENSP = pd.merge(xref, enspDF, on=['ENSP'])

    ## check ensembl pro sequence identity against canonical ukb seq from fasta
    print("identify ensembl protein sequences identical to UniProt reference")
    finaldf = identicalSequenceCheck(ukbDF, mergeENSP)

    # change colnames 
    if release == "v85":
            finaldf.columns = ['ENSP', 'UKBID', 'ENSPv', 
            'Length', 'proSequence', 
            'ENST', 'ENSTv','ENSG', 'ENSGv', 'Assembly', 'chromosome', 
            'start', 'stop', 'enspInXref', 'identical2UKBfasta']
            # simplify
            keepme = ['ENSP', 'UKBID', 'ENSPv', 'Length', 'proSequence', 
            'ENST', 'ENSTv','ENSG', 'ENSGv', 'Assembly', 'chromosome', 
            'start', 'stop']
            finaldf = finaldf[keepme].copy()
    else:
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
    print("Dropping rows with nonstandard chromosome names")
    finaldf= finaldf[finaldf['chromosome'].isin(chrlist)]

    print("********** Summary of saved script output **********")
    checkColumnValues(finaldf, 'chromosome')
    print("Ensembl Protein IDs")
    uniqueCount(finaldf, 'ENSP')
    print("Ensembl Transcript IDs")
    uniqueCount(finaldf, 'ENST')
    print("Ensembl Gene IDs")
    uniqueCount(finaldf, 'ENSG')
    print("UniProt IDs")
    uniqueCount(finaldf, 'UKBID')

    # save df
    finaldf.to_csv(outfileName, index=False)

if __name__ == '__main__':
    main()


