#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# processUniprotFasta.py

"""
UPDATES:
1. added flag options depending on source of fasta (uniprot, refseq, gencode)
    7.24.19
    7.29.19
"""
import argparse
from Bio import SeqIO
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("fastaFile", help="This is the Input file .fasta that if going to be converted into DF ")

parser.add_argument("parsedDFcsv", help="this is name of output file that is .fasta --> DF ")

parser.add_argument("dbTYPE", help="options are gencode fasta ,uniprot, or refseq fasta")

args = parser.parse_args()
filename = args.fastaFile
parsedOut = args.parsedDFcsv
dbtype = args.dbTYPE

# this function is in library.py work on it there
def fastaToDF(filename, dbtype):
    # parse sequence fasta file
    identifiers = [seq_record.id for seq_record in SeqIO.parse(filename, "fasta")]
    lengths = [len(seq_record.seq) for seq_record in SeqIO.parse(filename, "fasta")]
    proSeq = [seq_record.seq for seq_record in SeqIO.parse(filename, "fasta")]

    # splitting identifiers into 2 seperate list
    if dbtype == "uniprot":
        splitAcc = []
        splitEntry = []
        for id in identifiers:
            splitID = id.split('|')
            acc = splitID[1]
            splitAcc.append(acc)
            entryName = splitID[2]
            splitEntry.append(entryName)
        # converting lists to pandas Series
        s1 = pd.Series(splitAcc, name='ID')
        s2 = pd.Series(splitEntry, name= 'entryName')
        s3 = pd.Series(lengths, name='Length')
        s4 = pd.Series(proSeq, name='proSequence')
        # column header = ,ID,Length,proSequence
        series = [s1,s2,s3,s4]
        df = pd.concat(series, axis=1)
        # parsed_df = pd.DataFrame(df.ID.str.split('|',1).tolist(),columns = ['geneID','feature'])
        # Newdf = pd.concat([parsed_df,remains] ,1)
        # pd.concat([v, df], 1)
        return(df)

    if dbtype == "gencode":
        s_txid = []
        s_geneid = []
        s_genesymbol = []
        for id in identifiers:
            splitID = id.split('|')
            txid = splitID[0]
            s_txid.append(txid)
            geneid = splitID[1]
            s_geneid.append(geneid)
            genesymbol = splitID[5]
            s_genesymbol.append(genesymbol)
        # converting lists to pandas Series
        s1 = pd.Series(s_txid, name='ENST')
        s2 = pd.Series(s_geneid, name= 'ENSG')
        s3 = pd.Series(s_genesymbol, name='geneSymbol')
        s4 = pd.Series(lengths, name='Length')
        s5 = pd.Series(proSeq, name='proSequence')
        series = [s1, s2, s3, s4, s5]
        df = pd.concat(series, axis=1)
        return(df)

    #  >NP_000005.2 alpha-2-macroglobulin isoform a precursor [Homo sapiens]
    if dbtype == "refseq":
        s_txid = []
        for id in identifiers:
            splitID = id.split()
            txid = splitID[0]
            s_txid.append(txid)
        s1 = pd.Series(s_txid, name='RefSeq_NP')
        s4 = pd.Series(lengths, name='Length')
        s5 = pd.Series(proSeq, name='proSequence')
        series = [s1, s4, s5]
        df = pd.concat(series, axis=1)
        return(df)

def main():
    df = fastaToDF(filename, dbtype)
    print(df.columns)
    df.to_csv(parsedOut, index=False)

if __name__ == '__main__':
    main()
