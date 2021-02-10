#!/usr/bin/env python3
# -*- coding: utf-8 -*-
## author: maria f. palafox
# demonstration of script at identical_proSequence_checker.html

from maplib import *
import argparse

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
    print("*** Importing Ensembl fasta ***")
    enspDF = ENSPfasta2DF(ensemblFasta, release)
    print("Shape of Ensembl fasta dataframe: ",enspDF.shape)

    print("*** Importing UniProt fasta ***")
    ukbDF = uniprotkbFasta2csv(ukbfa)
    print("Shape of UniProt fasta dataframe: ",ukbDF.shape)
    listUKB = list(ukbDF['ID'])

    xref = pd.read_table(xrefile)
    if UKBdb == "SP":
        xref = xref[xref['db_name']== 'Uniprot/SWISSPROT'].copy() 
        print("Xref SP UniProt IDs kept")
    if UKBdb == "SPTREMBL":
        print("Xref SP and TREMBL UniProt IDs kept")

    ## filter ensembl xref file by ukbIDs list
    print("*** filtering Ensembl xref file by UniProt fasta stable IDs ***")
    xref = addcolumnconditionalDrop(listUKB, xref, 'xref', 'xrefInUKBfasta')

    ## filter ensembl fasta for ENSP in xref file (filtered by uniprot fasta)
    xref = xref[['protein_stable_id', 'xref']].copy()
    xref.columns = ['ENSP', 'xref']
    listENSP = list(set(xref['ENSP']))
    print("*** filtering  Ensembl peptide df by stable ENSP IDs from xref "
          "file ***")
    enspDF = addcolumnconditionalDrop(listENSP, enspDF, 'ENSP', 'enspInXref')

    ## merge xref file (with ENSP & UKB ID) with Ensembl fasta (with ENSP, ENST, ENSG, proSeq...) 
    print("*** combining Ensembl xref (ENSP - UKBID cross-ref) with Ensembl "
          "fasta DF (contains Ensembl biotype IDs both stable and versioned) ")
    mergeENSP = pd.merge(xref, enspDF, on=['ENSP'])

    ## check ensembl pro sequence identity against canonical ukb seq from fasta
    print("*** identify ensembl protein sequences identical to UniProt "
          "reference ***")
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


