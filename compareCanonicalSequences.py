#!/usr/bin/env python3
# -*- coding: utf-8 -*-
## author: maria f. palafox

"""
Compares UniProtKB FASTA files from 2 different database releases and marks the
amino acid and position differences between the peptides. Comparisons
are made using the index range of the shorter peptide and starting with the
first amino acid of both sequences.

Added columns description:
    - comparedFASTALen: Length of sequence with identical UniProtKB ID in the
    other FASTA file [A]

    - identical2Canonical: True/False for 100% seqeunce identity between the
    two versions of UniProtKB proteins

    - lengthDiff: Difference in length of protein from [B] FASTA - [A] FASTA

    - residueDiff: Amino acid difference between compared sequences using
    same sequence index

    - positionDiff: Position of differing amino acid between compared
    seqeunces using an index range based on the shorter sequence

Below is the sequence for UniProtKB ID O15392, Baculoviral IAP
repeat-containing protein 5, from 2018_06 and 2012_01 UniProtKB releases.

## [A] FASTA from 2018_06 UKB release for O15392 canonical sequence
MGAPTLPPAWQPFLKDHRISTFKNWPFLEGCACTPERMAEAGFIHCPTENEPDLAQCFFCFKELEGWEPDDDPIEEHKKHSSGCAFLSVKKQFEELTLGEFLKLDRERAKNKIAKETNNKKKEFEETA*K*VRRAIEQLAAMD

## [B] FASTA from 2012_01 UKB release for O15392 canonical sequence
MGAPTLPPAWQPFLKDHRISTFKNWPFLEGCACTPERMAEAGFIHCPTENEPDLAQCFFCFKELEGWEPDDDPIEEHKKHSSGCAFLSVKKQFEELTLGEFLKLDRERAKNKIAKETNNKKKEFEETA*E*KVRRAIEQLAAMD

NOTE:
    Comparison results are relative to the [B] input FASTA sequence. The shorter
    peptide length is set as the range of index positions to loop over both
    sequences. For protein O15392, position 129 in FASTA-A has amino acid 'K' and FASTA-B has amino acid 'E', therefore the difference is noted as 'K' for 129 position in the saved results file.

"""
from maplib import *
import argparse
parser = argparse.ArgumentParser(description="Compares sequences from "
                                             "UniProtKB releases to identify "
                                             "position and amino acid "
                                             "differences between ID-matched "
                                             "proteins.")
parser.add_argument("fastaA",
                    help="UniProtKB FASTA filename used as "
                                    "annotation reference proteome.")
parser.add_argument("fastaB",
                    help="UniProtKB FASTA filename used as meta-data "
                                  "for chemoproteomic-detected residues.")
parser.add_argument("idFilter",
                    help="A .csv file with first column containing UniProtKB "
                         "IDs used to filter FASTA files. Protein IDs "
                         "contained in this file and also in both FASTA "
						"will be used for sequence comparison.")
parser.add_argument("outfile",
                    help="name of new chemoproteomic results .csv w/ column "
                         "for FASTA sequence comparison results.")
args = parser.parse_args()
fastaA = args.fastaA
fastaB = args.fastaB
ids = args.idFilter
outfile = args.outfile

def compareCanonical(fastaA, fastaB, idFilter, outfile):
    dfA = uniprotkbFasta2csv(fastaA)
    dfB = uniprotkbFasta2csv(fastaB)
    ids = pd.read_csv(idFilter)
    # intersect between all files IDs in  first column
    idA = list(set(dfA.iloc[:,0]))
    idB = list(set(dfB.iloc[:,0]))
    ids = list(set(ids.iloc[:,0]))
    overlap1 = intersectionCheck(idA, idB)
    overlap2 = intersectionCheck(overlap1, ids)
    dfB = addcolumnconditionalDrop(overlap2, dfB, 'ID', 'filtercol')
    dfA = addcolumnconditionalDrop(overlap2, dfA, 'ID', 'filtercol')
    # make dict from ref IDs and sequences
    pepdict = dict(zip(dfA.ID, dfA.proSequence))
    # new columns
    otherLen = []
    identical = []
    BminusA = []
    diffcol = []
    posdiffcol = []
    # loop over reference FASTA rows
    for index, row in dfB.iterrows():
        pepB = row['proSequence'] # save ref pep var
        ukbID = row['ID'] # save ukb id
        pepA = pepdict[ukbID] # save pep in dictionary matching ukb id
        str(pepB)
        str(pepA)
        lenB = len(pepB)
        lenA = len(pepA)
        # 100% identical case
        if lenA == lenB:
            if pepA == pepB:
                otherLen.append(lenB)
                identical.append(True)
                BminusA.append(0)
                diffcol.append("0")
                posdiffcol.append("0")
            # same length not identical case
            else:
                # in reference to A FASTA
                output_list = [pepA[i] for i in range(len(pepB)) if pepB[i]
                               != pepA[i]]
                posoutput_list2 = [i+1 for i in range(len(pepB)) if pepB[i]
                                   != pepA[i]]
                out = ', '.join(output_list)
                otherLen.append(lenB)
                identical.append(False)
                BminusA.append(0)
                diffcol.append(out)
                posdiffcol.append(posoutput_list2)
        # A longer than B case
        if lenA > lenB:
            output_list = [pepA[i] for i in range(len(pepB)) if pepB[i] !=
                           pepA[i]]
            posoutput_list2 = [i+1 for i in range(len(pepB)) if pepB[i] !=
                               pepA[i]]
            out = ', '.join(output_list)
            otherLen.append(lenA)
            identical.append(False)
            BminusA.append((lenB-lenA)) # A reference -  order
            diffcol.append(out)
            posdiffcol.append(posoutput_list2)
        # A shorter than B case
        if lenA < lenB:
            output_list = [pepB[i] for i in range(len(pepA)) if pepA[i] !=
                           pepB[i]]
            posoutput_list2 = [i+1 for i in range(len(pepA)) if pepA[i] != pepB[i]]
            out = ', '.join(output_list)
            otherLen.append(lenA)
            identical.append(False)
            BminusA.append((lenB-lenA))
            diffcol.append(out)
            posdiffcol.append(posoutput_list2)
    # FASTA seq length compared to reference (A)
    dfB.loc[:, 'comparedFASTALen'] = otherLen
    # boolean
    dfB.loc[:, 'identical2Canonical'] = identical
    # length difference
        ## based on reference (A) always
    dfB.loc[:, 'lengthDiff'] = BminusA
    # residue difference
        ## based on shorter peptide
    dfB.loc[:, 'residueDiff'] = diffcol
    # position difference 
        ## based on shorter peptide
    dfB.loc[:, 'positionDiff'] = posdiffcol
    checkColumnValues(dfB, "identical2Canonical")
    # drop proSequence column
    # dfB.drop("proSequence", 1, inplace=True)
    dfB.to_csv(outfile, index=False)
    print("*** results saved as ", outfile, " ***")
    print("*** shape of results file: ", dfB.shape)
    print()
    print(dfB.head(3))
    return dfB

if __name__ == '__main__':
	compareCanonical(fastaA, fastaB ids, outfile)