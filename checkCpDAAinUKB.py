#!/usr/bin/env python3
# -*- coding: utf-8 -*-
## author: maria f. palafox

"""
Code used to check published residue-level
chemoproteomic results lacking meta-data for reference sequences against a
user-provided release of UniProtKB canonical sequences. The original search
database from the 2012 chemoproteomic studies were not accessible online,
so the following code  was used to find the best
available UniProtKB release with canonical sequences that best match
chemo-detected residue positions from published study results.

Prior to this quality control step, there were 4,119 protein IDs with detected
cysteines and lysines from previous studies. After checking residues against
an available UniProtKB release (2012_11, same year as search database used by
authors of chemoproteomic results), 26 protein IDs were dropped from our
analysis based on original CpDAA positions not overlapping positions in the
provided FASTA protein sequences. Total of 4,093 protein IDs remained for
mapping to the 2018 UniProtKB proteome we chose to serve as a reference for
all functional annotations and other chemoproteomic studies included.
"""

from maplib import *
import argparse
parser = argparse.ArgumentParser(description="Code used to check published "
			"residue-level chemoproteomic results lacking meta-data for reference "
            "sequences against a user-provided release of UniProtKB canonical sequences.")
parser.add_argument("fasta", help="FASTA filename used as reference proteome "
                                  "to check positions of "
                                  "chemoproteomic-detected residues from "
                                  "previous study.")
parser.add_argument("aaletter",
                    help="one letter amino acid (e.g. 'C' for Cysteine) to "
                         "make pos_ID key from using protein seqs in "
                         "fasta.")
parser.add_argument("cpdaa",
                    help="Residue-level chemoproteomic results .csv file w/ "
                         "added CpDAA key column named 'pos_ID' (e.g. "
                         "P11313_K170).")
parser.add_argument("outfile",
                    help="name of new chemoproteomic results .csv w/ column "
                         "for results of UniProtKB FASTA sequence check")
args = parser.parse_args()
fasta = args.fasta
aaletter = args.aaletter
cpdaa = args.cpdaa
outfile = args.outfile

def CpDAA_position_check(fasta, aaletter, cpdaa, outfile):
	# ukb fasta 2 csv
	fasta = uniprotkbFasta2csv(fasta)
	# making 'pos_ID' keys based on aaletter for all occurences of aaletter
	keyids = make_aapos_key(fasta, aaletter)
	# all AA keys as list
	aalist = list(keyids.pos_ID)
	# create list of detected keys
	cpdaa = pd.read_csv(cpdaa)
	cpdaaKey = list(cpdaa.pos_ID)
	ogcount = len(cpdaaKey)
	# overlap with positions from fasta sequences
	overlapAA = set(aalist).intersection(cpdaaKey)
	overlapcount  = len(overlapAA)
	notcount = ogcount -  overlapcount
	print("************************************")
	print("Total CpDAA positions from published study: ", ogcount)
	print("Total CpDAA positions overlapping provided fasta sequences: ",
	      overlapcount)
	print("Total CpDAA NOT found in fasta sequences: ", notcount)
	print("************************************")
	# new col in chemoproteomic results df for TRUE/FALSE if pos_ID was found
	#  in provided fasta file
	newls = list(overlapAA)
	finaldf = addcolumnconditional(newls, cpdaa, 'pos_ID',
	                                       'posID_in_checkedFASTA')
	print(finaldf.head(3))
	finaldf.to_csv(outfile, index=False)

if __name__ == '__main__':
	CpDAA_position_check(fasta, aaletter, cpdaa, outfile)