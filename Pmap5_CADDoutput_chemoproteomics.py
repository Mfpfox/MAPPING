#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Pmap_CADDoutput_chemoproteomics.py

import pandas as pd
import os
from ast import literal_eval

"""
# before running mapChemoproteomic filtered files with code below
def filterCADDoutput():
    inlist = ["chr10_CADD19_to_dbNSFP.csv","chr11_CADD19_to_dbNSFP.csv","chr12_CADD19_to_dbNSFP.csv","chr13_CADD19_to_dbNSFP.csv","chr14_CADD19_to_dbNSFP.csv","chr15_CADD19_to_dbNSFP.csv","chr16_CADD19_to_dbNSFP.csv","chr17_CADD19_to_dbNSFP.csv","chr18_CADD19_to_dbNSFP.csv","chr19_CADD19_to_dbNSFP.csv","chr1_CADD19_to_dbNSFP.csv","chr20_CADD19_to_dbNSFP.csv","chr21_CADD19_to_dbNSFP.csv","chr22_CADD19_to_dbNSFP.csv","chr2_CADD19_to_dbNSFP.csv","chr3_CADD19_to_dbNSFP.csv","chr4_CADD19_to_dbNSFP.csv","chr5_CADD19_to_dbNSFP.csv","chr6_CADD19_to_dbNSFP.csv","chr7_CADD19_to_dbNSFP.csv","chr8_CADD19_to_dbNSFP.csv","chr9_CADD19_to_dbNSFP.csv","chrX_CADD19_to_dbNSFP.csv","chrY_CADD19_to_dbNSFP.csv"]
    chrorder = ["chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr1","chr20","chr21","chr22","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"]
    for ii, ll in zip(inlist, chrorder):
        infile = pd.read_csv(ii)
        chrlabel = ll
        # parse out STOP CODON GAIN
        nostop = infile[infile['aaalt'] != 'X'].copy()
        nostop.reset_index(drop=True, inplace=True)
        # new column from CADD38 - CADD19 
        nostop['CADDdiff_38minus19'] = nostop['CADD_phred_hg38'] - nostop['CADD_phred_hg19'] 
        Foutname = chrlabel + "_CADD38_CADD19_noTer.csv"
        nostop.to_csv(Foutname, index=False)
        print("done with: ", Foutname)
        print(nostop.shape)
filterCADDoutput()
"""
def mapChemoproteomic(df):
    cysR = pd.read_csv('Cys_reactive_LMHset_1492.csv')
    cysP = pd.read_csv('Cys_panTar_5991.csv')
    lysR = pd.read_csv('Lys_reactive_LMHset_4504.csv')
    lysP = pd.read_csv('Lys_panTar_8220.csv')

    # cys dictionaries
    cysRdic = dict(zip(cysR.pos_ID, cysR.reactivity))
    cysValdic = dict(zip(cysR.pos_ID, cysR.Threshold))
    cysPdic = dict(zip(cysP.pos_ID, cysP.labelType))

    # cys map dict
    df['Cys_reactivity'] = df.pos_ID
    df['Cys_react_threshold'] = df.pos_ID
    df['Cys_target_label'] = df.pos_ID

    df.Cys_reactivity = df.Cys_reactivity.map(cysRdic)
    df.Cys_react_threshold = df.Cys_react_threshold.map(cysValdic)
    df.Cys_target_label = df.Cys_target_label.map(cysPdic)

    # lysine dictionaries
    lysRdic = dict(zip(lysR.pos_ID, lysR.reactivity))
    lysValdic = dict(zip(lysR.pos_ID, lysR.Threshold))
    lysPdic = dict(zip(lysP.pos_ID, lysP.labelType))

    # lys map dict
    df['Lys_reactivity'] = df.pos_ID
    df['Lys_react_threshold'] = df.pos_ID
    df['Lys_target_label'] = df.pos_ID

    df.Lys_reactivity = df.Lys_reactivity.map(lysRdic)
    df.Lys_react_threshold = df.Lys_react_threshold.map(lysValdic)
    df.Lys_target_label = df.Lys_target_label.map(lysPdic)
    return df


def main():
    inlist = ["chr10_CADD38_CADD19_noTer.csv","chr11_CADD38_CADD19_noTer.csv","chr12_CADD38_CADD19_noTer.csv","chr13_CADD38_CADD19_noTer.csv","chr14_CADD38_CADD19_noTer.csv","chr15_CADD38_CADD19_noTer.csv","chr16_CADD38_CADD19_noTer.csv","chr17_CADD38_CADD19_noTer.csv","chr18_CADD38_CADD19_noTer.csv","chr19_CADD38_CADD19_noTer.csv","chr1_CADD38_CADD19_noTer.csv","chr20_CADD38_CADD19_noTer.csv","chr21_CADD38_CADD19_noTer.csv","chr22_CADD38_CADD19_noTer.csv","chr2_CADD38_CADD19_noTer.csv","chr3_CADD38_CADD19_noTer.csv","chr4_CADD38_CADD19_noTer.csv","chr5_CADD38_CADD19_noTer.csv","chr6_CADD38_CADD19_noTer.csv","chr7_CADD38_CADD19_noTer.csv","chr8_CADD38_CADD19_noTer.csv","chr9_CADD38_CADD19_noTer.csv""chrX_CADD38_CADD19_noTer.csv","chrY_CADD38_CADD19_noTer.csv"]


    chrorder = ["chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr1","chr20","chr21","chr22","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"]

    for ii, ll in zip(inlist, chrorder):
        infile = pd.read_csv(ii)
        chrlabel = ll
        # TODO format matched_target
        infile.matched_target =\
            infile.matched_target.astype(str).str.replace('\[|\]|\'', '')

        # TODO add a/a column from list of p.X333Y
        s = infile.HGVSp_VEP
        ssplit = [x.split(";")[0] for x in s]
        ssplit = [x.replace("p.", "") for x in ssplit]
        oaa = [x[0:3] for x in ssplit]
        naa = [x[-3:] for x in ssplit]
        infile['lost_amino'] = oaa
        infile['gained_amino'] = naa
        infile['Amino_acids'] = infile['lost_amino'] +\
            '/' + infile['gained_amino']
        
        # parse out the target from non target positions
        falseTarget = infile[infile['matched_target'] == 'False'].copy()
        falseTarget.reset_index(drop=True, inplace=True)
        
        # creating pos_id for false target ck rows
        falseTarget['pos_ID_falseCKtarget'] = falseTarget['matched_UKBID'] + '_' + falseTarget['aaref'] + falseTarget['matched_aapos'].astype(str)
        # saving false target files
        Foutname = "CADDmapped_falsetarget_" + chrlabel + ".csv"
        falseTarget.to_csv(Foutname, index=False)
        print("done with: ", Foutname)
        # parse out the true target positions and add chemoproteomic details
        trueTarget = infile[infile['matched_target'] != 'False'].copy()
        trueTarget.reset_index(drop=True, inplace=True)

        trueTarget['pos_ID'] = trueTarget['matched_UKBID'] + '_' + \
            trueTarget['matched_target']
        chemopro = mapChemoproteomic(trueTarget)
        # saving trueTarget with chemoproteomics details
        Toutname = "CADDmapped_truetargets_" + chrlabel + ".csv"
        chemopro.to_csv(Toutname, index=False)
        print("done with: ", Toutname)

main()