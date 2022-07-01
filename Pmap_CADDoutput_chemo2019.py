#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Pmap_CADDoutput_chemo2019.py

import pandas as pd
import os
from ast import literal_eval

def mapChemoproteomic_Heta(df):
    # UPDATE for Heta's data:
    # mapChemoproteomic():  only cys isoTOP annotaiton
    # main(): no parsing into True False files, matched target col only remove '',
    # must have col pos_ID, reactivity, and Threshold
    cysR = pd.read_csv("preprocessed_HetaCYS2019_4112CpDC_2444UKBIDs.csv")
    # cys dictionaries
    cysRdic = dict(zip(cysR.pos_ID, cysR.reactivity))
    cysValdic = dict(zip(cysR.pos_ID, cysR.Threshold))
    # cys map dict
    df['Cys_reactivity'] = df.pos_ID
    df['Cys_react_threshold'] = df.pos_ID
    # map dict values
    df.Cys_reactivity = df.Cys_reactivity.map(cysRdic)
    df.Cys_react_threshold = df.Cys_react_threshold.map(cysValdic)
    return df


def main():
    # new filtered files with no stop and delta cadd38 -19
    inlist = ["chr10_CADD38_CADD37_noTer.csv","chr11_CADD38_CADD37_noTer.csv","chr12_CADD38_CADD37_noTer.csv",
              "chr13_CADD38_CADD37_noTer.csv",
              "chr14_CADD38_CADD37_noTer.csv","chr15_CADD38_CADD37_noTer.csv",
              "chr16_CADD38_CADD37_noTer.csv","chr17_CADD38_CADD37_noTer.csv",
              "chr18_CADD38_CADD37_noTer.csv","chr19_CADD38_CADD37_noTer.csv",
              "chr1_CADD38_CADD37_noTer.csv","chr20_CADD38_CADD37_noTer.csv",
              "chr21_CADD38_CADD37_noTer.csv","chr22_CADD38_CADD37_noTer.csv",
              "chr2_CADD38_CADD37_noTer.csv","chr3_CADD38_CADD37_noTer.csv",
              "chr4_CADD38_CADD37_noTer.csv","chr5_CADD38_CADD37_noTer.csv",
              "chr6_CADD38_CADD37_noTer.csv","chr7_CADD38_CADD37_noTer.csv",
              "chr8_CADD38_CADD37_noTer.csv","chr9_CADD38_CADD37_noTer.csv",
              "chrX_CADD38_CADD37_noTer.csv","chrY_CADD38_CADD37_noTer.csv"]    
    chrorder = ["chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                "chr1","chr20","chr21","chr22","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
                "chr9","chrX","chrY"]
    for ii, ll in zip(inlist, chrorder):
        infile = pd.read_csv(ii)
        chrlabel = ll
        # TODO format matched_target by removing quotes
        infile.matched_target = infile.matched_target.astype(str).str.replace('\'', '')
        # map key for chemoprot annotations
        infile['pos_ID'] = infile['matched_UKBID'] + '_' + infile['aaref'] + infile['matched_aapos'].astype(str)
        print(infile.head(1))
        
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
        # func call
        chemopro = mapChemoproteomic_Heta(infile)
        # saving scored data w/ chemoproteomics details
        Toutname = chrlabel + "_CADDmapped_targetsAnnotated" + ".csv"
        chemopro.to_csv(Toutname, index=False)
        print("done with: ", Toutname)
        print(chemopro.shape)
        checkColumnValues(chemopro, 'aaref')
        print(" ")
main()