# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import pandas as pd
import csv
from statistics import mean 

def add_ids(df):
    # add pos_id col
    dfcol = df.columns
    if 'pos_id' not in dfcol:
        # input names of columns for ukb id and aa pos
        df['pos_id'] = df['matched_UKBID'] + '_' + df['matched_aapos'].astype(str)
    # create aaref_pos columns
    if 'aapos_id' not in dfcol:
        df['aapos_id'] = df['aaref'] + '_' + df['matched_aapos'].astype(str)
    return df


def add_maxmean_col(df, scorecol):
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

def main():
    os.chdir("")
    file = ''
    df = pd.read_csv(file)
    # add_ids(df)
    # df columns = aapos_id matched_UKBID   CADD38_phred
    add_maxmean_col(df)
main()
