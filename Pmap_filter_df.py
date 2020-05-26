# importing all_funx into env eliminates need to import packages
# script filters based on seperate file containing single column of IDs

import sys
import os
sys.path.append("/Users/mariapalafox/Desktop/Toolbox")
from all_funx import *
# see also filter_by_col.py

def addcolumnconditionalDrop(mapList, df, dfcol, newcol):
    mendel = []
    for g in df[dfcol]:
        if g in mapList:
            mendel.append("True")
        else:
            mendel.append("False")
    df.loc[:, newcol] = mendel
    checkColumnValues(df, newcol)
    print("dropping false")
    df.drop(df[df[newcol] == "False"].index, inplace=True)
    df.reset_index(inplace=True, drop=True)
    print("df shape post drop: ", df.shape)
    return df

def ID_list(df, col):
    idls = list(df[col])
    print("len of ID list: ", len(idls))
    return idls

def filter_df():
    # UPDATE df ids idls colname and colname in addcoldrop input
    df = pd.read_csv("Lys_ever_labeled_set_9162.csv")
    ids = pd.read_csv("only_3901_accessions.csv")
    idls = ID_list(ids, 'ID')
    # Cys_ever_labeled_set_6315.csv
    # Lys_ever_labeled_set_9162.csv
    df2 = addcolumnconditionalDrop(idls, df, 'entry', 'in3901')
    df2.to_csv("Lys_ever_labeled_set_3901_filtered.csv", index=False)

filter_df()

