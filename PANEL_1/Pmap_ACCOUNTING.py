
from collections import Counter
import pandas as pd
import sys
import os
# CONFIRMED ALL pos_ID unique, no redundancy within each subset, and updated MAPPING_tables.xlsx with new counts

sys.path.append("/Users/mariapalafox/Desktop/Toolbox")
from all_funx import *
os.chdir("/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/CYS_LYS_2012_to2018UKBpositions/REFERENCE_2018_mappedCK/SUBSETS/")

def addcolumnconditionalDrop(mapList, df, dfcol, newcol): 
    mendel = []
    for g in df[dfcol]:
        if g in mapList: 
            mendel.append("True")
        else:
            mendel.append("False")
    df.loc[:, newcol] = mendel
    df.drop(df[df[newcol] == "False"].index, inplace = True)
    df.reset_index(inplace=True, drop=True)
    print(df.shape)
    return df

# using ref ID list to filter out subset/ensembl files for row with matching ID
def drop_IDs():
    IDdf = pd.read_csv("UKBIDs_2018_3963.csv")
    idlist = list(IDdf.ID)
    print("total IDs to be searched: ", len(idlist))
    # list of files to be cleaned
    oldfiles = ["Lys2018_reactive_4425.csv","Lys2018_pantarget_8115.csv","LYS2018_everlabeled_gene_valueCounts_2624.csv","Lys2018_ever_labeled_9046.csv","CYSLYS2018_everlabeled_gene_valueCounts_3964.csv","Cys2018_reactive_1431.csv","Cys2018_pantarget_5962.csv","CYS2018_everlabeled_gene_valueCounts_2915.csv","Cys2018_ever_labeled_6234.csv","ALLReact2018_grouped_entry_threshold_2751.csv","AllReact2018_5856.csv"]

    newfiles = ["Lys2018_reactive_4425.csv","Lys2018_pantarget_8115.csv","LYS2018_everlabeled_gene_valueCounts_2624.csv","Lys2018_ever_labeled_9046.csv","CYSLYS2018_everlabeled_gene_valueCounts_3964.csv","Cys2018_reactive_1431.csv","Cys2018_pantarget_5962.csv","CYS2018_everlabeled_gene_valueCounts_2915.csv","Cys2018_ever_labeled_6234.csv","ALLReact2018_grouped_entry_threshold_2751.csv","AllReact2018_5856.csv"]

    for i, d in enumerate(oldfiles):
        df = pd.read_csv(d)
        dropUnnamed(df)
        print("oldfile up: ", d)
        df2 = addcolumnconditionalDrop(idlist, df, 'ID', 'final_refID')
        outname = newfiles[i]
        df2.to_csv(outname, index=False)


# ===ACCOUNTING UPDATE with 2018 protein sequences===

def colvalues(df, col):
    countdf = Counter(list(df[col]))
    countdic = dict(countdf)
    print(countdic)

def update_accounting():
    # UPDATED EVER LABELED FILES MANUALLY BY CHECKING NEW LINE COUNTS
    reactset = ["Lys2018_reactive_4425.csv","Cys2018_reactive_1431.csv","AllReact2018_5856.csv"]
    panset = ["Lys2018_pantarget_8111.csv","Cys2018_pantarget_5962.csv"]

    for i, d in enumerate(reactset):
        print(d)
        df = pd.read_csv(d)
        colvalues(df, 'Threshold')
        print()
        if d == "Lys2018_reactive_4425.csv":
            klow = df[df["Threshold"] == "Low"].describe()
            klow.columns = ['K_LOW']
            kmed = df[df["Threshold"] == "Medium"].describe()
            kmed.columns = ['K_MED']
            khigh=df[df["Threshold"] == "High"].describe()
            khigh.columns=['K_HIGH']
            klm = klow.join(kmed)
            klmh = klm.join(khigh)
            print(klmh)
            klmh.to_csv("KLMH_description_table_2018.csv")

        if d == "Cys2018_reactive_1431.csv":
            clow=df[df["Threshold"] == "Low"].describe()
            clow.columns=['C_LOW']
            cmed = df[df["Threshold"] == "Medium"].describe()
            cmed.columns=['C_MED']
            chigh=df[df["Threshold"] == "High"].describe()
            chigh.columns=['C_HIGH']
            CLM = clow.join(cmed)
            CLMH = CLM.join(chigh)
            print(CLMH)
            CLMH.to_csv("CLMH_description_table_2018.csv")

    for i, p in enumerate(panset):
        dfp = pd.read_csv(p)
        print(p)
        checkColumnValues(dfp, "labelType")



# =====DOUBLE CHECK ALL POSITIONS IN EVER LABELED FILES ARE UNIQUE
def count_unique_IDs(df, col):
    idlen = df[col]
    idset = set(df[col])
    print("ID column len: ", len(idlen))
    print("ID column set: ", len(idset))

def check_everlabeled():
    files = ["Lys2018_ever_labeled_9042.csv","Lys2018_reactive_4425.csv","Lys2018_pantarget_8111.csv","Cys2018_reactive_1431.csv","Cys2018_pantarget_5962.csv","Cys2018_ever_labeled_6234.csv","AllReact2018_5856.csv"]
    for d in files:
        print(d)
        df = pd.read_csv(d)
        print("shape: ", df.shape)
        count_unique_IDs(df, 'pos_ID')
        print()
# CONFIRMED ALL UNIQUE
def check_everlabeled_ID():
    files = ["AllPanTarget_CK_14073.csv","Lys2018_ever_labeled_9042.csv","Lys2018_reactive_4425.csv","Lys2018_pantarget_8111.csv","Cys2018_reactive_1431.csv","Cys2018_pantarget_5962.csv","Cys2018_ever_labeled_6234.csv","AllReact2018_5856.csv"]
    for d in files:
        print(d)
        df = pd.read_csv(d)
        print("shape: ", df.shape)
        count_unique_IDs(df, 'ID')
        print()



