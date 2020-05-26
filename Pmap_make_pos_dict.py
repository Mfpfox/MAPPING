#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on
for Project Pmap, task A
@author: mfpfox
data: 11.6.19
creates dictionary of UKB ID: all cys and lys labeled positions
used to check 11_2012 sequences for CK chemoproteomic data positions
"""

"""

## create dictionary of {everLabeled entry : [list of CYS and LYS positions labeled]}
- GOAL: need to check if cys lys positions are matched in 11_2012 sequences from Cravatt lab

"""
import pandas as pd

# funx that turns list of labeled positions into sorted dictionary
# each key is index+1 of cys or lys that was labeled
# each value is chr of C or K
def get_pos_dictionary(pdcol):
    list_of_dicts = []
    for row in pdcol.tolist():
        listofkey=[]
        listofvalues=[]
        for i in row:
            v = i[0]
            listofvalues.append(v)
            k = i[1:]
            listofkey.append(k)
        di = dict(zip(listofkey, listofvalues))
        intdict = {int(oldkey): val for oldkey, val in di.items()}
        sdi = dict(sorted(intdict.items(), reverse=False))
        list_of_dicts.append(sdi)
    return list_of_dicts
def make_pos_dict():
    # import files
    cysever = "Cys2012_ever_labeled_6404.csv"
    lysever = "Lys2012_ever_labeled_9213.csv"
    cysever = pd.read_csv(cysever)
    lysever = pd.read_csv(lysever)
    # drop columns that i dont need
    cys = cysever.filter(['entry', 'pos'], axis=1)
    lys = lysever.filter(['entry', 'pos'], axis=1)
    # appending together at bottom, total shape is 
    appendedcyslys = cys.append(lys, ignore_index=True)
    print("appendedcyslys shape: ", appendedcyslys.shape)
    # creating list of positions labeled for each entry
    groupedCYSLYS = appendedcyslys.groupby('entry')['pos'].apply(list)
    groupedCYSLYS2 = pd.DataFrame(groupedCYSLYS, index=None)
    groupedCYSLYS2.reset_index(inplace=True)
    # creating column that is # of times entry has been labeled
    groupedCYSLYS2['labeled_pos_count'] = groupedCYSLYS2['pos'].str.len()
    # sorting df by decending # of times entry has been labeled
    groupedCYSLYS2_sort = groupedCYSLYS2.sort_values('labeled_pos_count',ascending=False)
    groupedCYSLYS2_sort.reset_index(inplace=True, drop=True)
    # calling funx by passing in df series (pos column)
    list_of_dicts = get_pos_dictionary(groupedCYSLYS2_sort['pos'])
    # turning list of dictionaries into a pandas column
    dictSeries = pd.Series(list_of_dicts)
    # adding dict col to original dataframe
    groupedCYSLYS2_sort['pos_dict'] = dictSeries
    # saving new file with positions of labeled cys and lys in dictionary form
    groupedCYSLYS2_sort.to_csv("dictOfPositions_allCysLys_.csv")
    return groupedCYSLYS2