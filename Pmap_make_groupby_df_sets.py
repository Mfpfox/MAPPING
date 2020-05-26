

def setme(s):
    ls = list(s)
    finalset = set(ls)
    return finalset

def groupby_df(df):
    # columns based on set of IDs or sequences from all releases
    ENSP_df = df.groupby('xref')['ENSP'].apply(setme)
    ENSPv_df = df.groupby('xref')['ENSPv'].apply(setme)
    ENST_df = df.groupby('xref')['ENST'].apply(setme)
    ENSTv_df = df.groupby('xref')['ENSTv'].apply(setme)
    ENSG_df = df.groupby('xref')['ENSG'].apply(setme)
    ENSGv_df = df.groupby('xref')['ENSGv'].apply(setme)
    stableID_key_df = df.groupby('xref')['stableID_key'].apply(setme)
    proSequence_df = df.groupby('xref')['proSequence'].apply(setme)
    # concat
    set_columns = pd.concat([ENSP_df, ENSPv_df, ENST_df, ENSTv_df, ENSG_df, ENSGv_df, stableID_key_df, proSequence_df], axis=1)
    # counting length of sets aka # unique values in each column
    set_columns['count_ENSP'] = set_columns['ENSP'].apply(lambda x: len(x))
    set_columns['count_ENSPv'] = set_columns['ENSPv'].apply(lambda x: len(x))
    set_columns['count_ENST'] = set_columns['ENST'].apply(lambda x: len(x))
    set_columns['count_ENSTv'] = set_columns['ENSTv'].apply(lambda x: len(x))
    set_columns['count_ENSG'] = set_columns['ENSG'].apply(lambda x: len(x))
    set_columns['count_ENSGv'] = set_columns['ENSGv'].apply(lambda x: len(x))
    set_columns['count_stableID_key'] = set_columns['stableID_key'].apply(lambda x: len(x))
    set_columns['count_proSequence'] = set_columns['proSequence'].apply(lambda x: len(x))
    # columns based on list of things NOT SET, for SD calculation
    Length_df = df.groupby('xref')['Length'].apply(list)
    prov_df = df.groupby('xref')['pro_ver'].apply(list)
    txv_df = df.groupby('xref')['tx_ver'].apply(list)
    genv_df = df.groupby('xref')['gen_ver'].apply(list)
    # concat
    list_columns = pd.concat([Length_df, prov_df, txv_df, genv_df], axis=1)
    # creating column for number of ENSP ids
    # CANT DO VARIANCE WITH ONLY 1 VALUE
    # list_columns['stdev_length'] = list_columns['Length'].apply(lambda x: stdev(x))
    # list_columns['stdev_prov'] = list_columns['pro_ver'].apply(lambda x: stdev(x))
    # list_columns['stdev_txv'] = list_columns['tx_ver'].apply(lambda x: stdev(x))
    # list_columns['stdev_genv'] = list_columns['gen_ver'].apply(lambda x: stdev(x))
    # merging set df and list df 
    final = pd.concat([set_columns, list_columns], axis=1)
    print("final df shape: ", final.shape)
    print("preview: ")
    print(final.head(1))
    return final

def main(): 
    v85 = pd.read_csv("v85_fasta_merge_xref_3953IDs_10183.csv")
    v92 = pd.read_csv("v92_fasta_merge_xref_3953IDs_10395.csv")
    v94 = pd.read_csv("v94_fasta_merge_xref_3953IDs_10612.csv")
    v96 = pd.read_csv("v96_fasta_merge_xref_3953IDs_10663.csv")
    v97 = pd.read_csv("v97_fasta_merge_xref_3953IDs_10564.csv")
    v85final = groupby_df(v85)
    v92final = groupby_df(v92)
    v94final = groupby_df(v94)
    v96final = groupby_df(v96)
    v97final = groupby_df(v97)



