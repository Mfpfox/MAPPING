def setme(s):
    ls = list(s)
    finalset = set(ls)
    return finalset

# infiles
v85 = pd.read_csv("v85_fasta_merge_xref_3955IDs_10188.csv")
v92 = pd.read_csv("v92_fasta_merge_xref_3955IDs_10399.csv")
v94 = pd.read_csv("v94_fasta_merge_xref_3955IDs_10616.csv")
v96 = pd.read_csv("v96_fasta_merge_xref_3955IDs_10667.csv")
v97 = pd.read_csv("v97_fasta_merge_xref_3955IDs_10568.csv")


# concat all releases together, stacked on top
ax0_3963 = pd.concat([v85, v92, v94, v96, v97])

'', '', '', '', '', 'ENSG', 'Length',
       '', '', 'xref'

# apply funx takes series
ENSPv_all_3963 = all_3963.groupby('xref')['ENSPv'].apply(setme)
ENSP_all_3963 = all_3963.groupby('xref')['ENSP'].apply(setme)
ENSTv_all_3963 = all_3963.groupby('xref')['ENSTv'].apply(setme)
ENST_all_3963 = all_3963.groupby('xref')['ENST'].apply(setme)
ENSGv_all_3963 = all_3963.groupby('xref')['ENSGv'].apply(setme)
ENSG_all_3963 = all_3963.groupby('xref')['ENSG'].apply(setme)
Length_all_3963 = all_3963.groupby('xref')['Length'].apply(setme)
stableID_key_all_3963 = all_3963.groupby('xref')['stableID_key'].apply(setme)
proSequence_all_3963 = all_3963.groupby('xref')['proSequence'].apply(setme)

# create df
panelA = pd.concat([ENSPv_all_3963, ENSP_all_3963, ENSTv_all_3963, ENST_all_3963, ENSGv_all_3963, ENSG_all_3963, Length_all_3963, stableID_key_all_3963, proSequence_all_3963], axis=1)