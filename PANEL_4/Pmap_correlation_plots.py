# only looking at ML scores for detected cys and lys
# specific missense changes
# 

os.chdir("/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/CADDmapped/ALL_CONSEQUENCES/")
combo = pd.read_csv("MERGE_COMBO_cadd38corrected_1327386.csv")

# DROP NAN ROWS IF USING ALL SCORES IN CORRELATION PLOT
combo.dropna(inplace=True)
combo.drop(['Unnamed: 0'], axis=1,inplace=True)

"""
# missense counts before drop nan
   Amino_acids   Count
0      Lys/Asn  299040
1      Lys/Arg  149520
2      Lys/Thr  149520
3      Lys/Gln  149519
4      Lys/Glu  149519
5      Lys/Met   82554
6      Cys/Ser   80214
7      Lys/Ile   66966
8      Cys/Arg   40107
9      Cys/Tyr   40107
10     Cys/Gly   40107
11     Cys/Phe   40107
12     Cys/Trp   40106

# missense counts following drop
   Amino_acids   Count
0      Lys/Asn  203230
1      Lys/Thr  101812
2      Lys/Gln  101759
3      Lys/Arg  101721
4      Lys/Glu  101706
5      Cys/Ser   56677
6      Lys/Met   55754
7      Lys/Ile   45842
8      Cys/Phe   28351
9      Cys/Gly   28349
10     Cys/Arg   28341
11     Cys/Tyr   28322
12     Cys/Trp   28251

Cys AA : 4227
Lys AA : 5992
cys lost from dropnan 1830
lys lost from dropnan 2876
total missense : 
detected cys missense:  (29541, 30)
detected lys misssense:  (41850, 30)
total ck missense: 71391
"""


# function for making spearman heat maps
def make_corr_heat(df, savename):
    corr = df.corr(method='spearman')
    mask = np.zeros_like(corr)
    mask[np.triu_indices_from(mask)] = True
    with sn.axes_style("white"):
        f, ax = plt.subplots(figsize=(7, 5))
        ax =sn.heatmap(corr, annot=True, annot_kws={"size": 6},vmin=0, vmax=1, mask=mask, square=True, cmap="PuBu")
        plt.subplots_adjust(top=1, bottom=0.5)
        #ax.set_ylim(len(8)-0.5, -0.5)
        # fix for mpl bug that cuts off top/bottom of seaborn viz
        b, t = plt.ylim() # discover the values for bottom and top
        b += 0.5 # Add 0.5 to the bottom
        t -= 0.5 # Subtract 0.5 from the top
        plt.ylim(b, t) # update the ylim(bottom, top) values
    plt.savefig(savename, dpi=300, bbox_inches = "tight")
    plt.show()

def make_corr(df, missense, group):
    # 'Cys/Trp', 'detected'
    df = df[df['Amino_acids'] == missense].copy()
    df = df[df['GROUP'] == group].copy()
    print("total missense in df: ", df.shape)
    print("total specific AA count in group: ", len(df.pos_ID.unique()))
    # subset df
    df = df[['M-CAP_score','REVEL_score', 
       'MPC_score', 'PrimateAI_score', 'DANN_score', 
       'fathmm-MKL_coding_score', 'RawScore_hg19',
       'CADD38_raw', 'GROUP']].copy()
    # order by 37 group and 38 group
    colorder = ['RawScore_hg19', 'fathmm-MKL_coding_score',
            'M-CAP_score','REVEL_score', 'PrimateAI_score', 'MPC_score',  
            'DANN_score', 'CADD38_raw', 'GROUP']
    # setting new col order
    df = df[colorder]
    df.columns = ['CADD37', 'fathmmMKL',
            'MCAP','REVEL', 'PrimateAI', 'MPC',  
            'DANN', 'CADD38', 'GROUP']
    return df

# scores scores for specific missense groups
KRplot_det = make_corr(combo, 'Lys/Arg', 'detected')
KQplot_det = make_corr(combo, 'Lys/Gln', 'detected')
KEplot_det = make_corr(combo, 'Lys/Glu', 'detected')
KTplot_det = make_corr(combo, 'Lys/Thr', 'detected')
KNplot_det = make_corr(combo, 'Lys/Asn', 'detected')
KMplot_det = make_corr(combo, 'Lys/Met', 'detected')
KIplot_det = make_corr(combo, 'Lys/Ile', 'detected')

CSplot_det = make_corr(combo, 'Cys/Ser', 'detected')
CGplot_det = make_corr(combo, 'Cys/Gly', 'detected')
CRplot_det = make_corr(combo, 'Cys/Arg', 'detected')
CYplot_det = make_corr(combo, 'Cys/Tyr', 'detected')
CFplot_det = make_corr(combo, 'Cys/Phe', 'detected')
CWplot_det = make_corr(combo, 'Cys/Trp', 'detected')


dfs = [CWplot_det, KEplot_det, CWplot_not, KEplot_not]
dfnames = ['corr_CWdetected_4209aa.pdf', 'corr_KEdetected_5979aa.pdf', 'corr_CWnotdet_24042aa.pdf', 'corr_KEnotdet_95727aa.pdf']
for x,y in zip(dfs, dfnames):
    make_corr_heat(x, y)