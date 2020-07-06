# Welcome to MAP project repo 
This project is a collaboration between the [Backus lab](https://www.backuslab.com/) and [Arboleda lab](https://www.arboledalab.org/) at UCLA. Our paper, **From Chemoproteomic-Detected Amino Acids to Genomic Coordinates: Insights into Precise Multi-omic Data Integration**, can be found on [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.07.03.186007v1). 

#### Code is grouped by paper sections (1-5)

**[1]**

- Pmap_chemoproteomic_2012.py
- Pmap_filter_df.py
- Pmap_found_notfound_IDcheck.py
- processUniprotFasta.py

**[2]**

- Pmap_isoform_counting.py
- Pmap_filteron_posID.py
- Pmap_mismap_score_UKBvsUKB.py
- Pmap_protein_sequence_comparison.py

**[3]**

- Pmap_calc_multimapping.py
- Pmap_ensembl_fasta.py
- Pmap_make_Ensembl_multimap_df.py
- Pmap_MISMAP_score_funx.py
- Pmap_mismap_score_UKBvsENSP.py
- Pmap_multimapping_fig.py

**[4]**

- Pmap_accFilter_description.py
- Pmap_CADD_max_mean.py
- Pmap_correlation_plots.py
- Pmap_dbNSFP_CADD_merge_DECT19.py
- Pmap_dbNSFP_parselabeledID.py
- Pmap_missense_annotations_QC.py
- Pmap_pos_overlap_CADD.py
- Pmap1_parseID_correction.py
- Pmap2_dbNSFP_CysLysfilter.py
- Pmap3_findtargetmatch_dbNSFP.py
- Pmap4_dbNSFP_mapped_CADD19.py
- Pmap5_CADDoutput_chemoproteomics.py
- oddsRatio_detected_vs_notdet.R

**[5]**

- boxplot_reactivity.R
- linePlot_G6PD.R
- pymol_CADD_coloring.py

**[Mapping module]**

- maplib.py

---

python packages used: 
```python
import os
import sys
import numpy as np
import pandas as pd
import csv
from Bio import SeqIO
from ast import literal_eval
import difflib
from statistics import mean
import textdistance
```

R packages used: 
```R
library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(plotly)
```


