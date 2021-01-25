## ChemoProteomic-Detected Amino Acids (CpDAAs) mapped to missense pathogenicity predictions and ClinVar pathogenic variants

**Visit the [CpDAA shiny app](http://mfpalafox.shinyapps.io/CpDAA) to visualize deleteriousness scores and ClinVar mutation overlap for cysteine and lysine residues in 4,526 chemoproteomic-detected proteins. Link to manuscript [here](https://www.biorxiv.org/content/10.1101/2020.07.03.186007v1).** 

---
---
---

### Protocol for mapping CpDAA:

#### [STEP 1] Pre-processing of published cysteine and lysine chemoproteomic datasets
*Combining quantitative chemoproteomics studies for meta-analysis is hindered by inter-study variability in experimental design and data normalization and quantification pipelines. Additionally, differences in available meta-data and file formats from chemoproteomic studies limits automated large-scale integration. Therefore, the legacy cysteine isoTOP-ABPP and the dataset from 2019 were kept separate for reactivity and missense pathogenicity analysis.*

❏    For available residue-level chemoproteomic result files, filter out peptides with multiple amino acids marked as modified (e.g. MA**C**AL**C**AEK)

❏    Filter detected residue data by applying rigorous peptide filtering, such as setting requirement for minimum # of times residue must be detected in replicate samples

❏    For reactivity profiling results, average ratio (R) values across samples in order to assign one R value per detected residue, `avg_reactivity()` in maplib.py

❏    Filter out residues with average R values < 0.5 

❏    Label all residues from reactivity profiling experiments as Low, Medium, or High based on user-defined R value range, `add_thresholds()` in maplib.py

❏    Label filtered chemoproteomic residue-level data by position ID key (e.g. ‘P11413_K170’) in order to account for all unique CpDAA positions in UniProtKB proteins, `UKBposID()` in maplib.py)

#### [STEP 2] Mapping residue-level chemoproteomic data to an annotation reference proteome 
*Instead of merging and reprocessing raw results from several chemoproteomic studies, we re-mapped all residue-level data to a version of the UniProtKB human proteome serving as our functional annotation reference in order to facilitate comparisons between annotated sets of detected cysteines and lysines. For some chemoproteomic studies in our analysis, the original search database FASTA file was not accessible online, and instead we used a version of UniProtKB canonical sequences released in the same year (release 2012_11) as the version used in the study. The python scripts below were used to **(A)** find the closest available UniProtKB release matching positions of CpDAA in chemoproteomic results file and **(B)** compare canonical sequences from two UniProtKB releases and note sequence differences.* 

**checkCpDAAinUKB.py**

❏    Input parameters: 

  ❏    UniProtKB FASTA filename
  
  ❏    single letter amino acid type (e.g. C for Cysteine)
  
  ❏    residue-level chemoproteomic results filename that includes column for CpDAA keys (added in STEP 1)
  
  ❏    .csv filename for output residue-level data with additional column (boolean) for CpDAA key search results
  
❏    Reads in FASTA using `uniprotkbFasta2csv()` in maplib.py

❏    Calls `make_aapos_key()` in maplib.py on protein sequence column to create key IDs with the  protein ID and positions of all amino acids matching the input parameter for type  (e.g. ‘P11413_C17’) 

❏    Reads in residue-level results, makes list from column of position ID keys

❏    Intersects position keys from FASTA with detected residue keys, adds boolean column to chemoproteomic results dataframe based on overlap

❏    Saves sequence checked CpDAA results as .csv file

**compareCanonicalSequences.py**

❏    Input parameters: 

  ❏    UniProtKB FASTA filename A
  
  ❏    UniProtKB FASTA filename B 
  
  ❏    .csv filename containing column of UniProtKB IDs to filter FASTA by
  
  ❏    .csv filename to save script output (modifies FASTA B by adding sequence comparison results)
  
❏    Reads in FASTA using `uniprotkbFasta2csv()` in maplib.py

❏    Reads .csv with IDs

❏    Selects first column of IDs from each data frame and creates set of IDs shared between all sources

❏    Filters FASTA data frames to only contain IDs and sequences found in shared ID set

❏    Creates dictionary of {UniProtKB ID : proSequence} based on FASTA-A, loops over rows of FASTA-B data frames to compare

❏    Protein length differences (B - A)

❏    Amino acid differences

❏    Positions of differing amino acids (1-index)

❏    Modifies the FASTA-B data frame with columns for comparison results

#### [STEP 3] Prep for dbNSFP mapping with UniProtKB BED file (from 2018 to match sequence versions)
*Goal of annotating UniProt stable IDs by GRCh38 chromosome coordinates using genomic coordinate mappings from [McGarvey et al. 2019](https://doi.org/10.1002/humu.23738).*

❏    Search for UniProt stable IDs in BED file, drop protein IDs not contained in BED (IDs lacking GRCh38 mapping)

❏    Drop protein IDs mapping to more than one chromosome

❏    Parse UniProt IDs based on mapped chromosome column, `makeCHRfiles()` in maplib.py

❏    Create posID key for all undetected equivalent amino acids in list of UniProt IDs, `makeAAposID()` in maplib.py

❏    Labels residue-level rows with ‘detected’ or ‘undetected’ flag, based on whether residue was identified in chemoproteomic studies

#### [STEP 4] Identify Ensembl protein IDs with identical sequences in order to map transcript-dependent annotations (see Table X with mapping terms) such as gene constraint metrics
*Identifies Ensembl IDs from a particular database release that are associated with sequences identical to a reference sequence from UniProtKB. Outputs a custom cross-reference key validated for equivalent protein sequences between UniProtKB and Ensembl releases. Includes columns for Ensembl IDs, UniProtKB IDs, HGNC info, protein length, genomic coordinates pulled from Ensembl FASTA (assembly, chr, cds start, cds stop).*

**identicalProSequenceChecker.py**

❏    Input parameters:

  ❏    Ensembl peptide FASTA
  
  ❏    Ensembl-UniProt cross-reference (xref) file from Ensembl
  
  ❏    Ensembl release number that matched FASTA and xref file
  
  ❏    UniProtKB FASTA
  
  ❏    UniProtKB database, options 'SP' for SWISSPROT or 'SPTREMBL' for SWISSPROT & TREMBL
  
❏    Filters out non-standard chromosome names

❏    Filters out Ensembl rows that do not have protein sequence identical to UniProtKB sequence reference

❏    Maps info columns from Ensembl FASTA with UniProtKB sequence info to create custom xref key

❏    Saves custom xref key as “{Ensembl release #}CheckedProSeq_ENSPtoUKB.csv”


#### [STEP 5] Map CpDAA to dbNSFP pathogenicity scores, CADD v1.4 GRCh37 scores, and ClinVar
*The following scripts are not modular (hardcoded file names) and were used to map chemoproteomic annotations from 3 published studies. The tasks performed by each script are described below.*

**Pmap1_parseID_correction.py**

❏    Filters dbNSFP with UKB ID chr dictionary, keep rows with UniProtKB stable ID match

❏    Adds columns 'matched_UKBID' 'matched_aapos' 'matched_index' from parsing dbNSFP row lists

**Pmap2_dbNSFP_CysLysfilter.py**

❏    Filters dbNSFP ID matched files for only AA reference column values that match AA type of CpDAA data (our study looks at both cysteines and lysines)

❏    Make position dictionary file of CpDAA data, `get_pos_dictionary()` in maplib.py

**Pmap3_findtargetmatched_dbNSFP.py**

❏    Makes GRCh37- and GRCh38-based genomic coordinate key IDs

❏    Drops variant rows with missing genomic coordinates for either genome reference assembly (ensure compatibility between both assemblies)

❏    Maps CpDAA annotations from chemoproteomic studies to dbNSFP data using position dictionary file made by `get_pos_dictionary()` in maplib.py

❏    Saves files with dbNSFP and CpDAA annotations

**Pmap4_dbNSFP_mapped_CADD37.py**

❏    Reads in CADD v1.4 model GRCh37 chromosome-parsed files 

❏    Maps CADD file columns to annotated files from `Pmap3_findtargetmatched_dbNSFP.py` using GRCh37-based coordinate key IDs

❏    Saves .csv of data frame with dbNSFP, CADD, and CpDAA annotations

❏    Remove missense consequence gain of stop (‘X’ amino acid)

❏    Add new column for difference between CADD models (GRCh38 scores - GRCh37 scores)

❏    Parse HGVS column and add formatted missense column (same format used in ClinVar file)

❏    QC 7 overlapping nonsynonymous SNVs per Cys and Lys codon (only non-stop gain missense)

❏    Drop proteins with missing cysteine or lysine residue-level data (final set of proteins has all cysteines and lysines accounted for) 

❏    Map ClinVar pathogenic & likely pathogenic variants to annotated data using GRCh37-based coordinate key IDs and confirm missense consequence and protein position matches between sources

❏    Save missense-level data as .csv (shiny app input file)

❏    Take max and mean score summaries for all detected and undetected cysteine and lysine posID keys, `addMaxMean()` in maplib.py, save residue-level data as .csv (shiny app input file)

