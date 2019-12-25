# [1] run --tab 
./vep -i data/MAP_CK/BED_DETECTED_HGVSchange_104475.txt -o VEP_DETECTED_HGVSchange_104475.txt \
--cache --offline --tab --assembly GRCh37 --force_overwrite --port 3337 \
--gencode_basic --check_existing --symbol --canonical --tsl --appris --biotype --total_length \
--uniprot --protein --ccds --domains --xref_refseq --coding_only --pick \
--af_gnomad --af --max_af --gene_phenotype --numbers --regulatory --af_1kg \
--plugin Blosum62 \


# [2]
./vep -i data/MAP_CK/BED_DETECTED_HGVSchange_104475.txt -o VEP_DETECTED_HGVSchange_104475_v2.txt \
--cache --offline --tab --assembly GRCh37 --force_overwrite --port 3337 \
--gencode_basic --symbol --canonical --total_length \
--uniprot --protein --domains --coding_only --pick \
--plugin Blosum62 \
--plugin MPC,/Users/mariapalafox/.vep/Plugins/fordist_constraint_official_mpc_values_v2.txt.gz \
--plugin dbNSFP,data/dbNSFP4.0a_hg19/dbNSFP4.0a_hg19.gz,hg19_chr,hg19_pos\(1-based\),codon_degeneracy,codonpos,M-CAP_score,M-CAP_pred,REVEL_score,MutPred_score,MutPred_protID,MutPred_AAchange,MutPred_Top5features,GERP++_NR,GERP++_RS,PrimateAI_score,PrimateAI_pred,clinvar_id



-----------------------------------------------------------------
REMOVED: 
--fork 6\
q
| ./filter_vep -filter "Consequence is missense_variant" -i VEP_DETECTED_HGVSchange_104475.txt -o VEP_detected_filtered_missenseonly.txt

--fasta data/ref/ --hgvs --hgvsg 

--vcf

--plugin LoFtool,LoFtool_scores.txt \
--plugin ExACpLI,ExACpLI_values.txt 

(to generalize https://github.com/Ensembl/VEP_plugins/issues/108)