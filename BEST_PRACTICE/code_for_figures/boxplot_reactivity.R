# PANEL 5 boxplots of CADD38 max codon scores of Low/Medium/High reactivity cysteine and lysine residues

library(tidyverse)
library(ggplot2)
library(readr)
library(scales)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(tidyr)

det_scores <- read_csv("/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/CADDmapped/ALL_CONSEQUENCES/Rfig_SCORES_max_mean_detected_14925.csv") %>%
  mutate(Cys_react_threshold = factor(Cys_react_threshold, levels = c(
    "Low",
    "Medium",
    "High"))) %>%
  mutate(Lys_react_threshold = factor(Lys_react_threshold, levels = c(
    "Low",
    "Medium",
    "High"))) %>%
  mutate(Cys_target_label = factor(Cys_target_label, levels = c(
    "panReactive",
    "Target"))) %>%
  mutate(Lys_target_label = factor(Lys_target_label, levels = c(
    "panReactive",
    "Target"))) 

names(det_scores)<-str_replace_all(names(det_scores), c("-" = "","_" = "." ))

cLMH <- det_scores %>%  
  filter(!is.na(Cys.reactivity))
    # 1401 LMH cysteine total

kLMH <- det_scores %>% 
   filter(!is.na(Lys.reactivity))
    # 4363 LMH lysine total 


### CYSTEINE MAX CADD38 boxplot ###
my_comparisons <- list(c("Low", "Medium"), c("Medium", "High"), c("Low", "High"))
mystats2 <- compare_means(CADD38.phred.max ~ Cys.react.threshold, data=cLMH, method="wilcox.test", p.adjust.method ="BH", paired=FALSE, comparisons = my_comparisons) 
mystats2 <- mystats2 %>% mutate(y.position = c(40, 50, 45))
cstat <- ggboxplot(cLMH, x = "Cys.react.threshold", y = "CADD38.phred.max",fill = "Cys.react.threshold", alpha=0.3, notch=TRUE) + 
  stat_compare_means(label.y = 55) +    # Add global p-value  
  stat_pvalue_manual(mystats2, label ="p.adj") +
  scale_fill_manual(values=c("#00C1AA", "#529EFF", "gold")) +  
  theme_minimal() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x="Cysteine reactivity thresholds", y= "CADD hg38 max phred scores / codon") +  
  theme(legend.position='none') + 
  theme(plot.title = element_text(hjust=0.5), axis.title.x = element_text(size=12, color="black", margin=margin(t=15, b=5)),axis.title.y = element_text(size=12, color="black", margin=margin(t=0, r=15, b=0)), axis.text=element_text(size=10, color="black", margin=margin(t=20, b=20)))  
ggsave("boxplot_cysteine_LMH_pairwise_padjBH_global_kruskal_1401aa.pdf", width = 6, height=6, dpi= 300) 


### LYSINE MAX CADD38 boxplot stats ###
my_comparisons <- list(c("Low", "Medium"), c("Medium", "High"), c("Low", "High"))
mystats3 <- compare_means(CADD38.phred.max ~ Lys.react.threshold, data=kLMH, method="wilcox.test", 
                          p.adjust.method ="BH", paired=FALSE, comparisons = my_comparisons) 
mystats3 <- mystats3 %>% mutate(y.position = c(40, 50, 45))
kstat <- ggboxplot(kLMH, x = "Lys.react.threshold", y = "CADD38.phred.max",fill = "Lys.react.threshold", 
                   alpha=0.3, notch=TRUE) + stat_compare_means(label.y = 60) +    # Add global p-value  
  stat_pvalue_manual(mystats3, label ="p.adj") + scale_fill_manual(values=c("#00C1AA", "#529EFF", "gold")) +  
  theme_minimal() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), 
                          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x="Lysine reactivity thresholds", y= "CADD hg38 max phred scores / codon") +  
  theme(legend.position='none') + 
  theme(plot.title = element_text(hjust=0.5), axis.title.x = element_text(size=12, color="black", margin=margin(t=15, b=5)),axis.title.y = element_text(size=12, color="black", margin=margin(t=0, r=15, b=0)), axis.text=element_text(size=10, color="black", margin=margin(t=20, b=20)))  
ggsave("boxplot_lysine_LMH_pairwise_padjBH_global_kruskal_4363aa.pdf", width = 6, height=6, dpi= 300) 

