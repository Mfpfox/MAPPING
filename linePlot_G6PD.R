# Amino acid position plot in figure 5

### notes: 
# - cadd scores for GRCh38 come from dbNSFP accession fltered files
# - cadd version v1.4 both models
# - dbnsfp version 4.0a
# - only missense CADD scores (filtered out nonsense)

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(plotly)
library(tidyr)

df <- read_csv('AAposition_plots/aapos_g6pd_350_cadd.csv')

# na values are currently coded as '-'
df <- na_if(df, '-')

names(df)<-str_replace_all(names(df), c(" " = "_")) 

df$Position <- as.integer(df$Position)

### code for figure 5 line plot ###

lineaa <- df %>%
  ggplot(aes(x=Position, y=max_CADD38)) + 
  scale_x_continuous(expand=c(0,0), limits=c(0, 350), breaks=seq(0, 350, 50))  + scale_y_continuous(expand=c(0,0), limits=c(0,45),breaks=seq(0, 45, 5)) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "none") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +  
  ylab("CADD38 PHRED scores") +  xlab("G6PD amino acid positions") + 
  theme(plot.title = element_text(hjust=0.5), axis.title.x = element_text(size=13, color="black", margin=margin(t=15, b=5)),axis.title.y = element_text(size=11, color="black", margin=margin(t=0, r=15, b=0)), axis.text=element_text(size=11, color="black",margin=margin(t=15, b=5))) + 
  geom_line(data = df, aes(x = Position, y = max_CADD38), size=.6,color = "black") +
  geom_vline(xintercept=c(205,171,97,89,82,47), color="blue", size=.2, alpha=0.9) +  
  geom_vline(xintercept=c(294,158,13), color="red", size=.2, alpha=0.9) +
  geom_hline(yintercept=25, linetype="dashed", color='grey',size=1)

ggsave("G6PDline_350aa_maxCADD38.pdf", width = 6, height=3, dpi= 300)

