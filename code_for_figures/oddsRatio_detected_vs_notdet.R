# PANEL 4 MISSENSE C>W K>E ODDS RATIO of Detected vs Undetected 
# using CADD, DANN, FATHMM scores w/ thresholds

library(tidyverse)
library(ggplot2)
library(readr)
library(scales)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(tidyr)

### data import ###

all_scores <- read_csv("ALL_CONSEQUENCES/MERGE_COMBO_cadd38corrected_1327386.csv") %>%
  mutate(GROUP = factor(GROUP, levels = c(
    "detected",
    "notdetected"
  ))) # files has all missense for detected and not detected cysteine and lysines

names(all_scores)<-str_replace_all(names(all_scores), c("-" = "","_" = "." ))


### Filtering ###

posC <- all_scores %>% 
  filter(aaref == 'C') %>% 
  filter(aaalt == 'W') # all cysteine positions total: 40,107


posK <- all_scores %>% 
  filter(aaref == 'K') %>% 
  filter(aaalt == 'E') # all lysine positions total: 149,520


### Contingency tables based on thresholds ###

# CADD38 phred >= 25
posC$CADD25 <- ifelse(posC$CADD.phred.hg38 >= 25, "True", "False") 
posC$CADD25 <- factor(posC$CADD25, levels = c(
  "True",
  "False"
))

posK$CADD25 <- ifelse(posK$CADD.phred.hg38 >= 25, "True", "False") 
posK$CADD25 <- factor(posK$CADD25, levels = c(
  "True",
  "False"
))

# DANN >= 0.98
posC$DANN98 <- ifelse(posC$DANN.score >= 0.98, "True", "False") 
posC$DANN98 <- factor(posC$DANN98, levels = c(
  "True",
  "False"
))

posK$DANN98 <- ifelse(posK$DANN.score >= 0.98, "True", "False") 
posK$DANN98 <- factor(posK$DANN98, levels = c(
  "True",
  "False"
))

# FATHMM >= 0.95
posC$FATHMM95 <- ifelse(posC$fathmmMKL.coding.score >= 0.95, "True", "False") 
posC$FATHMM95 <- factor(posC$FATHMM95, levels = c(
  "True",
  "False"
))

posK$FATHMM95 <- ifelse(posK$fathmmMKL.coding.score >= 0.95, "True", "False") 
posK$FATHMM95 <- factor(posK$FATHMM95, levels = c(
  "True",
  "False"
))

# CYSTEINE
letters = c('CADD25', 'FATHMM95', 'DANN98')
i = 1
N = 3
rrows <- list()
for (i in i:N){
  low = letters[i]
  testX <- table(posC[[low]], posC$GROUP)
  matrixX <- matrix(testX, nrow=1)
  rrows[[i]] <- matrixX
}
final<- do.call(rbind, rrows)
finaldf <- as.data.frame(final)
write_csv(finaldf, "C_W_contingency_panel4.csv")

## LYSINE
letters = c('CADD25', 'FATHMM95', 'DANN98')
i = 1
N = 3
rrows <- list()
for (i in i:N){
  low = letters[i]
  testX <- table(posK[[low]], posK$GROUP)
  matrixX <- matrix(testX, nrow=1)
  rrows[[i]] <- matrixX
}
final<- do.call(rbind, rrows)
finaldf <- as.data.frame(final)
write_csv(finaldf, "K_E_contingency_panel4.csv")

# C>W analysis
analysis <- read.csv("C_W_contingency_panel4.csv", header=FALSE, sep=",") %>% 
  setNames(c("group","V1","V2","V3","V4")) %>% 
  nest(-group) %>% 
  mutate(matrix=map(data, ~matrix(unlist(.x), nrow=2))) %>% 
  mutate(fisher = map(matrix, ~fisher.test(.x))) %>% 
  mutate(stats = map(fisher, ~broom::glance(.x)))

stats1 <- analysis %>% 
  unnest(stats) %>% 
  select(group, p.value, odds=estimate, conf.low, conf.high)
write_csv(stats1, "C_W_fisher_panel4.csv")

# K>E analysis
analysis <- read.csv("K_E_contingency_panel4.csv", header=FALSE, sep=",") %>% 
  setNames(c("group", "V1","V2","V3","V4")) %>% 
  nest(-group) %>% 
  mutate(matrix=map(data, ~matrix(unlist(.x), nrow=2))) %>% 
  mutate(fisher = map(matrix, ~fisher.test(.x))) %>% 
  mutate(stats = map(fisher, ~broom::glance(.x)))

stats2 <- analysis %>% 
  unnest(stats) %>% 
  select(group, p.value, odds=estimate, conf.low, conf.high)

write_csv(stats2, "K_E_fisher_panel4.csv")


### FIGURE MAKING CODE ###

  # C>W fisher group	p.value
  #CADD25	3.40E-22 yes ***
  #FATHMM95	0.019891847 *
  #DANN98	6.69E-26 ***

  # K>E
  #group	p.value
  #CADD25	1.37E-78 ***
  #FATHMM95	1.97E-21 ***
  #DANN98	1.10E-28 ***

boxLabels = c("CADDhg38 >= 25", "FATHMMmkl >= 0.95", "DANN >= 0.98","CADDhg38 >= 25", "FATHMMmkl >= 0.95", "DANN >= 0.98")
y2 = length(boxLabels):1
outname <- "ODDS_CW_KE_DETvsNOT_3840_3scores.pdf"
df2 <- data.frame(
  boxOdds = c(0.762167374, 0.921329137, 0.689701922, 1.511460233, 1.240922511, 1.427656956),
  boxCILow = c(0.72123527, 0.859314503, 0.644416395, 1.446995783, 1.186290675, 1.337242924),
  boxCIHigh = c(0.8053926, 0.987277438, 0.738476628, 1.579000173, 1.298278156, 1.525348324))

p2 <- ggplot(df2, aes(x = boxOdds, y = y2)) + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .2, height = 0.5, color = "gray50") +
  geom_point(size = 2.5 , color = "black", shape=15) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) + scale_y_continuous(breaks = y2, labels = boxLabels) + 
  scale_x_continuous() +
  ylab("") +   xlab("Odds missense at CpDAA more deleterious than\nat not detected AA (n=3,840 proteins)") +
  xlim(0,2)  +
  theme(axis.title.x = element_text(size=8, color="black", margin=margin(t=15, b=15)), axis.text=element_text(size=8, color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2
ggsave(outname, width = 4.5, height=2, dpi= 300) 

# ---------------------------------------------------------

boxLabels = c("CADDhg38 >= 25", "FATHMMmkl >= 0.95", "DANN >= 0.98","CADDhg38 >= 25", "FATHMMmkl >= 0.95", "DANN >= 0.98")
y2 = length(boxLabels):1
outname <- "ODDS_CW_KI_DETvsNOT_3840_3scores.pdf"
df2 <- data.frame(
  boxOdds = c(0.762167374, 0.921329137, 0.689701922,1.801176997, 1.548448051, 1.749913143),
  boxCILow = c(0.72123527, 0.859314503, 0.644416395, 1.665568695, 1.438843127, 1.492835862),
  boxCIHigh = c(0.8053926, 0.987277438, 0.738476628,1.949192938,1.667232172,2.062922539))

p2 <- ggplot(df2, aes(x = boxOdds, y = y2)) + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .2, height = 0.5, color = "gray50") +
  geom_point(size = 2.5 , color = "black", shape=15) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) + scale_y_continuous(breaks = y2, labels = boxLabels) + 
  scale_x_continuous() +
  ylab("") +   xlab("Odds missense at CpDAA more deleterious than\nat not detected AA (n=3,840 proteins)") +
  xlim(0,2.5)  +
  theme(axis.title.x = element_text(size=8, color="black", margin=margin(t=15, b=15)), axis.text=element_text(size=8, color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2
ggsave(outname, width = 4.5, height=2, dpi= 300) 