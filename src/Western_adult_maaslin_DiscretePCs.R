setwd("/Users/tobynbranck/Documents/Effect_size/analysis/")
library(ggpubr)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(permute)
library(vegan)
library(tidyverse)
library(lme4)
library(factoextra)
library(ade4)
library(FactoMineR)
library(psych)
library(stringr)
library(reshape2)
require(cowplot)

hmp2_maaslin2_sig_results_cd <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_DISCRETE_diseaseStrat_CD/all_results.tsv", header=T, check.names=T))
hmp2_maaslin2_sig_results_cd$metadata <- paste("HMP2_cd", hmp2_maaslin2_sig_results_cd$metadata, sep="_")

hmp2_maaslin2_sig_results_uc <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_DISCRETE_diseaseStrat_UC/all_results.tsv", header=T, check.names=T))
hmp2_maaslin2_sig_results_uc$metadata <- paste("HMP2_uc", hmp2_maaslin2_sig_results_uc$metadata, sep="_")

hmp2_maaslin2_sig_results_nonibd <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_DISCRETE_diseaseStrat_a_nonIBD/all_results.tsv", header=T, check.names=T))
hmp2_maaslin2_sig_results_nonibd$metadata <- paste("HMP2_nonIBD", hmp2_maaslin2_sig_results_nonibd$metadata, sep="_")

mlvs_maaslin2_sig_results <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_maaslin2_PCs/all_results.tsv", header=T, check.names=T))
mlvs_maaslin2_sig_results$metadata <- paste("MLVS", mlvs_maaslin2_sig_results$metadata, sep="_")

testt = rbind(hmp2_maaslin2_sig_results_uc,hmp2_maaslin2_sig_results_cd,hmp2_maaslin2_sig_results_nonibd,mlvs_maaslin2_sig_results)
testt = testt[- grep("antibiotic", testt$metadata),]
testt = testt[- grep("diagnosis2", testt$metadata),]
testt = testt[- grep("age", testt$metadata),]
testt = testt[- grep("WeekBin", testt$metadata),]
testt = testt[- grep("WeekNum", testt$metadata),]
testt = testt[- grep("SampNum", testt$metadata),]
testt$metadata = gsub("HMP2_cd_RC1","HMP2 CD dietary pattern 1",testt$metadata)
testt$metadata = gsub("HMP2_cd_RC2","HMP2 CD dietary pattern 2",testt$metadata)
testt$metadata = gsub("HMP2_uc_RC1","HMP2 UC dietary pattern 1",testt$metadata)
testt$metadata = gsub("HMP2_uc_RC2","HMP2 UC dietary pattern 2",testt$metadata)
testt$metadata = gsub("HMP2_nonIBD_RC1","HMP2 nonIBD dietary pattern 1",testt$metadata)
testt$metadata = gsub("HMP2_nonIBD_RC2","HMP2 nonIBD dietary pattern 2",testt$metadata)
testt$metadata = gsub("MLVS_RC1","MLVS dietary pattern 1",testt$metadata)
testt$metadata = gsub("MLVS_RC2","MLVS dietary pattern 2",testt$metadata)

ggplot(testt,aes(x=reorder(metadata,-abs(coef)), y=abs(coef),fill=metadata)) +
  geom_boxplot() +
  labs(y=expression(beta~'coefficients'), x="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) + 
  scale_fill_manual(values=c("indianred3","indianred3","azure3","azure3","lemonchiffon3","lemonchiffon3","azure4","azure4")) +
  theme(legend.position=c(.85,.6)) +
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.title.y = element_text(size = 10)) +
  theme(legend.title=element_blank()) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.key=element_blank()) + 
  theme(legend.position = 'none')
