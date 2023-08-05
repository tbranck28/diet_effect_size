library("RNOmni")
library("ggplot2")
library("grid")
library("vegan")
library(tidyverse)
library(lme4)
library(factoextra)
library(ade4)
library(FactoMineR)
library("dbplyr")
library("plyr")
library("phyloseq")
library(gtools)
library("viridis")
library("colorspace")
library("cowplot")
library("reshape2")
library("Maaslin2")
library(randomForest)
require(caTools)
library("ape")
#library(psych)
library(rstatix)
library(gridExtra)
library(ggthemes)
library(ggpubr)
library(pheatmap)

setwd("/Users/tobynbranck/Documents/Effect_size/analysis/")

#sample profiles processed through our pipeline
test = data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/infant_backhed/metaphlan_taxonomic_profiles.tsv", header=T, check.names = F)
names(test) <- gsub("Fixed_", "", names(test), fixed = TRUE)
names(test) <- gsub("_taxonomic_profile", "", names(test), fixed = TRUE)
names(test) <- gsub("# ", "", names(test), fixed = TRUE)

test = test[grep("s__", test$taxonomy),]
names(test) = gsub(pattern = "_.*", replacement = "", x = names(test))
test$taxonomy <- gsub('^.*?s__', '', test$taxonomy)
test$taxonomy = gsub("_", " ", test$taxonomy)

#divide by 100 to put on 0-1 scale
rownames(test) = test$taxonomy
test$taxonomy = NULL
bugnames = rownames(test)
test = test/100
rownames(test) = bugnames

#less strict filter for prevalence / abundance; 10^-5 in at least 25% of the samples
test$taxonomy = bugnames
filttest = test
filttest = data.frame(filttest)
rownames(filttest) = filttest$taxonomy
filttest$taxonomy = NULL

filttest = as.data.frame(t(filttest))

meta = data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/infant_backhed/mmc2.csv", header=T, check.names = F)
meta_id = data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/infant_backhed/studyID_person_association.csv", header=T, check.names = F)

meta_test = cbind(Individual = meta_id$Individual[ match(meta$`STUDY ID`, meta_id$`STUDY ID`) ],meta)
names(meta_test)[names(meta_test) == 'feeding practice first week'] = "_B"
names(meta_test)[names(meta_test) == 'feeding practice 4M'] = "_4M"
names(meta_test)[names(meta_test) == 'Any breastfeeding 12 M'] = "_12M"
meta_test$`Antibiotic treatment to infant 0-4M times` = ifelse(meta_test$`Antibiotic treatment to infant 0-4M times`==0, "no", "yes")
meta_test$`Antibiotic treatment to infant 4-12M times` = ifelse(meta_test$`Antibiotic treatment to infant 4-12M times`==0, "no", "yes")

#remove samples that have no response for feeding strategy
meta_test$feeding_history = paste(meta_test$'_B',meta_test$'_4M')
ggplot(meta_test, aes(x=as.factor(meta_test$feeding_history))) +
  geom_bar(fill='red') +  labs(x='feeding history') + scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + theme(text = element_text(size=5))
ggsave("feedinghistory_distribution.png")

age_table = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/filereport_read_run_PRJEB6456_tsv.txt", header=T, check.names = F))
age_table = age_table[,c("run_accession","sample_alias")]

age_table = age_table[!grepl("_M", age_table$sample_alias),]
accessions_keep = age_table$run_accession
filttest = filttest[rownames(filttest) %in% accessions_keep,]

filttest$Individual <- age_table$sample_alias[match(rownames(filttest), age_table$run_accession)]
names(filttest) <- gsub(x = names(filttest), pattern = "\\[", replacement = "_")  
names(filttest) <- gsub(x = names(filttest), pattern = "\\]", replacement = "_")
names(filttest) <- gsub(x = names(filttest), pattern = " ", replacement = "_")
#write.table(filttest, "/Users/tobynbranck/Documents/Effect_size/analysis/backhed_species.csv", sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

filttest_baseline = filttest[grep("_B", filttest$Individual), ]
filttest_four = filttest[grep("_4M", filttest$Individual), ]
filttest_twelve = filttest[grep("_12M", filttest$Individual), ]

#filtering by age group
filttest_baseline_IDs = filttest_baseline$Individual
filttest_four_IDs = filttest_four$Individual
filttest_twelve_IDs = filttest_twelve$Individual

rownames(filttest_baseline) = filttest_baseline_IDs
rownames(filttest_four) = filttest_four_IDs
rownames(filttest_twelve) = filttest_twelve_IDs

filttest_baseline = as.data.frame(t(filttest_baseline))
filttest_baseline <- mutate_all(filttest_baseline, function(x) as.numeric(as.character(x)))
filttest_four = as.data.frame(t(filttest_four))
filttest_four <- mutate_all(filttest_four, function(x) as.numeric(as.character(x)))
filttest_twelve = as.data.frame(t(filttest_twelve))
filttest_twelve <- mutate_all(filttest_twelve, function(x) as.numeric(as.character(x)))

filttest_baseline = filttest_baseline[rowSums(filttest_baseline > 0.00001) >= dim(filttest_baseline)[2]*.25, ]
filttest_four = filttest_four[rowSums(filttest_four > 0.00001) >= dim(filttest_four)[2]*.25, ]
filttest_twelve = filttest_twelve[rowSums(filttest_twelve > 0.00001) >= dim(filttest_twelve)[2]*.25, ]

filttest_baseline = as.data.frame(t(filttest_baseline[1:(dim(filttest_baseline)[1]-1),]))
filttest_four = as.data.frame(t(filttest_four[1:(dim(filttest_four)[1]-1),]))
filttest_twelve = as.data.frame(t(filttest_twelve[1:(dim(filttest_twelve)[1]-1),]))

#creating dataframes of metadata, one for each timepoint
baseline_meta = meta_test[,c("Individual","_B","GENDER","Delivery mode", "Age at sample Newborn (days)")]
baseline_meta$Individual <- paste0(baseline_meta$Individual, "_B")
baseline_meta = as.data.frame(baseline_meta)

four_meta = meta_test[,c("Individual","_4M","GENDER","Delivery mode","Age at Sampling \n4 M","Antibiotic treatment to infant 0-4M times")]
four_meta$Individual <- paste0(four_meta$Individual, "_4M")
four_meta = as.data.frame(four_meta)

twelve_meta = meta_test[,c("Individual","_12M","GENDER","Delivery mode","Age at Sampling 12M","feeding_history","Antibiotic treatment to infant 4-12M times")]
twelve_meta$Individual <- paste0(twelve_meta$Individual, "_12M")
twelve_meta = as.data.frame(twelve_meta)

filttest_baseline = filttest_baseline[order(rownames(filttest_baseline)), ]
filttest_four = filttest_four[order(rownames(filttest_four)), ]
filttest_twelve = filttest_twelve[order(rownames(filttest_twelve)), ]

baseline_meta = baseline_meta %>% drop_na()
four_meta = four_meta %>% drop_na()
twelve_meta = twelve_meta %>% drop_na()

baseline_meta = baseline_meta[order(baseline_meta$Individual), ]
four_meta = four_meta[order(four_meta$Individual), ]
twelve_meta = twelve_meta[order(twelve_meta$Individual), ]

baseline_meta = baseline_meta[!(is.na(baseline_meta$'_B') | baseline_meta$'_B'==""), ]
baseline_meta = baseline_meta[!baseline_meta$'_B' == "Formula feeding", ]
baseline_intersect <- intersect(rownames(filttest_baseline), baseline_meta$Individual)
filttest_baseline <- subset(filttest_baseline, rownames(filttest_baseline) %in% baseline_intersect)
baseline_meta <- subset(baseline_meta, Individual %in% baseline_intersect)

four_meta = four_meta[!(is.na(four_meta$`_4M`) | four_meta$`_4M`==""), ]
four_intersect <- intersect(rownames(filttest_four), four_meta$Individual)
filttest_four <- subset(filttest_four, rownames(filttest_four) %in% four_intersect)
four_meta <- subset(four_meta, Individual %in% four_intersect)

twelve_intersect <- intersect(rownames(filttest_twelve), twelve_meta$Individual)
filttest_twelve <- subset(filttest_twelve, rownames(filttest_twelve) %in% twelve_intersect)
twelve_meta <- subset(twelve_meta, Individual %in% twelve_intersect)

###### heatmap for time ######
#get union of top abundant species for each time point
baseline_vector = c(colnames(filttest_baseline))
four_vector = c(colnames(filttest_four))
twelve_vector = c(colnames(filttest_twelve))
bug_union = union(baseline_vector,union(four_vector,twelve_vector))

#subset original abundance table based on the union of bugs
filttest_for_heat = filttest
rownames(filttest_for_heat) = filttest_for_heat$Individual
names.use <- names(filttest_for_heat)[(names(filttest_for_heat) %in% bug_union)]
filttest_for_heat <- filttest_for_heat[, names.use]
filttest_for_heat <- mutate_all(filttest_for_heat, function(x) as.numeric(x))
filttest_for_heat$Group <- sub("^[^_]*_", "", rownames(filttest_for_heat))

rc <- ifelse(filttest_for_heat$Group == "B", "#440154FF",
             ifelse(filttest_for_heat$Group == "4M", "#31688EFF",
                    ifelse(filttest_for_heat$Group == "12M",
                           "#35B779FF", "#FDE725FF")))
filttest_for_heat$Group <- ordered(filttest_for_heat$Group, levels = c("B", "4M", "12M"))
heatmap_object = heatmap(as.matrix(filttest_for_heat[1:(length(filttest_for_heat)-1)]),
        RowSideColors = rc, Rowv = filttest_for_heat$Group, revC = TRUE, cexCol = 3/10,
        margins = c(5, 5))

#saving species df and metadata after 
filttest_for_heat2 = filttest_for_heat
base_comb = baseline_meta
four_comb = four_meta
twelve_comb = twelve_meta
names(base_comb)[2] = 'feeding'
names(four_comb)[2] = 'feeding'
names(twelve_comb)[2] = 'feeding'
names(base_comb)[5] = 'age'
names(four_comb)[5] = 'age'
names(twelve_comb)[5] = 'age'
base_comb$GENDER = NULL
four_comb$GENDER = NULL
twelve_comb$GENDER = NULL
base_comb$`Delivery mode` = NULL
four_comb$`Delivery mode` = NULL
twelve_comb$`Delivery mode` = NULL
four_comb$`Antibiotic treatment to infant 0-4M times`=NULL
twelve_comb$`Antibiotic treatment to infant 4-12M times` = NULL
twelve_comb$feeding_history = NULL
combined_meta = rbind(base_comb, four_comb, twelve_comb)
rownames(combined_meta) = combined_meta$Individual
#combined_meta$Individual = NULL
combined_meta = subset(combined_meta, rownames(combined_meta) %in% rownames(filttest_for_heat2))
filttest_for_heat2 = subset(filttest_for_heat2, rownames(filttest_for_heat2) %in% rownames(combined_meta))
write.csv(filttest_for_heat2, "/Users/tobynbranck/Documents/Effect_size/analysis/backhed_species_for_joint_ord.csv")
write.csv(combined_meta, "/Users/tobynbranck/Documents/Effect_size/analysis/backhed_meta_joint_ord.csv")

par(lend = 1)
legend("top", lwd = 8, box.lty=0,
       legend = c("Baseline", "4 Month", "12 Month"),
       col = c("#440154FF", "#31688EFF", "#35B779FF"),cex=.5)
annotation_row = subset(filttest_for_heat, select = c("Group"))
pdf("figure3/backhed_heatmap.pdf", width =12, height = 14)
pheatmap(t(filttest_for_heat[1:(length(filttest_for_heat)-1)]), cluster_cols = T, cluster_rows = T, annotation_col = annotation_row, fontsize_col = 6,show_rownames=T,cex=1,angle_col = "45")
dev.off()

######### beta diversity ordination plot###########
# FOUR MONTHS
#filttest_four[is.na(filttest_four)] <- 0
filttest_four_dist_matrix = vegdist(filttest_four, method="bray",diag=TRUE)
four_pcoa = pcoa(filttest_four_dist_matrix)
four_pcoa_df = as.data.frame(four_pcoa$vectors[,1:2])

NMS <-
  metaMDS(filttest_four_dist_matrix,
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)
data.scores = as.data.frame(scores(NMS))
data.scores$Feeding = four_meta$'_4M'
FOUR_nmds = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(alpha=8/10,size = 4, aes(colour = Feeding))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  labs(x = "NMDS1", colour = "4 Month Feeding", y = "NMDS2")  + 
  scale_colour_manual(values = c('lightskyblue4', 'maroon2','blue')) 

FOUR_pcplot = ggplot(four_pcoa_df, aes(x = Axis.1, y = Axis.2)) + 
  geom_point(alpha=8/10,size = 4, aes(colour = four_meta$'_4M'))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x=paste0("PC1: ",round(four_pcoa$values$Relative_eig[1],digits=3)*100,"%"), colour = "4 Month Feeding", y=paste0("PC2: ",round(four_pcoa$values$Relative_eig[2],digits=3)*100,"%"))  + 
  scale_colour_manual(values = c('lightskyblue4', 'maroon2','blue')) 

# BASELINE
filttest_baseline_dist_matrix = vegdist(filttest_baseline, method="bray",diag=TRUE)
baseline_pcoa = pcoa(filttest_baseline_dist_matrix)
baseline_pcoa_df = as.data.frame(baseline_pcoa$vectors[,1:2])

NMS_baseline <-
  metaMDS(filttest_baseline_dist_matrix,
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)
data.scores_baseline = as.data.frame(scores(NMS_baseline))
data.scores_baseline$Feeding = baseline_meta$'_B'
baseline_nmds = ggplot(data.scores_baseline, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(alpha=8/10,size = 4, aes(colour = Feeding))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  labs(x = "NMDS1", colour = "Baseline Feeding", y = "NMDS2")  + 
  scale_colour_manual(values = c('lightskyblue4', 'maroon2','blue')) 

baseline_pcplot = ggplot(baseline_pcoa_df, aes(x = Axis.1, y = Axis.2)) + 
  geom_point(alpha=8/10,size = 4, aes(colour = baseline_meta$'_B'))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x=paste0("PC1: ",round(baseline_pcoa$values$Relative_eig[1],digits=3)*100,"%"), colour = "Baseline Feeding", y=paste0("PC2: ",round(baseline_pcoa$values$Relative_eig[2],digits=3)*100,"%"))  + 
  scale_colour_manual(values = c('lightskyblue4', 'maroon2','blue')) 

# TWELVE
filttest_twelve_dist_matrix = vegdist(filttest_twelve, method="bray",diag=TRUE)
twelve_pcoa = pcoa(filttest_twelve_dist_matrix)
twelve_pcoa_df = as.data.frame(twelve_pcoa$vectors[,1:2])

NMS_twelve <-
  metaMDS(filttest_twelve_dist_matrix,
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)
data.scores_twelve = as.data.frame(scores(NMS_twelve))
data.scores_twelve$Feeding = twelve_meta$'_12M'
twelve_nmds = ggplot(data.scores_twelve, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(alpha=8/10,size = 4, aes(colour = Feeding))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  labs(x = "NMDS1", colour = "12 Month Feeding", y = "NMDS2")  + 
  scale_colour_manual(values = c('lightskyblue4', 'maroon2','blue')) 

twelve_pcplot = ggplot(twelve_pcoa_df, aes(x = Axis.1, y = Axis.2)) + 
  geom_point(alpha=8/10,size = 4, aes(colour = twelve_meta$'_12M'))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x=paste0("PC1: ",round(twelve_pcoa$values$Relative_eig[1],digits=3)*100,"%"), colour = "12 Month Feeding", y=paste0("PC2: ",round(twelve_pcoa$values$Relative_eig[2],digits=3)*100,"%"))  + 
  scale_colour_manual(values = c('lightskyblue4', 'maroon2','blue')) 

######### adonis section ###########
rownames(baseline_meta) = baseline_meta$Individual
names(baseline_meta)[5] = 'age'
names(baseline_meta)[2] = 'feeding'
baseline_month_adonis = adonis(filttest_baseline ~ baseline_meta$GENDER + baseline_meta$`Delivery mode` + baseline_meta$age + baseline_meta$feeding, data = baseline_meta, method = "bray",na.rm=TRUE)

rownames(four_meta) = four_meta$Individual
names(four_meta)[5] = 'age'
names(four_meta)[2] = 'feeding'
four_month_adonis = adonis(filttest_four ~ four_meta$`Antibiotic treatment to infant 0-4M times` + four_meta$GENDER + four_meta$`Delivery mode` + four_meta$age + four_meta$feeding, data = four_meta, method = "bray",na.rm=TRUE)

rownames(twelve_meta) = twelve_meta$Individual
names(twelve_meta)[5] = 'age'
names(twelve_meta)[2] = 'feeding'
twelve_month_adonis = adonis(filttest_twelve ~ twelve_meta$`Antibiotic treatment to infant 4-12M times` + twelve_meta$GENDER + twelve_meta$`Delivery mode` + twelve_meta$age + twelve_meta$feeding, data = twelve_meta, method = "bray",na.rm=TRUE)
twelve_month_adonis_FeedHist = adonis(filttest_twelve ~ twelve_meta$`Antibiotic treatment to infant 4-12M times` + twelve_meta$GENDER + twelve_meta$`Delivery mode` + twelve_meta$age + twelve_meta$feeding_history, data = twelve_meta, method = "bray",na.rm=TRUE)

combined_adonis_df = as.data.frame(rbind(baseline_month_adonis$aov.tab,four_month_adonis$aov.tab,twelve_month_adonis$aov.tab))
######### maaslin2 section ##########
baseline_meta$Individual = NULL
four_meta$Individual= NULL
twelve_meta$Individual = NULL
names(baseline_meta)[3] = 'Delivery_mode'
names(four_meta) <- gsub(x = names(four_meta), pattern = " ", replacement = "_")
names(four_meta) <- gsub(x = names(four_meta), pattern = "-", replacement = "_")
names(twelve_meta) <- gsub(x = names(twelve_meta), pattern = " ", replacement = "_")
names(twelve_meta) <- gsub(x = names(twelve_meta), pattern = "-", replacement = "_")
twelve_meta$feeding_history = NULL

#Maaslin2(input_data = filttest_baseline, input_metadata = baseline_meta, output = "/Users/tobynbranck/Documents/Effect_size/analysis/backhed_maaslin2_baseline_8.31",reference=c("feeding,exclusively breastfeeding"),max_significance = 1.0)
#Maaslin2(input_data = filttest_four, input_metadata = four_meta, output = "/Users/tobynbranck/Documents/Effect_size/analysis/backhed_maaslin2_four_8.31",reference=c("feeding,exclusively breastfeeding"),max_significance = 1.0)
#Maaslin2(input_data = filttest_twelve, input_metadata = twelve_meta, output = "/Users/tobynbranck/Documents/Effect_size/analysis/backhed_maaslin2_twelve_8.31",reference=c("feeding,any breastfeeding"),max_significance = 1.0)

write.csv(filttest_baseline, "/Users/tobynbranck/Documents/Effect_size/analysis/BackhedBaseline_bugs.csv")
write.csv(baseline_meta, "/Users/tobynbranck/Documents/Effect_size/analysis/BackhedBaseline_meta.csv")
write.csv(filttest_four, "/Users/tobynbranck/Documents/Effect_size/analysis/BackhedFour_bugs.csv")
write.csv(four_meta, "/Users/tobynbranck/Documents/Effect_size/analysis/BackhedFour_meta.csv")
write.csv(filttest_twelve, "/Users/tobynbranck/Documents/Effect_size/analysis/BackhedTwelve_bugs.csv")
write.csv(twelve_meta, "/Users/tobynbranck/Documents/Effect_size/analysis/BackhedTwelve_meta.csv")

###section testing only non-zero values###
filttest_twelve_ = filttest_twelve
filttest_twelve_$Ruminococcus_gnavus = NULL

abund_table_list = list(filttest_baseline,filttest_four,filttest_twelve)
meta_list = list(baseline_meta,four_meta,twelve_meta)
df_list = list()
for(s in seq(abund_table_list)){
  print(abund_table_list[s])
  abund_table = abund_table_list[s]
  meta_table = meta_list[s]
  
  bug_list = list()
  chi_sq_list = list()
  chi_pval_list = list()
  w_stat_list = list()
  w_pval_list = list()
  fisher_odds_list = list()
  fisher_pval_list = list()
  abund_table = as.data.frame(abund_table)
  meta_table = as.data.frame(meta_table)
  for(i in 1:ncol(abund_table)) {
    bug = names(abund_table)[i]
    bug_list = append(bug_list,bug)
    abund_table_no0 = abund_table
    abund_table_no0[abund_table_no0 == 0] <- 0
    abund_table_no0[abund_table_no0 != 0] <- 1
    abund_table_no0 = cbind(abund_table_no0,meta_table$feeding)
    names(abund_table_no0)[ncol(abund_table_no0)] = 'feeding'
    
    fisher_input = table(abund_table_no0$feeding,abund_table_no0[,bug])
    print(bug)
    print(fisher_input)
    fishertest = fisher.test(fisher_input)
    print(fishertest)
    fisher_odds = fishertest$estimate
    fisher_pval = fishertest$p.value
    fisher_odds_list = append(fisher_odds_list,fisher_odds)
    fisher_pval_list = append(fisher_pval_list,fisher_pval)
    print(fisher_pval_list)
  }
  if (length(fisher_odds_list)==0) {
    fisher_odds_list[1:ncol(abund_table_no0)-1] = NA
  }
  df = do.call(rbind,Map(data.frame, Bug=bug_list,Fisher_Odds=fisher_odds_list,Fisher_Pval=fisher_pval_list))
  df = as.data.frame(df)
  x = s
  write.csv(df,sprintf("fisher_results_%f.csv",x))
  df_list[[s]] = df
}
#wilcoxon code
abund_table_no0 = abund_table
abund_table_no0[abund_table_no0 == 0] <- NA
abund_table_no0 = cbind(abund_table_no0,meta_table$feeding)
names(abund_table_no0)[27] = 'feeding'
bf_values = abund_table_no0[,bug][abund_table_no0$feeding=='exclusively breastfeeding']
bf_values = bf_values[!is.na(bf_values)]
mixed_values = abund_table_no0[,bug][abund_table_no0$feeding=='mixed feeding']
mixed_values = mixed_values[!is.na(mixed_values)]
w_test = wilcox.test(bf_values,mixed_values)
w_stat = w_test$statistic
w_pval = w_test$p.value
w_stat_list = append(w_stat_list,w_stat)
w_pval_list = append(w_pval_list,w_pval)


#remove mixed feeding to better compare to cohen's d in backhed at al.
four_meta_NOmix = four_meta[!(four_meta$feeding=="mixed feeding"),]
#Maaslin2(input_data = filttest_four, input_metadata = four_meta_NOmix, output = "/Users/tobynbranck/Documents/Effect_size/analysis/backhed_maaslin2_four_8.31_NOmix",reference=c("feeding,exclusively breastfeeding"),max_significance = 1.0)

##### combining adonis results #####
adonis_results_combined = as.data.frame(rbind(baseline_month_adonis$aov.tab[4,5:6],four_month_adonis$aov.tab[5,5:6],twelve_month_adonis$aov.tab[5,5:6]))
rownames(adonis_results_combined) = gsub("\\$.*","",rownames(adonis_results_combined))
write.csv(adonis_results_combined, "/Users/tobynbranck/Documents/Effect_size/analysis/backhed_adonis_results.csv")

### wilcoxon rank sum - extra text not used for generating figures/manuscript components###
# four_month_wilx = cbind(filttest_four,four_meta$feeding)
# four_month_dfs_list = split(four_month_wilx, f=four_month_wilx$`four_meta$feeding`)
# exclusively_bf = as.data.frame(four_month_dfs_list[1])
# exclusively_ff = as.data.frame(four_month_dfs_list[2])
# names(exclusively_bf) = sub("exclusively.breastfeeding.", "", names(exclusively_bf))
# names(exclusively_ff) = sub("exclusively.formula.feeding.", "", names(exclusively_ff))
# exclusively_bf$four_meta.feeding = NULL
# exclusively_ff$four_meta.feeding = NULL
# 
# p.vals = list()
# bugs = list()
# cohens = list()
# for (i in 1:ncol(exclusively_bf)) {
#   combined = append(exclusively_bf[,i],exclusively_ff[,i])
#   label_col_bf = c(rep("bf",61))
#   label_col_ff = c(rep("ff",11))
#   label_col = append(label_col_bf,label_col_ff)
#   df = as.data.frame(cbind(label_col,combined))
#   bugs = append(bugs, names(exclusively_bf[i]))
#   w = wilcox.test(exclusively_bf[,i],exclusively_ff[,i], alternative = "two.sided")
#   p.vals = append(p.vals, w$p.value)
#   df$combined = as.numeric(df$combined)
#   cohenDstat = cohens_d(df,formula = combined ~ label_col, ref.group = "bf")
#   cohens = append(cohens,cohenDstat$effsize)
#   #cohenDstat = cohen.d(df,"label_col")
#   #cohens = append(cohens,cohenDstat$cohen.d[2])
#   #cohenDstat = cohen.d(df$combined,df$label_col)
#   #cohens = append(cohens,cohenDstat$estimate)
# }
# wilcoxon_results = do.call(rbind,Map(data.frame, P_value=p.vals,CohenD = cohens))
# rownames(wilcoxon_results) = bugs
# wilcoxon_results = wilcoxon_results[order(-wilcoxon_results$CohenD),]
# 
# wilcoxon_results %>%
#   rownames_to_column(var = "bugs") %>%
#   ggplot(aes(x = P_value, y = CohenD)) +
#   geom_point(color = "blue") +
#   geom_text(aes(label = rownames(wilcoxon_results)),size=3, hjust=0, nudge_x = 0.01) +
#   xlab("P value (Wilcoxon rank-sum)") +
#   ylab("Effect size (Cohen's d)") +
#   geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
#   geom_hline(yintercept=0.00, linetype="dashed", color = "red") +
#   theme_classic()
# write.csv(wilcoxon_results, file = "cohensd_results.csv")
# 
# ####BASELINE AND TWELVE MONTH WILCOXON AND COHEN'S D TESTING#####
# ####BASELINE#######
# names(baseline_meta)[2] = 'feeding'
# baseline_wilx = cbind(filttest_baseline,baseline_meta$feeding)
# baseline_dfs_list = split(baseline_wilx, f=baseline_wilx$`baseline_meta$feeding`)
# baseline_bf = as.data.frame(baseline_dfs_list[1])
# baseline_mixed = as.data.frame(baseline_dfs_list[2])
# names(baseline_bf) = sub("exclusively.breastfeeding.", "", names(baseline_bf))
# names(baseline_mixed) = sub("mixed.feeding.", "", names(baseline_mixed))
# baseline_bf$baseline_meta.feeding = NULL
# baseline_mixed$baseline_meta.feeding = NULL
# 
# p.vals_baseline = list()
# bugs_baseline = list()
# cohens_baseline = list()
# for (i in 1:ncol(baseline_bf)) {
#   combined_baseline = append(baseline_bf[,i],baseline_mixed[,i])
#   label_col_bf_baseline = c(rep("bf",59))
#   label_col_mixed_baseline = c(rep("mixed",21))
#   label_col_baseline = append(label_col_bf_baseline,label_col_mixed_baseline)
#   df_baseline = as.data.frame(cbind(label_col_baseline,combined_baseline))
#   bugs_baseline = append(bugs_baseline, names(baseline_bf[i]))
#   w_baseline = wilcox.test(baseline_bf[,i],baseline_mixed[,i], alternative = "two.sided")
#   p.vals_baseline = append(p.vals_baseline, w_baseline$p.value)
#   df_baseline$combined_baseline = as.numeric(df_baseline$combined_baseline)
#   cohenDstat_baseline = cohens_d(df_baseline,formula = combined_baseline ~ label_col_baseline, ref.group = "bf")
#   cohens_baseline = append(cohens_baseline,cohenDstat_baseline$effsize)
#   #cohenDstat = cohen.d(df,"label_col")
#   #cohens = append(cohens,cohenDstat$cohen.d[2])
#   #cohenDstat_baseline = cohen.d(df_baseline$combined_baseline,df_baseline$label_col_baseline)
#   #cohens_baseline = append(cohens_baseline,cohenDstat_baseline$estimate)
# }
# wilcoxon_results_baseline = do.call(rbind,Map(data.frame, P_value=p.vals_baseline,CohenD = cohens_baseline))
# rownames(wilcoxon_results_baseline) = bugs_baseline
# wilcoxon_results_baseline = wilcoxon_results_baseline[order(-wilcoxon_results_baseline$CohenD),]
# 
# wilcoxon_results_baseline %>%
#   rownames_to_column(var = "bugs") %>%
#   ggplot(aes(x = P_value, y = CohenD)) +
#   geom_point(color = "blue") +
#   geom_text(aes(label = rownames(wilcoxon_results_baseline)),size=3, hjust=0, nudge_x = 0.01) +
#   xlab("P value (Wilcoxon rank-sum)") +
#   ylab("Effect size (Cohen's d)") +
#   geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
#   geom_hline(yintercept=0.00, linetype="dashed", color = "red") +
#   theme_classic()
# write.csv(wilcoxon_results_baseline, file = "cohensd_results_baseline.csv")
# 
# #####TWELVE MONTHS#####
# twelve_wilx = cbind(filttest_twelve,twelve_meta$feeding)
# twelve_dfs_list = split(twelve_wilx, f=twelve_wilx$`twelve_meta$feeding`)
# twelve_bf = as.data.frame(twelve_dfs_list[1])
# twelve_nobf = as.data.frame(twelve_dfs_list[2])
# names(twelve_bf) = sub("any.breastfeeding.", "", names(twelve_bf))
# names(twelve_nobf) = sub("no.breastfeeding.feeding.", "", names(twelve_nobf))
# twelve_bf$twelve_meta.feeding = NULL
# twelve_nobf$twelve_meta.feeding = NULL
# 
# p.vals_twelve = list()
# bugs_twelve = list()
# cohens_twelve = list()
# for (i in 1:ncol(twelve_bf)) {
#   combined_twelve = append(twelve_bf[,i],twelve_nobf[,i])
#   label_col_bf_twelve = c(rep("bf",14))
#   label_col_nobf_twelve = c(rep("nobf",78))
#   label_col_twelve = append(label_col_bf_twelve,label_col_nobf_twelve)
#   df_twelve = as.data.frame(cbind(label_col_twelve,combined_twelve))
#   bugs_twelve = append(bugs_twelve, names(twelve_bf[i]))
#   w_twelve = wilcox.test(twelve_bf[,i],twelve_nobf[,i], alternative = "two.sided")
#   p.vals_twelve = append(p.vals_twelve, w_twelve$p.value)
#   df_twelve$combined_twelve = as.numeric(df_twelve$combined_twelve)
#   cohenDstat_twelve = cohens_d(df_twelve,formula = combined_twelve ~ label_col_twelve, ref.group = "bf")
#   cohens_twelve = append(cohens_twelve,cohenDstat_twelve$effsize)
#   #cohenDstat_twelve = cohen.d(df_twelve$combined_twelve,df_twelve$label_col_twelve)
#   #cohens_twelve = append(cohens_twelve,cohenDstat_twelve$estimate)
# }
# wilcoxon_results_twelve = do.call(rbind,Map(data.frame, P_value=p.vals_twelve,CohenD = cohens_twelve))
# rownames(wilcoxon_results_twelve) = bugs_twelve
# wilcoxon_results_twelve = wilcoxon_results_twelve[order(-wilcoxon_results_twelve$CohenD),]
# 
# wilcoxon_results_twelve %>%
#   rownames_to_column(var = "bugs") %>%
#   ggplot(aes(x = P_value, y = CohenD)) +
#   geom_point(color = "blue") +
#   geom_text(aes(label = rownames(wilcoxon_results_twelve)),size=3, hjust=0, nudge_x = 0.01) +
#   xlab("P value (Wilcoxon rank-sum)") +
#   ylab("Effect size (Cohen's d)") +
#   geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
#   geom_hline(yintercept=0.00, linetype="dashed", color = "red") +
#   theme_classic()
# write.csv(wilcoxon_results_twelve, file = "cohensd_results_twelve.csv")

### bug abundance line plots ###
#combine dataframes from baseline, 4, and 12 timepoints (separated because filtered
# by time point)
common <- intersect(names(filttest_baseline), names(filttest_four))
common = intersect(common,names(filttest_twelve))
all_bugs = rbind(filttest_baseline[,common], filttest_four[,common],filttest_twelve[,common])

common_meta = intersect(names(baseline_meta),names(four_meta))
common_meta = intersect(common_meta,names(twelve_meta))
all_meta = rbind(baseline_meta[,common_meta],four_meta[,common_meta],twelve_meta[,common_meta])
all_meta$time = gsub(".*_", "\\1", rownames(all_meta))
all_bugs$time = gsub(".*_", "\\1", rownames(all_bugs))
all_bugs$person = gsub("_.*", "\\1", rownames(all_bugs))

all_bugs$time = as.character(all_bugs$time)
all_bugs$time = factor(all_bugs$time, levels=c("B", "4M", "12M"))
all_meta$feeding[all_meta$feeding == "exclusively breastfeeding"] = "any breastfeeding"
all_meta$feeding[all_meta$feeding == "mixed feeding"] = "any breastfeeding"
all_meta$feeding[all_meta$feeding == "exclusively formula feeding"] = "no breastfeeding"
all_bugs$feeding = all_meta$feeding

test = all_bugs %>% 
  group_by(person) %>% 
  mutate(flag2 = if_else(lag(time) == 'B' & time == '4M' &  lag(feeding) != feeding,1, 0, 0))

test = test %>% 
  group_by(person) %>% 
  mutate(flag3 = if_else(lag(time) == '4M' & time == '12M' &  lag(feeding) != feeding,1, 0, 0))

test = test %>%
  mutate(Feeding_Change = case_when(
    flag2 == 1 ~ "B to 4 change",
    flag3 == 1 ~ "4 to 12 change"
  ))
test$Feeding_Change = test$Feeding_Change %>% replace_na("no change")
test$Feeding_Change = factor(test$Feeding_Change)
test = test %>% ungroup()
test = as.data.frame(test)

test$time = as.character(test$time)
test$time[test$time == "B"] = 0
test$time[test$time == "4M"] = 4
test$time[test$time == "12M"] = 12
test$time = as.numeric(test$time)

split_test = split(test,test$time)
df_baseline = as.data.frame(split_test[1])
df_four = as.data.frame(split_test[2])
df_twelve = as.data.frame(split_test[3])
names(df_baseline)[19] = "person"
names(df_four)[19] = "person"
names(df_twelve)[19] = "person"

first_segment = merge(df_baseline, df_four, by="person", all = T)
second_segment = merge(df_four, df_twelve, by="person", all = T)

myplots <- vector('list', ncol(test[,1:17]))
for (i in 1:ncol(test[,1:17])) {
  myplots[[i]] <- local({
    i <- i
    bugname = names(test[i])
    firstcol = gsub(" ", "", paste("X0.",bugname))
    secondcol = gsub(" ", "", paste("X4.",bugname))
    thirdcol = gsub(" ", "", paste("X12.",bugname))
    plt = ggplot(test, aes(time,test[,i])) +
      geom_point(aes(color=feeding), size=5, alpha=1) +
      scale_x_continuous(breaks = c(0, 4, 12)) +    
      geom_segment(data = first_segment,
                   aes(x = first_segment$X0.time,
                       y = first_segment[,firstcol],
                       xend = first_segment$X4.time,
                       yend = first_segment[,secondcol],
                       col = first_segment$X4.Feeding_Change)) +
      geom_segment(data = second_segment,
                   aes(x = second_segment$X4.time,
                       y = second_segment[,secondcol],
                       xend = second_segment$X12.time,
                       yend = second_segment[,thirdcol],
                       col = second_segment$X12.Feeding_Change)) +
      scale_color_manual(values=c("black","pink","black","purple","#E1E1E1"))+
      ylab(gsub("_", " ", bugname)) +
      theme_hc() +
      theme(legend.position = "none") +
      theme(axis.title.x=element_blank())
    print(plt)
})
}

cairo_pdf("figure3/backhed_abund_over_time.pdf", width=12, height = 16)
ggarrange(plotlist = myplots, ncol=3,nrow=6)
dev.off()

cairo_pdf("figure3/bacteroides_dorei_AbundOverTime.pdf", width=10, height = 6)
myplots[1]
dev.off()

cairo_pdf("figure3/collins_aero_AbundOverTime.pdf", width=10, height = 6)
myplots[12]
dev.off()

cairo_pdf("figure3/bifido_longum_AbundOverTime.pdf", width=10, height = 6)
myplots[15]
dev.off()

