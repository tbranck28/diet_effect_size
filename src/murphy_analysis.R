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
library(rstatix)
library(pheatmap)

ids = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/Infants_murphy/ids.txt", header=T, check.names = F))
meta = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/Infants_murphy/metadata.csv", header=T, check.names = F))
bugs = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/Infants_murphy/metaphlan_taxonomic_profiles.tsv", header=T, check.names = F))
meta_SRA = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/Infants_murphy/SraRunTable.txt", header=T, check.names = F))

names(bugs) <- gsub("_taxonomic_profile", "", names(bugs), fixed = TRUE)
names(bugs) <- gsub("re", "", names(bugs), fixed = TRUE)
names(bugs) <- gsub("# ", "", names(bugs), fixed = TRUE)

bugs = bugs[grep("s__", bugs$taxonomy),]
#names(bugs) = gsub(pattern = "_.*", replacement = "", x = names(bugs))
bugs$taxonomy <- gsub('^.*?s__', '', bugs$taxonomy)
bugs$taxonomy = gsub("_", " ", bugs$taxonomy)
rownames(bugs) = bugs$taxonomy
bugs$taxonomy = NULL
bugnames = rownames(bugs)
bugs = bugs/100
rownames(bugs) = bugnames

#less strict filter for prevalence / abundance; 10^-5 in at least 25% of the samples
bugs$taxonomy = bugnames
#filteredbugs = bugs[rowSums(bugs > 0.00001) >= dim(bugs)[2]*.25, ]
filteredbugs = bugs
filteredbugs = data.frame(filteredbugs)
rownames(filteredbugs) = filteredbugs$taxonomy
filteredbugs$taxonomy = NULL
filteredbugs = as.data.frame(t(filteredbugs))
write.table(filteredbugs, "/Users/tobynbranck/Documents/Effect_size/analysis/murphy_species.csv", sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

meta_test = merge(meta, ids, by.x='Studyid', by.y='Patient.ID')
meta_SRA_subset = meta_SRA[, c("Run", "Library Name","host_sex","eth","Antibiotics_before_3_months","Antibiotics_before_6_months")]
meta_test2 = merge(meta_test, meta_SRA_subset, by.x = 'Otago.ID', by.y = 'Library Name') 
meta_test2 = meta_test2[meta_test2$studygroup == "placeb", ]
# split metadata dataframe by exclusively breastfeeding; subset meta_test2 by samples who's bf_duration >= 3 months and
# ageanyformula and solid food = 0
meta_test2$bfduration <- as.numeric(meta_test2$bfduration) * 4

#coding baseline feeding. not including "solid" food because earliest solid food was 10 weeks
meta_test2 = meta_test2 %>%
  mutate(baseline_feeding = case_when(bfduration>=0 & ageanyformula>0 ~ "exclusively breastfed", 
                                      bfduration>=0 & is.na(ageanyformula) ~ "exclusively breastfed", 
                                      bfduration>=0 & ageanyformula==0 ~ "mixed_feeding",
                                      is.na(bfduration) ~ "no_breastfeeding"))
# three month feeding (12 weeks)
meta_test2 = meta_test2 %>%
  mutate(three_month_feeding = case_when(bfduration>=12 & ageanyformula>12 & solid>12 ~ "exclusively breastfed",
                                         bfduration>=12 & is.na(ageanyformula) & solid>12 ~ "exclusively breastfed", 
                                         bfduration>=12 & ageanyformula<=12 & solid>12 ~ "mixed_feeding",
                                         bfduration>=12 & ageanyformula<=12 & solid<=12 ~ "mixed_bf_solid_formula",
                                         bfduration<12 ~ "no_breastfeeding",
                                         is.na(bfduration) ~ "no_breastfeeding"))

# twelve month feeding (48 weeks) ####combining mixed_bf_formula_solid and mixed_bf_solid into mixed_bf_solid_formula
meta_test2 = meta_test2 %>%
  mutate(twelve_month_feeding = case_when(bfduration>=48 & ageanyformula>48 & solid>48 ~ "exclusively breastfed",
                                          bfduration>=48 & is.na(ageanyformula) & solid>48 ~ "exclusively breastfed", 
                                          bfduration>=48 & ageanyformula<=48 & solid>48 ~ "mixed_feeding",
                                          bfduration>=48 & ageanyformula<=48 & solid<=48 ~ "mixed_bf_solid_formula",
                                          bfduration>=48 & ageanyformula>48 & solid<=48 ~ "mixed_bf_solid_formula",
                                          bfduration>=48 & is.na(ageanyformula) & solid<=48 ~ "mixed_bf_solid_formula",
                                          bfduration<48 ~ "no_breastfeeding",
                                          is.na(bfduration) ~ "no_breastfeeding"))
# twenty-four month feeding (96 weeks)
meta_test2 = meta_test2 %>%
  mutate(twentyfour_month_feeding = case_when(bfduration>=96 & ageanyformula>96 & solid>96 ~ "exclusively breastfed",
                                              bfduration>=96 & is.na(ageanyformula) & solid>96 ~ "exclusively breastfed", 
                                              bfduration>=96 & ageanyformula<=96 & solid>96 ~ "mixed_feeding",
                                              bfduration>=96 & ageanyformula<=96 & solid<=96 ~ "mixed_bf_solid_formula",
                                              bfduration>=96 & ageanyformula>96 & solid<=96 ~ "mixed_bf_solid_formula",
                                              bfduration>=96 & is.na(ageanyformula) & solid<=96 ~ "mixed_bf_solid_formula",
                                              bfduration<96 ~ "no_breastfeeding",
                                              is.na(bfduration) ~ "no_breastfeeding"))

##### baseline samples #####
baseline_meta = meta_test2 %>% dplyr::filter(Time == 0)
baseline_meta = baseline_meta[order(baseline_meta$Run),]
#distribution of feeding types: exclusively breastfed = 146, mixed_feeding = 28, no breastfeeding = 2, removing the 2 samples of "no breastfeeding"
baseline_meta <-baseline_meta[!(baseline_meta$baseline_feeding=="no_breastfeeding"),]
filtbugs_baseline = filteredbugs[order(rownames(filteredbugs)),]
intersect_baseline = intersect(baseline_meta$Run, rownames(filtbugs_baseline))
baseline_meta = subset(baseline_meta, Run %in% intersect_baseline)
filtbugs_baseline = subset(filtbugs_baseline,rownames(filtbugs_baseline) %in% intersect_baseline)
filtbugs_baseline = filtbugs_baseline[,colSums(filtbugs_baseline > 0.00001) >= dim(filtbugs_baseline)[1]*.25]
rownames(baseline_meta) = baseline_meta$Run
baseline_meta$Run = NULL

baseline_adonis = adonis(filtbugs_baseline ~ baseline_meta$host_sex + baseline_meta$caesar + baseline_meta$baseline_feeding, data = baseline_meta, method = "bray",na.rm=TRUE)

#maaslin2 for baseline
baseline_month_meta_subset = baseline_meta[,c("host_sex","baseline_feeding","caesar")]
Maaslin2(input_data = filtbugs_baseline, input_metadata = baseline_month_meta_subset, output = "/Users/tobynbranck/Documents/Effect_size/analysis/murphy_maaslin2_baseline_8.31",reference = c("baseline_feeding,exclusively breastfed"),max_significance = 1)
write.csv(filtbugs_baseline, "/Users/tobynbranck/Documents/Effect_size/analysis/MurphyBaseline_bugs.csv")
write.csv(baseline_month_meta_subset, "/Users/tobynbranck/Documents/Effect_size/analysis/MurphyBaseline_meta.csv")

##### 3 month samples #####
three_month_meta <- meta_test2 %>% dplyr::filter(Time == 3)
three_month_meta = three_month_meta[order(three_month_meta$Run), ]
#distribution of feeding types: exclusively breastfed = 41, mixed_feeding = 24, no_breastfeeding = 9, 
#mixed_bf_formual_solid = 1; removing "mixed_bf_formula_solid" because only 1 samples
three_month_meta <-three_month_meta[!(three_month_meta$three_month_feeding=="mixed_bf_solid_formula"),]
filtbugs = filteredbugs[order(rownames(filteredbugs)),]
intersect_three <- intersect(three_month_meta$Run, rownames(filtbugs))
three_month_meta <- subset(three_month_meta, Run %in% intersect_three)
filtbugs <- subset(filtbugs, rownames(filtbugs) %in% intersect_three)
filtbugs = filtbugs[,colSums(filtbugs > 0.00001) >= dim(filtbugs)[1]*.25]
rownames(three_month_meta) = three_month_meta$Run
three_month_meta$Run = NULL

three_month_adonis = adonis(filtbugs ~ three_month_meta$Antibiotics_before_3_months + three_month_meta$host_sex + three_month_meta$caesar + three_month_meta$three_month_feeding, data = three_month_meta, method = "bray",na.rm=TRUE)

###removing the below to show that weird ordination shape due to bugs with many zeros - extra text for exploration
# filtbugs$`Aeriscardovia aeriphila`=NULL
# filtbugs$`Collinsella aerofaciens` = NULL
# filtbugs$`Bifidobacterium bifidum`= NULL
# filtbugs$`Bifidobacterium breve`=NULL
# filtbugs$`Bifidobacterium longum`=NULL
# filttest_three_dist_matrix = vegdist(filtbugs, method="bray",diag=TRUE)
# three_pcoa = pcoa(filttest_three_dist_matrix)
# three_pcoa_df = as.data.frame(three_pcoa$vectors[,1:2])
# 
# NMS <-
#   metaMDS(filttest_three_dist_matrix ,
#           distance = "bray",
#           k = 3,
#           maxit = 999, 
#           trymax = 500,
#           wascores = TRUE)
# data.scores = as.data.frame(scores(NMS))
# data.scores$Feeding = three_month_meta$three_month_feeding
# three_nmds = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
#   geom_point(alpha=8/10,size = 4, aes(colour = Feeding))+ 
#   theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
#         axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
#         legend.text = element_text(size = 12, face ="bold", colour ="black"), 
#         legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
#         axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
#         legend.title = element_text(size = 14, colour = "black", face = "bold"), 
#         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#         legend.key=element_blank()) +
#   labs(x = "NMDS1", colour = "3 Month Feeding", y = "NMDS2")  + 
#   scale_colour_manual(values = c('lightskyblue4', 'maroon2','blue',"red")) 
# 
# three_pcplot = ggplot(three_pcoa_df, aes(x = Axis.1, y = Axis.2)) + 
#   geom_point(alpha=8/10,size = 4, aes(colour = three_month_meta$three_month_feeding))+ 
#   theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
#         axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
#         legend.text = element_text(size = 12, face ="bold", colour ="black"), 
#         legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
#         axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
#         legend.title = element_text(size = 14, colour = "black", face = "bold"), 
#         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#         legend.key=element_blank()) + 
#   geom_text(aes(label = rownames(three_pcoa_df)),size=3, hjust=0, nudge_x = 0.01) +
#   labs(x=paste0("PC1: ",round(three_pcoa$values$Relative_eig[1],digits=3)*100,"%"), colour = "3 Month Feeding", y=paste0("PC2: ",round(three_pcoa$values$Relative_eig[2],digits=3)*100,"%"))  + 
#   scale_colour_manual(values = c('lightskyblue4', 'maroon2','blue',"red")) 
# 
# three_pcoa_df_sub = three_pcoa_df[three_pcoa_df$Axis.2 < -.2, ]
# row.names.remove <- c("SRR7351955", "SRR7411292")
# three_pcoa_df_sub = three_pcoa_df_sub[!(row.names(three_pcoa_df_sub) %in% row.names.remove), ]
# filtbugs_sub = filtbugs[(row.names(filtbugs) %in% row.names(three_pcoa_df_sub)),]

#maaslin2 for three months
three_month_meta_subset = three_month_meta[,c("host_sex","Antibiotics_before_3_months","three_month_feeding","caesar")]
Maaslin2(input_data = filtbugs, input_metadata = three_month_meta_subset, output = "/Users/tobynbranck/Documents/Effect_size/analysis/murphy_maaslin2_3_month_8.31",reference=c("three_month_feeding,exclusively breastfed"),max_significance = 1)
write.csv(filtbugs, "/Users/tobynbranck/Documents/Effect_size/analysis/MurphyThreeMonth_bugs.csv")
write.csv(three_month_meta_subset, "/Users/tobynbranck/Documents/Effect_size/analysis/MurphyThreeMonth_meta.csv")

# 12 month samples 
twelve_month_meta <- meta_test2 %>% dplyr::filter(Time == 12)
twelve_month_meta = twelve_month_meta[order(twelve_month_meta$Run), ]
#distribution of feeding types: mixed_bf_formula_solid = 47, mixed_bf_solid = 25, no_breastfeeding = 122, 
filtbugs_twelve = filteredbugs[order(rownames(filteredbugs)),]
intersect_twelve <- intersect(twelve_month_meta$Run, rownames(filtbugs_twelve))
twelve_month_meta <- subset(twelve_month_meta, Run %in% intersect_twelve)
filtbugs_twelve <- subset(filtbugs_twelve, rownames(filtbugs_twelve) %in% intersect_twelve)
filtbugs_twelve = filtbugs_twelve[,colSums(filtbugs_twelve > 0.00001) >= dim(filtbugs_twelve)[1]*.25]
rownames(twelve_month_meta) = twelve_month_meta$Run
twelve_month_meta$Run = NULL

twelve_month_adonis = adonis(filtbugs_twelve ~ twelve_month_meta$Antibiotics_before_6_months + twelve_month_meta$host_sex + twelve_month_meta$caesar + twelve_month_meta$twelve_month_feeding, data = twelve_month_meta, method = "bray",na.rm=TRUE)

#maaslin2 for twelve months
twelve_month_meta_subset = twelve_month_meta[,c("host_sex","twelve_month_feeding","Antibiotics_before_6_months","caesar")]
Maaslin2(input_data = filtbugs_twelve, input_metadata = twelve_month_meta_subset, output = "/Users/tobynbranck/Documents/Effect_size/analysis/murphy_maaslin2_twelve_months_8.31",reference=c("twelve_month_feeding,mixed_bf_solid_formula"),max_significance = 1)
write.csv(filtbugs_twelve, "/Users/tobynbranck/Documents/Effect_size/analysis/MurphyTwelveMonth_bugs.csv")
write.csv(twelve_month_meta_subset, "/Users/tobynbranck/Documents/Effect_size/analysis/MurphyTwelveMonth_meta.csv")

# 24 month samples 
twentyfour_month_meta <- meta_test2 %>% dplyr::filter(Time == 24)
twentyfour_month_meta = twentyfour_month_meta[order(twentyfour_month_meta$Run), ]
#distribution of feeding types: mixed_bf_formula_solid = 1, mixed_bf_solid = 3, no_breastfeeding = 116; NOT enough variation in feeding type to analyze

##putting all together
adonis_results_combined = as.data.frame(rbind(baseline_adonis$aov.tab[2,5:6],three_month_adonis$aov.tab[3,5:6],twelve_month_adonis$aov.tab[2,5:6]))
rownames(adonis_results_combined) = gsub("\\$.*","",rownames(adonis_results_combined))
write.csv(adonis_results_combined, "/Users/tobynbranck/Documents/Effect_size/analysis/infants_murphy_adonis_results.csv")

###### heatmap for time ######
#get union of top abundant species for each time point
baseline_vector = c(colnames(filtbugs_baseline))
three_vector = c(colnames(filtbugs))
twelve_vector = c(colnames(filtbugs_twelve))
bug_union = union(baseline_vector, union(three_vector,twelve_vector))

#subset original abundance table based on the union of bugs
filttest_for_heat = filteredbugs
names.use <- names(filttest_for_heat)[(names(filttest_for_heat) %in% bug_union)]
filttest_for_heat <- filttest_for_heat[, names.use]
filttest_for_heat <- mutate_all(filttest_for_heat, function(x) as.numeric(x))
age = meta_SRA[,c("Run","Host_Age")]
rownames(age) = age$Run
age$Run = NULL
filttest_for_heat = as.data.frame(merge(filttest_for_heat,age,by=0))
rownames(filttest_for_heat) = filttest_for_heat$Row.names
filttest_for_heat$Row.names = NULL

#filttest_for_heat$Group <- sub("^[^_]*_", "", rownames(filttest_for_heat))
filttest_for_heat<-filttest_for_heat[!(filttest_for_heat$Host_Age=="24_months"),]
annotation_row = subset(filttest_for_heat, select = c("Host_Age"))
pdf("figure3/murphy_heatmap.pdf", width = 12, height = 14)
pheatmap(t(filttest_for_heat[1:(length(filttest_for_heat)-2)]), cluster_cols = T, cluster_rows = T, annotation_col = annotation_row, fontsize_col = 6,show_rownames=T,cex=1,angle_col = "45")
dev.off()

filttest_for_heat2 = filttest_for_heat

base_comb = as.data.frame(baseline_meta[c("Studyid","baseline_feeding")])
three_comb = as.data.frame(three_month_meta[c("Studyid","three_month_feeding")])
twelve_comb = as.data.frame(twelve_month_meta[c("Studyid","twelve_month_feeding")])
names(base_comb)[2] = 'feeding'
names(three_comb)[2] = 'feeding'
names(twelve_comb)[2] = 'feeding'
combined_meta = rbind(base_comb, three_comb, twelve_comb)
names(combined_meta)[1] = 'Individual'
combined_meta$age = NA
combined_meta$age[1:86] = 'baseline'
combined_meta$age[87:160] = 'three months'
combined_meta$age[161:251] = 'twelve months'
combined_meta = subset(combined_meta, rownames(combined_meta) %in% rownames(filttest_for_heat2))
filttest_for_heat2 = subset(filttest_for_heat2, rownames(filttest_for_heat2) %in% rownames(combined_meta))
write.csv(filttest_for_heat2, "/Users/tobynbranck/Documents/Effect_size/analysis/murphy_species_for_joint_ord.csv")
write.csv(combined_meta, "/Users/tobynbranck/Documents/Effect_size/analysis/murphy_meta_joint_ord.csv")


######### beta diversity ordination plot###########
## THREE MONTHS
#filttest_four[is.na(filttest_four)] <- 0
filttest_three_dist_matrix = vegdist(filtbugs, method="bray",diag=TRUE)
three_pcoa = pcoa(filttest_three_dist_matrix)
three_pcoa_df = as.data.frame(three_pcoa$vectors[,1:2])

NMS <-
  metaMDS(filttest_three_dist_matrix,
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)
data.scores = as.data.frame(scores(NMS))
data.scores$Feeding = three_month_meta$three_month_feeding
three_nmds = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(alpha=8/10,size = 4, aes(colour = Feeding))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  labs(x = "NMDS1", colour = "3 Month Feeding", y = "NMDS2")  + 
  scale_colour_manual(values = c('lightskyblue4', 'maroon2','blue',"red")) 

three_pcplot = ggplot(three_pcoa_df, aes(x = Axis.1, y = Axis.2)) + 
  geom_point(alpha=8/10,size = 4, aes(colour = three_month_meta$three_month_feeding))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x=paste0("PC1: ",round(three_pcoa$values$Relative_eig[1],digits=3)*100,"%"), colour = "3 Month Feeding", y=paste0("PC2: ",round(three_pcoa$values$Relative_eig[2],digits=3)*100,"%"))  + 
  scale_colour_manual(values = c('lightskyblue4', 'maroon2','blue',"red")) 

# BASELINE
filttest_baseline_dist_matrix = vegdist(filtbugs_baseline, method="bray",diag=TRUE)
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
data.scores_baseline$Feeding = baseline_meta$baseline_feeding
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
  geom_point(alpha=8/10,size = 4, aes(colour = baseline_meta$feeding))+ 
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
filttest_twelve_dist_matrix = vegdist(filtbugs_twelve, method="bray",diag=TRUE)
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
data.scores_twelve$Feeding = twelve_month_meta$twelve_month_feeding
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
  geom_point(alpha=8/10,size = 4, aes(colour = twelve_month_meta$twelve_month_feeding))+ 
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

### wilcoxon rank sum ###
three_month_wilx = cbind(filtbugs,three_month_meta$three_month_feeding)
three_month_dfs_list = split(three_month_wilx, f=three_month_wilx$'three_month_meta$three_month_feeding')
exclusively_bf = as.data.frame(three_month_dfs_list[1])
exclusively_ff = as.data.frame(three_month_dfs_list[3])
names(exclusively_bf) = sub("exclusively.breastfed.", "", names(exclusively_bf))
names(exclusively_ff) = sub("no_breastfeeding.", "", names(exclusively_ff))
exclusively_bf$three_month_meta.three_month_feeding = NULL
exclusively_ff$three_month_meta.three_month_feeding = NULL

p.vals = list()
bugs = list()
cohens = list()
for (i in 1:ncol(exclusively_bf)) {
  combined = append(exclusively_bf[,i],exclusively_ff[,i])
  label_col_bf = c(rep("bf",61))
  label_col_ff = c(rep("ff",11))
  label_col = append(label_col_bf,label_col_ff)
  df = as.data.frame(cbind(label_col,combined))
  bugs = append(bugs, names(exclusively_bf[i]))
  w = wilcox.test(exclusively_bf[,i],exclusively_ff[,i], alternative = "two.sided")
  p.vals = append(p.vals, w$p.value)
  df$combined = as.numeric(df$combined)
  cohenDstat = cohens_d(df,formula = combined ~ label_col, ref.group = "bf")
  cohens = append(cohens,cohenDstat$effsize)
  #cohenDstat = cohen.d(df,"label_col")
  #cohens = append(cohens,cohenDstat$cohen.d[2])
  #cohenDstat = cohen.d(df$combined,df$label_col)
  #cohens = append(cohens,cohenDstat$estimate)
}
wilcoxon_results = do.call(rbind,Map(data.frame, P_value=p.vals,CohenD = cohens))
rownames(wilcoxon_results) = bugs
wilcoxon_results = wilcoxon_results[order(-wilcoxon_results$CohenD),]

wilcoxon_results %>%
  rownames_to_column(var = "bugs") %>%
  ggplot(aes(x = P_value, y = CohenD)) +
  geom_point(color = "blue") +
  geom_text(aes(label = rownames(wilcoxon_results)),size=3, hjust=0, nudge_x = 0.01) +
  xlab("P value (Wilcoxon rank-sum)") +
  ylab("Effect size (Cohen's d)") +
  geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
  geom_hline(yintercept=0.00, linetype="dashed", color = "red") +
  theme_classic()
write.csv(wilcoxon_results, file = "cohensd_results_murhpy_3months.csv")

####BASELINE AND TWELVE MONTH WILCOXON AND COHEN'S D TESTING#####
####BASELINE#######
names(baseline_meta)[21] = 'feeding'
baseline_wilx = cbind(filtbugs_baseline,baseline_meta$feeding)
baseline_dfs_list = split(baseline_wilx, f=baseline_wilx$`baseline_meta$feeding`)
baseline_bf = as.data.frame(baseline_dfs_list[1])
baseline_mixed = as.data.frame(baseline_dfs_list[2])
names(baseline_bf) = sub("exclusively.breastfed.", "", names(baseline_bf))
names(baseline_mixed) = sub("mixed.feeding.", "", names(baseline_mixed))
baseline_bf$baseline_meta.feeding = NULL
baseline_mixed$baseline_meta.feeding = NULL

p.vals_baseline = list()
bugs_baseline = list()
cohens_baseline = list()
for (i in 1:ncol(baseline_bf)) {
  combined_baseline = append(baseline_bf[,i],baseline_mixed[,i])
  label_col_bf_baseline = c(rep("bf",59))
  label_col_mixed_baseline = c(rep("nobf",21))
  label_col_baseline = append(label_col_bf_baseline,label_col_mixed_baseline)
  df_baseline = as.data.frame(cbind(label_col_baseline,combined_baseline))
  bugs_baseline = append(bugs_baseline, names(baseline_bf[i]))
  w_baseline = wilcox.test(baseline_bf[,i],baseline_mixed[,i], alternative = "two.sided")
  p.vals_baseline = append(p.vals_baseline, w_baseline$p.value)
  df_baseline$combined_baseline = as.numeric(df_baseline$combined_baseline)
  cohenDstat_baseline = cohens_d(df_baseline,formula = combined_baseline ~ label_col_baseline, ref.group = "bf")
  cohens_baseline = append(cohens_baseline,cohenDstat_baseline$effsize)
  #cohenDstat = cohen.d(df,"label_col")
  #cohens = append(cohens,cohenDstat$cohen.d[2])
  #cohenDstat_baseline = cohen.d(df_baseline$combined_baseline,df_baseline$label_col_baseline)
  #cohens_baseline = append(cohens_baseline,cohenDstat_baseline$estimate)
}
wilcoxon_results_baseline = do.call(rbind,Map(data.frame, P_value=p.vals_baseline,CohenD = cohens_baseline))
rownames(wilcoxon_results_baseline) = bugs_baseline
wilcoxon_results_baseline = wilcoxon_results_baseline[order(-wilcoxon_results_baseline$CohenD),]

wilcoxon_results_baseline %>%
  rownames_to_column(var = "bugs") %>%
  ggplot(aes(x = P_value, y = CohenD)) +
  geom_point(color = "blue") +
  geom_text(aes(label = rownames(wilcoxon_results_baseline)),size=3, hjust=0, nudge_x = 0.01) +
  xlab("P value (Wilcoxon rank-sum)") +
  ylab("Effect size (Cohen's d)") +
  geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
  geom_hline(yintercept=0.00, linetype="dashed", color = "red") +
  theme_classic()
write.csv(wilcoxon_results_baseline, file = "cohensd_results_murphy_baseline.csv")

#####TWELVE MONTHS#####
twelve_wilx = cbind(filtbugs_twelve,twelve_month_meta$twelve_month_feeding)
twelve_dfs_list = split(twelve_wilx, f=twelve_wilx$`twelve_month_meta$twelve_month_feeding`)
twelve_mixed = as.data.frame(twelve_dfs_list[1])
twelve_nobf = as.data.frame(twelve_dfs_list[2])
names(twelve_mixed) = sub("mixed_bf_solid_formula.", "", names(twelve_mixed))
names(twelve_nobf) = sub("no_breastfeeding.", "", names(twelve_nobf))
twelve_mixed$twelve_month_meta.twelve_month_feeding = NULL
twelve_nobf$twelve_month_meta.twelve_month_feeding = NULL

p.vals_twelve = list()
bugs_twelve = list()
cohens_twelve = list()
for (i in 1:ncol(twelve_mixed)) {
  combined_twelve = append(twelve_mixed[,i],twelve_nobf[,i])
  label_col_bf_twelve = c(rep("bf",59))
  label_col_nobf_twelve = c(rep("nobf",21))
  label_col_twelve = append(label_col_bf_twelve,label_col_nobf_twelve)
  df_twelve = as.data.frame(cbind(label_col_twelve,combined_twelve))
  bugs_twelve = append(bugs_twelve, names(twelve_mixed[i]))
  w_twelve = wilcox.test(twelve_mixed[,i],twelve_nobf[,i], alternative = "two.sided")
  p.vals_twelve = append(p.vals_twelve, w_twelve$p.value)
  df_twelve$combined_twelve = as.numeric(df_twelve$combined_twelve)
  cohenDstat_twelve = cohens_d(df_twelve,formula = combined_twelve ~ label_col_twelve, ref.group = "bf")
  cohens_twelve = append(cohens_twelve,cohenDstat_twelve$effsize)
  #cohenDstat_twelve = cohen.d(df_twelve$combined_twelve,df_twelve$label_col_twelve)
  #cohens_twelve = append(cohens_twelve,cohenDstat_twelve$estimate)
}
wilcoxon_results_twelve = do.call(rbind,Map(data.frame, P_value=p.vals_twelve,CohenD = cohens_twelve))
rownames(wilcoxon_results_twelve) = bugs_twelve
wilcoxon_results_twelve = wilcoxon_results_twelve[order(-wilcoxon_results_twelve$CohenD),]

wilcoxon_results_twelve %>%
  rownames_to_column(var = "bugs") %>%
  ggplot(aes(x = P_value, y = CohenD)) +
  geom_point(color = "blue") +
  geom_text(aes(label = rownames(wilcoxon_results_twelve)),size=3, hjust=0, nudge_x = 0.01) +
  xlab("P value (Wilcoxon rank-sum)") +
  ylab("Effect size (Cohen's d)") +
  geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
  geom_hline(yintercept=0.00, linetype="dashed", color = "red") +
  theme_classic()
write.csv(wilcoxon_results_twelve, file = "cohensd_results_murphy_twelve.csv")

###getting number individuals in study (for samples that were used)