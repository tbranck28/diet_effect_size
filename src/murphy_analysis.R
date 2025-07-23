#This script analyzes one of infant microbiome datasets (Murphy et al.) to
#understand the effects of diet (breastfeeding, intro of solid food) on microbiome
#resulting number of samples after pairing taxonomic profiles with dietary metadata: 247 total (baseline: 83; three months: 73; twelve months: 91)
#resulting number of individuals: 99

library(RNOmni)
library(ggplot2)
library(grid)
library(vegan)
library(tidyverse)
library(lme4)
library(factoextra)
library(ade4)
library(FactoMineR)
library(dbplyr)
library(plyr)
library(phyloseq)
library(gtools)
library(viridis)
library(colorspace)
library(cowplot)
library(reshape2)
library(Maaslin2)
library(randomForest)
require(caTools)
library(ape)
library(rstatix)
library(gridExtra)
library(ggthemes)
library(ggpubr)
library(pheatmap)

#set directory paths and load taxonomy and dietary profiles
input_dir = file.path("/home","tobynb","diet_test","input","Murphy_infants","Initial",fsep = "/")
output_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")

#load metadata
ids = as.data.frame(data.table::fread(paste0(input_dir,"/ids.txt"), header=T, check.names = F)) #acquired from corresponding author
meta = as.data.frame(data.table::fread(paste0(input_dir,"/metadata.csv"), header=T, check.names = F)) #acquired from corresponding author
bugs = as.data.frame(data.table::fread(paste0(input_dir,"/metaphlan_taxonomic_profiles.tsv"), header=T, check.names = F))
meta_SRA = as.data.frame(data.table::fread(paste0(input_dir,"/SraRunTable.txt"), header=T, check.names = F)) #acquired from SRA

#formatting sample IDs
names(bugs) <- gsub("_taxonomic_profile", "", names(bugs), fixed = TRUE)
names(bugs) <- gsub("re", "", names(bugs), fixed = TRUE)
names(bugs) <- gsub("# ", "", names(bugs), fixed = TRUE)

#grep for species-level taxonomy and formats taxonomic names; puts relative abundance values on a 0-1 scale
bugs = bugs[grep("s__", bugs$taxonomy),]
bugs$taxonomy <- gsub('^.*?s__', '', bugs$taxonomy)
bugs$taxonomy = gsub("_", " ", bugs$taxonomy)

bugs[,2:dim(bugs)[2]] = bugs[,2:dim(bugs)[2]]/100
rownames(bugs) = bugs$taxonomy
bugs$taxonomy = NULL
bugnames = rownames(bugs)
bugs = as.data.frame(t(bugs))

#adding age @ sample collection and additional ID columns to metadata
meta_test = merge(meta, ids, by.x='Studyid', by.y='Patient.ID')

#subsetting SRA metadata to include only the relevant columns and merging with other metadata table (meta_test); keeping only the placebo group
meta_SRA_subset = meta_SRA[, c("Run", "Library Name","host_sex","eth","Antibiotics_before_3_months","Antibiotics_before_6_months")]
meta_test2 = merge(meta_test, meta_SRA_subset, by.x = 'Otago.ID', by.y = 'Library Name') 
meta_test2 = meta_test2[meta_test2$studygroup == "placeb", ]

#Coding diet status at baseline, 3 months, 12 months, and 24 months old. bfduration = breastfeeding duration. ageanyformula = age at start of formula. solid = solid food.
#recoding breastfeeding duration to be on week scale
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

##### analysis of baseline samples #####
#filter for baseline samples and order rows based on sample names
baseline_meta = meta_test2 %>% dplyr::filter(Time == 0)
baseline_meta = baseline_meta[order(baseline_meta$Run),]

#distribution of feeding types: exclusively breastfed = 146, mixed_feeding = 28, no breastfeeding = 2 - removing the 2 samples of "no breastfeeding"
baseline_meta <-baseline_meta[!(baseline_meta$baseline_feeding=="no_breastfeeding"),]

#ensure the removal one of the duplicate samples in a pair if occurs from the same individual / same time point.
#very few here, 3 duplicate pairs at baseline; 1 duplicate pair at 3 months; none at 12 months
#if the majority of individuals had duplicate samples within the same time point, I would keep the samples and run
#the analyses accordingly i.e., strata within individual, random effects etc. Choosing one sample from each duplicate pair to drop.
n_occur <- data.frame(table(baseline_meta$Studyid))
n_occur[n_occur$Freq > 1,]
droplist_baseline = c('SRR7352040','SRR7351797','SRR7351697')
baseline_meta = baseline_meta[!baseline_meta$Run %in% droplist_baseline, ]

#order rows of taxonomy table
filtbugs_baseline = bugs[order(rownames(bugs)),]

#keep only the samples that have baseline metadata
intersect_baseline = intersect(baseline_meta$Run, rownames(filtbugs_baseline))
baseline_meta = subset(baseline_meta, Run %in% intersect_baseline)
filtbugs_baseline = subset(filtbugs_baseline,rownames(filtbugs_baseline) %in% intersect_baseline)

#filter microbes for prevalence and abundance
filtbugs_baseline = filtbugs_baseline[,colSums(filtbugs_baseline > 0.00001) >= dim(filtbugs_baseline)[1]*.25]
rownames(baseline_meta) = baseline_meta$Run
baseline_meta$Run = NULL
identical(rownames(baseline_meta),rownames(filtbugs_baseline)) #double-check rows are ordered identically

#adonis2 function
baseline_adonis = adonis2(filtbugs_baseline ~ baseline_meta$host_sex + baseline_meta$caesar + baseline_meta$baseline_feeding, data = baseline_meta, method = "bray", na.rm=TRUE, by="margin")

#maaslin2 for baseline samples
baseline_month_meta_subset = baseline_meta[,c("host_sex","baseline_feeding","caesar")]
Maaslin2(input_data = filtbugs_baseline, input_metadata = baseline_month_meta_subset, output = paste0(output_dir,"/murphy_maaslin2_baseline_8.31"),reference = c("baseline_feeding,exclusively breastfed"),max_significance = 1)
write.csv(filtbugs_baseline, paste0(output_dir,"/MurphyBaseline_bugs.csv"))
write.csv(baseline_month_meta_subset, paste0(output_dir,"/MurphyBaseline_meta.csv"))

##### analysis of 3 month samples #####
#filter for three month samples and order rows based on sample names
three_month_meta <- meta_test2 %>% dplyr::filter(Time == 3)
three_month_meta = three_month_meta[order(three_month_meta$Run), ]

#distribution of feeding types: exclusively breastfed = 41, mixed_feeding = 24, no_breastfeeding = 9, 
#mixed_bf_formual_solid = 1; removing "mixed_bf_formula_solid" because only 1 sample
three_month_meta <-three_month_meta[!(three_month_meta$three_month_feeding=="mixed_bf_solid_formula"),]

#ensure the removal one of the duplicate samples in a pair if occurs from the same individual / same time point.
n_occur_3 <- data.frame(table(three_month_meta$Studyid))
n_occur_3[n_occur_3$Freq > 1,]
droplist_three_month = c('SRR7351953')
three_month_meta = three_month_meta[!three_month_meta$Run %in% droplist_three_month, ]

#order rows of taxonomy table
filtbugs = bugs[order(rownames(bugs)),]

#keep only the three month samples
intersect_three <- intersect(three_month_meta$Run, rownames(filtbugs))
three_month_meta <- subset(three_month_meta, Run %in% intersect_three)
filtbugs <- subset(filtbugs, rownames(filtbugs) %in% intersect_three)

#filter microbes for prevalence and abundance
filtbugs = filtbugs[,colSums(filtbugs > 0.00001) >= dim(filtbugs)[1]*.25]
rownames(three_month_meta) = three_month_meta$Run
three_month_meta$Run = NULL
identical(rownames(three_month_meta),rownames(filtbugs)) #double-check rows are ordered identically

three_month_adonis = adonis2(filtbugs ~ three_month_meta$Antibiotics_before_3_months + three_month_meta$host_sex + three_month_meta$caesar + three_month_meta$three_month_feeding, data = three_month_meta, method = "bray",na.rm=TRUE, by="margin")

#maaslin2 for three months
three_month_meta_subset = three_month_meta[,c("host_sex","Antibiotics_before_3_months","three_month_feeding","caesar")]
Maaslin2(input_data = filtbugs, input_metadata = three_month_meta_subset, output = paste0(output_dir,"/murphy_maaslin2_3_month_8.31"),reference=c("three_month_feeding,exclusively breastfed"),max_significance = 1)
write.csv(filtbugs, paste0(output_dir,"/MurphyThreeMonth_bugs.csv"))
write.csv(three_month_meta_subset, paste0(output_dir,"/MurphyThreeMonth_meta.csv"))

##### analysis of 12 month samples #####
#filter for baseline samples and order rows based on sample names
twelve_month_meta <- meta_test2 %>% dplyr::filter(Time == 12)
twelve_month_meta = twelve_month_meta[order(twelve_month_meta$Run), ]
#distribution of feeding types: mixed_bf_formula_solid = 47, mixed_bf_solid = 25, no_breastfeeding = 122, 

#ensure the removal one of the duplicate samples in a pair if occurs from the same individual / same time point
n_occur_12 <- data.frame(table(twelve_month_meta$Studyid))
n_occur_12[n_occur_12$Freq > 1,]

#order rows of taxonomy table
filtbugs_twelve = bugs[order(rownames(bugs)),]

#keep only the three month samples
intersect_twelve <- intersect(twelve_month_meta$Run, rownames(filtbugs_twelve))
twelve_month_meta <- subset(twelve_month_meta, Run %in% intersect_twelve)
filtbugs_twelve <- subset(filtbugs_twelve, rownames(filtbugs_twelve) %in% intersect_twelve)

#filter microbes for prevalence and abundance
filtbugs_twelve = filtbugs_twelve[,colSums(filtbugs_twelve > 0.00001) >= dim(filtbugs_twelve)[1]*.25]
rownames(twelve_month_meta) = twelve_month_meta$Run
twelve_month_meta$Run = NULL

#adonis2 function
twelve_month_adonis = adonis2(filtbugs_twelve ~ twelve_month_meta$Antibiotics_before_6_months + twelve_month_meta$host_sex + twelve_month_meta$caesar + twelve_month_meta$twelve_month_feeding, data = twelve_month_meta, method = "bray",na.rm=TRUE,by="margin")

#maaslin2 for twelve months
twelve_month_meta_subset = twelve_month_meta[,c("host_sex","twelve_month_feeding","Antibiotics_before_6_months","caesar")]
Maaslin2(input_data = filtbugs_twelve, input_metadata = twelve_month_meta_subset, output = paste0(output_dir,"/murphy_maaslin2_twelve_months_8.31"),reference=c("twelve_month_feeding,mixed_bf_solid_formula"),max_significance = 1)
write.csv(filtbugs_twelve, paste0(output_dir,"/MurphyTwelveMonth_bugs.csv"))
write.csv(twelve_month_meta_subset, paste0(output_dir,"/MurphyTwelveMonth_meta.csv"))

##### analysis of 24 month samples #####
#filter for baseline samples and order rows based on sample names - NOT enough variation in feeding type to analyze
twentyfour_month_meta <- meta_test2 %>% dplyr::filter(Time == 24)
twentyfour_month_meta = twentyfour_month_meta[order(twentyfour_month_meta$Run), ]
#distribution of feeding types: mixed_bf_formula_solid = 1, mixed_bf_solid = 3, no_breastfeeding = 116; NOT enough variation in feeding type to analyze

#combining adonis results to save as a file
adonis_results_combined = as.data.frame(rbind(baseline_adonis[3,3:5],three_month_adonis[4,3:5],twelve_month_adonis[4,3:5]))
rownames(adonis_results_combined) = gsub("\\$.*","",rownames(adonis_results_combined))
write.csv(adonis_results_combined, paste0(output_dir,"/infants_murphy_adonis_results.csv"))

###### heatmap for time ######
#get union of top abundant species for each time point
baseline_vector = c(colnames(filtbugs_baseline))
three_vector = c(colnames(filtbugs))
twelve_vector = c(colnames(filtbugs_twelve))
bug_union = union(baseline_vector, union(three_vector,twelve_vector))

#subset original abundance table based on the union of bugs
filttest_for_heat = bugs
names.use <- names(filttest_for_heat)[(names(filttest_for_heat) %in% bug_union)]
filttest_for_heat <- filttest_for_heat[, names.use]
filttest_for_heat <- mutate_all(filttest_for_heat, function(x) as.numeric(x))
age = meta_SRA[,c("Run","Host_Age")]
rownames(age) = age$Run
age$Run = NULL
filttest_for_heat = as.data.frame(merge(filttest_for_heat,age,by=0))
rownames(filttest_for_heat) = filttest_for_heat$Row.names
filttest_for_heat$Row.names = NULL

filttest_for_heat<-filttest_for_heat[!(filttest_for_heat$Host_Age=="24_months"),]
annotation_row = subset(filttest_for_heat, select = c("Host_Age"))
pdf(paste0(output_dir,"/murphy_heatmap.pdf"), width = 12, height = 14)
pheatmap(t(filttest_for_heat[1:(length(filttest_for_heat)-2)]), cluster_cols = T, cluster_rows = T, cex=1, annotation_col = annotation_row, fontsize_col = 6,show_rownames=T,angle_col = "45")
dev.off()

#formatting all-time point metadata table and taxonomy table to save for figure preparation
filttest_for_heat2 = filttest_for_heat
base_comb = as.data.frame(baseline_meta[c("Studyid","baseline_feeding","Time")])
three_comb = as.data.frame(three_month_meta[c("Studyid","three_month_feeding","Time")])
twelve_comb = as.data.frame(twelve_month_meta[c("Studyid","twelve_month_feeding","Time")])
names(base_comb)[2] = 'feeding'
names(three_comb)[2] = 'feeding'
names(twelve_comb)[2] = 'feeding'
combined_meta = rbind(base_comb, three_comb, twelve_comb)
names(combined_meta)[1] = 'Individual'
names(combined_meta)[3] = 'age'
combined_meta = combined_meta %>%
  mutate(age = case_when(
    age == 0 ~ "baseline",
    age == 3 ~ "three months",
    age == 12 ~ "twelve months"))
combined_meta = subset(combined_meta, rownames(combined_meta) %in% rownames(filttest_for_heat2))
filttest_for_heat2 = subset(filttest_for_heat2, rownames(filttest_for_heat2) %in% rownames(combined_meta))
write.csv(filttest_for_heat2, paste0(output_dir,"/murphy_species_for_joint_ord.csv"))
write.csv(combined_meta, paste0(output_dir,"/murphy_meta_joint_ord.csv"))