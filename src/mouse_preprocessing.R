#This script analyzes the mouse microbiome dataset to
#understand the effects of diet on microbiome in a non-human gut environment.
#resulting number of samples after pairing taxonomic profiles with dietary metadata: 227 (62 complex diet; 165 refined diet)
#resulting number of individuals: 29

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

#set directory paths
input_dir = file.path("/home","tobynb","diet_test","input","Faith_mice","Initial",fsep = "/")
output_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")

#load taxonomy profiles and format the taxonomy table
tax = data.table::fread(paste0(input_dir,"/metaphlan_taxonomic_profiles.tsv"), header=T, check.names = F)
names(tax) <- gsub("# ", "", names(tax), fixed = TRUE)
tax = tax[grep("s__", tax$taxonomy),]
names(tax) = gsub(pattern = "_.*", replacement = "", x = names(tax))
tax$taxonomy <- gsub('^.*?s__', '', tax$taxonomy)
tax$taxonomy = gsub("_", " ", tax$taxonomy)

#plot to observe distribution of sum of abundance by sample - for samples with < 100, how to handle? why are some so low (in 80-98 range)??
qplot(colSums(tax[,-1]), geom="histogram")
#sample SRR097225 has all zeros
tax$SRR097225 = NULL

#divide by 100 to put on 0-1 scale
bugnames = tax$taxonomy
tax = tax[,2:dim(tax)[2]]/100
rownames(tax) = bugnames

#filter microbes for prevalence / abundance; 10^-5 in at least 25% of the samples
tax$taxonomy = bugnames
filttest = tax[rowSums(tax > 0.00001) >= dim(tax)[2]*.25, ]
filttest = data.frame(filttest)
rownames(filttest) = filttest$taxonomy
filttest$taxonomy = NULL
names(filttest) <- sub("^X", "", names(filttest))

filttest = as.data.frame(t(filttest))
samplist = rownames(filttest)
#write.table(filttest, paste0(output_dir,"/mouse_species.csv"), sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

#### Metadata curation section ###
meta = data.table::fread(paste0(input_dir,"/mouse_metadata.csv"), header=T, check.names = F) #From SRA table

#removing irrelevant samples (must have diet info, must be a part of the study that were colonized by the same mouse community)
meta = meta[meta$diet != "", ]
meta = meta[meta$sample_source=='Fecal pellet from a mouse (c57Bl6) colonized with 10-gut bacteria',]

#Extracting sample metadata (individual ID, study day, and diet period) from the 'Library Name' column in SRA metadata and add these extracted categories to the metadata
lib_name_strings = as.data.frame(stringr::str_split_fixed(meta$`Library Name`, "_", 3))
meta$individual = stringr::str_split_fixed(lib_name_strings$V3, " ", 2)[,1]
meta$sampling_day = lib_name_strings$V2
meta$diet_period = stringr::str_split_fixed(lib_name_strings$V1, " ", 2)[,2]
meta$diet_type = as.data.frame(stringr::str_split_fixed(meta$diet, "Diet ", 2))[2]
write.csv(meta,paste0(output_dir,"/mouse_all_metadata.csv")) #saving to have formatted metadata that has both complex and refined diet samples

#splitting the metadata into two different dataframes, one for the complex diets, one for the simple diets
simple_diet = meta %>%
  filter(str_detect(diet_type, "TD."))
complex_diet = meta[meta$diet_type=="",]

#keep samples with at least 500,000 depth
kneadd = data.table::fread(paste0(input_dir,"/kneaddata_read_count_table.tsv"), header=T, check.names = F)
kneadd = kneadd[kneadd$`final single`>=500000,]
simple_diet = simple_diet[simple_diet$Run %in% kneadd$Sample,]
complex_diet = complex_diet[complex_diet$Run %in% kneadd$Sample,]

#there are two refined diet sub-experiments denoted by diet_period
simple_diet_exp1 = simple_diet[simple_diet$diet_period %in% c('D2','D3','D4'),]
simple_diet_exp2 = simple_diet[simple_diet$diet_period %in% c('E1','E2','E3'),]

#subset the dataframe to include only relevant columns
simple_diet_exp1 = simple_diet_exp1[,c('Run','individual','sampling_day','diet_period','diet_type')]
simple_diet_exp2 = simple_diet_exp2[,c('Run','individual','sampling_day','diet_period','diet_type')]
complex_diet = complex_diet[,c('Run','Library Name','individual','sampling_day','diet_period','diet')]

#bring in simple and complex diet formulas, these were acquired from the Faith et al. supplement
simple_diets_formula = as.data.frame(data.table::fread(paste0(input_dir,"/mouse_simple_diets.txt"), sep = '\t', header = TRUE, check.names = F))
complex_diets_formula = as.data.frame(data.table::fread(paste0(input_dir,"/mouse_complex_diet.txt"), sep = '\t', header = TRUE, check.names = F))

#reshaping the complex diet dataframe
melted = melt(complex_diets_formula, id=c("V1","g/kg"))
foods <- c('apple sauce', 'peaches', 'chicken', 'sweet potatoes','oatmeal', 'peas', 'rice', 'beef')
melted[foods] <- lapply(foods, function(f) as.integer(melted$value == f))
melted$value = NULL

#remove even number weeks (these are the "reset" weeks where mice were given the refined diet to stabilize the microbiome)
melted <- melted[!grepl("week2|week4|week6|week8|week10", melted$V1), ]

# replace 1 in the row by the g/kg value
melted = melted %>% 
  dplyr::mutate_if(. %in% melted[,4:11], ~ dplyr::recode(., "1" = melted$`g/kg`, "0" = '0'))
melted[,4:11] <- sapply(melted[,4:11],as.numeric)

#collapse so that there is 1 row per week, per mouse
complex_diets_formula = aggregate(cbind(melted$`apple sauce`,melted$peaches,melted$chicken,melted$`sweet potatoes`,melted$oatmeal,melted$peas,melted$rice,melted$beef) ~ V1 + variable, data = melted, FUN = sum, na.rm = TRUE)
names(complex_diets_formula) = c("Week","Mouse","apple sauce","peaches","chicken","sweet potatoes","oatmeal","peas","rice","beef")

#recode week variables and individual variables
complex_diet <- complex_diet %>%
  mutate(
    diet_period = recode(diet_period,
                         'E9' = 'week1',
                         'E11' = 'week3',
                         'E13' = 'week5',
                         'E15' = 'week7',
                         'E17' = 'week9',
                         'E19' = 'week11'),
    individual = recode(individual,
                        'm1' = 'mouse1',
                        'm2' = 'mouse2',
                        'm3' = 'mouse3',
                        'm4' = 'mouse4',
                        'm5' = 'mouse5',
                        'm6' = 'mouse6',
                        'm7' = 'mouse7',
                        'm8' = 'mouse8')
  )


#merge individual and diet week columns to be able to get sample name 
complex_diet$identifier <- paste(complex_diet$individual,complex_diet$diet_period)
complex_diets_formula$identifier = paste(complex_diets_formula$Mouse, complex_diets_formula$Week)

#merge metadata diet intake with the rest of the metadata info
complex_diet <- complex_diet %>%
  left_join(
    complex_diets_formula %>%
      select(identifier, `apple sauce`, peaches, chicken, `sweet potatoes`,
             oatmeal, peas, rice, beef),
    by = "identifier"
  )

complex_diet_metadata = complex_diet[,c("Run","individual","sampling_day","diet_period","apple sauce","peaches","chicken","sweet potatoes","oatmeal","peas","rice","beef")]
rownames(complex_diet_metadata) = complex_diet_metadata$Run
complex_diet_metadata$Run = NULL

#save complex diet metadata
write.csv(complex_diet_metadata,paste0(output_dir,"/mouse_complex_diet_forSimul.csv"))

#checking my output
# complex_diet_metadata$Run = rownames(complex_diet_metadata)
# complex_diet_metadata = complex_diet_metadata[mixedorder(complex_diet_metadata$Run), ]
# mouse_complex_metadata = mouse_complex_metadata[mixedorder(mouse_complex_metadata$...1),]
# identical(mouse_complex_metadata$...1,complex_diet_metadata$Run)
# identical(mouse_complex_metadata$`apple sauce`,complex_diet_metadata$`apple sauce`)

#get common samples in filtered taxonomy table and complex_diet_metadata
filttest_complex = filttest[match(row.names(complex_diet_metadata), row.names(filttest)), ]
setdiff(rownames(filttest_complex), rownames(complex_diet_metadata))

#### complex diet exploratory plots
#plot complex diet compositions by mouse and week (bar plots) ####
for_bar = melt(complex_diet_metadata, id=c("individual","sampling_day","diet_period"))
for_bar$concat = paste(for_bar$individual,for_bar$diet_period,for_bar$variable)
for_bar = for_bar[order(for_bar[,'individual']),]
for_bar = for_bar[!duplicated(for_bar$concat),]
for_bar$sampling_day = NULL
for_bar$concat = paste(for_bar$individual,for_bar$diet_period)
colors = c( "#A54657",  "#582630", "#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A","#86BBD8")
pdf(file = paste0(output_dir,"/mouse_complex_diet_barplots.pdf"),width=12,height=9)
ggplot(for_bar, aes(x = factor(diet_period, levels=c("week1","week3","week5","week7","week9","week11")), fill = variable, y = value)) + 
  geom_bar(stat = "identity", colour = "black") + 
  theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1, face= "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold", colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "g/kg", fill = "diet variable") + 
  facet_wrap(~individual) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

#### plot heatmap of complex diet variables ####
mat_col = complex_diet_metadata[,"individual"]
rownames(mat_col) = rownames(complex_diet_metadata)
mat_colors = list(individual = RColorBrewer::brewer.pal(8,"Set3"))
names(mat_colors$individual) = unique(complex_diet_metadata$individual)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
#scaling variables prior to heatmap generation
z_diets_complex = as.data.frame(scale(complex_diet_metadata[,4:11]))
z_diets_complex = do.call(data.frame,lapply(z_diets_complex, function(x) replace(x, is.infinite(x),0)))
rownames(z_diets_complex) = rownames(complex_diet_metadata)
## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

pheatmap(z_diets_complex, annotation_row = mat_col, annotation_colors = mat_colors, border_color = NA,show_rownames = FALSE,annotation_names_row=FALSE, color = RColorBrewer::brewer.pal(6,"Purples"))


#PCoA plotting, ordinations overlayed by covariates
mouse_dist_matrix = vegdist(filttest_complex, method="bray",diag=TRUE)
mouse_pcoa = pcoa(mouse_dist_matrix)
mouse_pcoa_df = as.data.frame(mouse_pcoa$vectors[,1:2])
for (col in seq_along(complex_diet_metadata)){
  print(col)
  var = names(complex_diet_metadata)[col]
  print(var)
  plot = ggplot(mouse_pcoa_df, aes(x = Axis.1, y = Axis.2, color=factor(complex_diet_metadata[[col]]))) + 
    geom_point(alpha = 5/10, size=3) +
    labs(colour = "g/kg") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(x=paste0("PC1: ",round(mouse_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
         y=paste0("PC2: ",round(mouse_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
    ggtitle(paste("complex diet; colored by",var)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = c(.9, .25)) +
    theme(legend.key.size = unit(0.4, "cm")) +
    theme(legend.title = element_text(size = 8), 
          legend.text = element_text(size = 8))
  #ggsave(plot, file=paste0("/home/tobynb/diet_test/output_test/mouse_complex_ord_", var,".png"), width = 14, height = 10, units = "cm")
}

#ordination of samples labeled by individual mouse and annotated by sampling day
pdf(file = paste0(output_dir,"/mouse_sampsummary_ord.pdf"), width=10, height = 12)
ggplot(mouse_pcoa_df, aes(x = Axis.1, y = Axis.2, color=as.factor(complex_diet_metadata$individual))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Mouse cohort sample summary") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(.9, .25)) + 
  geom_line(aes(group = interaction(complex_diet_metadata$individual, complex_diet_metadata$diet_period)),color="black",size=.3) +
  geom_point(aes(shape = factor(complex_diet_metadata$sampling_day)), size = 3) +
  guides(col=guide_legend("individual"),
         shape=guide_legend("Sampling Day")) +
  labs(x=paste0("PC1: ",round(mouse_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(mouse_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
  geom_polygon(aes(group = interaction(complex_diet_metadata$individual, complex_diet_metadata$diet_period)),alpha=0.2) + theme(legend.title = element_text(size = 7),                                                                                                                                 legend.text = element_text(size = 7))
dev.off()  

#to show distribution of individuals' samples across diets and days after diet switch
ggplot(complex_diet_metadata, aes(y=interaction(individual,diet_period), x=as.factor(sampling_day)), color=as.factor(complex_diet_metadata$individual)) + 
  geom_point(stat="identity", aes(colour = factor(individual))) +
  scale_fill_viridis(discrete = T, option = "E") +
  ggtitle("") +
  theme(legend.position="none") +
  xlab("")

#number of samples per individual
numSamps_perInd = complex_diet_metadata %>% count("individual")
mean(numSamps_perInd$freq)
sum(numSamps_perInd$freq)

######complex diet analysis #######################

#adonis with raw, quantitative intake amounts
#initialize vectors which will store the adonis result
adonis_output = vector()
adonis_res_rsq = vector()
adonis_res_pval = vector()
#isolate diet variable metadata
complex_diet_metadata_adon = complex_diet_metadata[,4:11]
rownames(complex_diet_metadata_adon) = rownames(complex_diet_metadata)
#ensures the rownames for both matrices are identical
identical(rownames(filttest_complex),rownames(complex_diet_metadata_adon))
#loops through each complex diet variable and conducts an adonis test
for (col in seq_along(complex_diet_metadata_adon)){
  var_name = names(complex_diet_metadata_adon)[col]
  adonis.univ = adonis2(filttest_complex ~ complex_diet_metadata_adon[,var_name], data = complex_diet_metadata_adon, method = "bray", strata = complex_diet_metadata$individual,by="margin")
  print(MicEco::adonis_OmegaSq(adonis.univ, partial = TRUE))
  print(adonis.univ)
  adonis_output[col] = adonis.univ
  adonis_res_rsq[col] = adonis.univ[1,]$R2
  adonis_res_pval[col] = adonis.univ[1,]$`Pr(>F)`
  }
univar_res_rgnav = rbind(adonis_res_pval, adonis_res_rsq)
univar_res_rgnav = as.data.frame(t(univar_res_rgnav))
names(univar_res_rgnav) = c("P-Value", "R2")
univar_res_rgnav$`P-Value` = as.numeric(univar_res_rgnav$`P-Value`)
univar_res_rgnav$R2 = as.numeric(univar_res_rgnav$R2)
univar_res_rgnav$p_adj = p.adjust(univar_res_rgnav$`P-Value`, "fdr")
univar_res_rgnav$stars = cut(univar_res_rgnav$p_adj, c(0, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "-", ""))
univar_res_rgnav$diet_variable = colnames(complex_diet_metadata[,4:11])
univar_res_rgnav = univar_res_rgnav[order(univar_res_rgnav$R2),]

#saving complex_mouse metadata
write.table(complex_diet_metadata,paste0(output_dir,"/mouse_complex_metadata.csv"), sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(filttest_complex,paste0(output_dir,"/mouse_complex_species.csv"), sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

#maaslin2 analysis for complex diets
metadata_input_complex = complex_diet_metadata
metadata_input_complex$sampling_day = NULL
metadata_input_complex$diet_period = NULL

#removing spaces in column names
colnames(metadata_input_complex)[2] <- "apple_sauce"
colnames(metadata_input_complex)[5] <- "sweet_potatoes"

metadata_input_complex = as.data.frame(metadata_input_complex)
rownames(metadata_input_complex) = rownames(complex_diet_metadata)

#note: no significant results when run with quantitative diet vars - now categorizing into none, mid, or high categories
#Maaslin2(input_data = filttest_complex, input_metadata = metadata_input_complex[,2:9], output = paste0(output_dir,"/mouse_complexdiet_maaslin2"),random_effects = metadata_input_complex$individual)

metadata_input_complex[] <- lapply(metadata_input_complex, as.character)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "666.7", replacement = "High", fixed = TRUE)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "421.1", replacement = "High", fixed = TRUE)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "250", replacement = "Mid", fixed = TRUE)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "222.2", replacement = "Mid", fixed = TRUE)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "105.3", replacement = "Mid", fixed = TRUE)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "55.6", replacement = "Mid", fixed = TRUE)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "52.6", replacement = "Mid", fixed = TRUE)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "0", replacement = "None", fixed = TRUE)
write.table(metadata_input_complex,paste0(output_dir,"/mouse_complex_metadata_binned.csv"), sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

#Maaslin2(input_data = filttest_complex, input_metadata = metadata_input_complex, output = paste0(output_dir,"/mouse_complexdiet_maaslin2_binned"), fixed_effects = c('apple_sauce','peaches','chicken','sweet_potatoes','oatmeal','peas','rice','beef'),random_effects = c('individual'), reference=c("apple_sauce,None","peaches,None","chicken,None","sweet_potatoes,None","oatmeal,None","peas,None","rice,None","beef,None")) 

#adonis section using categorical variables
#initialize vectors which will store the adonis result
adonis_output = vector()
adonis_res_rsq = vector()
adonis_res_pval = vector()
#isolate diet variable metadata
complex_diet_metadata_adon = complex_diet_metadata[,4:11]
rownames(complex_diet_metadata_adon) = rownames(complex_diet_metadata)
#ensures the rownames for both matrices are identical
identical(rownames(filttest_complex),rownames(complex_diet_metadata_adon))
#loops through each complex diet variable and conducts an adonis test
for (col in seq_along(complex_diet_metadata_adon)){
  var_name = names(complex_diet_metadata_adon)[col]
  adonis.univ = adonis2(filttest_complex ~ complex_diet_metadata_adon[,var_name], data = complex_diet_metadata_adon, method = "bray", strata = complex_diet_metadata$individual,by="margin")
  print(MicEco::adonis_OmegaSq(adonis.univ, partial = TRUE))
  print(adonis.univ)
  adonis_output[col] = adonis.univ
  adonis_res_rsq[col] = adonis.univ[1,]$R2
  adonis_res_pval[col] = adonis.univ[1,]$`Pr(>F)`
}
univar_res_rgnav = rbind(adonis_res_pval, adonis_res_rsq)
univar_res_rgnav = as.data.frame(t(univar_res_rgnav))
names(univar_res_rgnav) = c("P-Value", "R2")
univar_res_rgnav$`P-Value` = as.numeric(univar_res_rgnav$`P-Value`)
univar_res_rgnav$R2 = as.numeric(univar_res_rgnav$R2)
univar_res_rgnav$p_adj = p.adjust(univar_res_rgnav$`P-Value`, "fdr")
univar_res_rgnav$stars = cut(univar_res_rgnav$p_adj, c(0, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "-", ""))
univar_res_rgnav$diet_variable = colnames(complex_diet_metadata[,4:11])
univar_res_rgnav = univar_res_rgnav[order(univar_res_rgnav$R2),]
write.csv(univar_res_rgnav,paste0(output_dir,"/mouse_complex_adonis.csv"))

############################################################################
##### Below is the code used for the refined diets shown in supplement #####
############################################################################

#format simple diets formula table
simple_diets_formula = as.data.frame(t(simple_diets_formula))
colnames(simple_diets_formula) = simple_diets_formula[1,]
simple_diets_formula = simple_diets_formula[-1,]
simple_diets_formula = simple_diets_formula[,-1]
simple_diets_formula$'diet_type' = rownames(simple_diets_formula)

#isolating the diets from the first simple diet sub-experiment and splitting into macro and micro nutrients
simple_diet1 = merge(simple_diet_exp1, simple_diets_formula, by = 'diet_type')
simple_diet1_micro = simple_diet1[,1:18]
simple_diet1_macro = cbind(simple_diet1[,1:5],simple_diet1[,19:25])
simple_diet1_micro = as.data.frame(simple_diet1_micro)
#remove "SRR097225"#
simple_diet1_micro = simple_diet1_micro[simple_diet1_micro$Run != "SRR097225", , drop = FALSE]
rownames(simple_diet1_micro) = simple_diet1_micro$Run
simple_diet1_micro$Run = NULL
simple_diet1_macro = as.data.frame(simple_diet1_macro)
rownames(simple_diet1_macro) = simple_diet1_macro$Run
simple_diet1_macro$Run = NULL

#isolating the diets from the second simple diet sub-experiment and splitting into macro and micro nutrients
simple_diet2 = merge(simple_diet_exp2, simple_diets_formula, by = 'diet_type')
simple_diet2_micro = simple_diet2[,1:18]
simple_diet2_macro = cbind(simple_diet2[,1:5],simple_diet2[,19:25])
simple_diet2_micro = as.data.frame(simple_diet2_micro)
rownames(simple_diet2_micro) = simple_diet2_micro$Run
simple_diet2_micro$Run = NULL
simple_diet2_macro = as.data.frame(simple_diet2_macro)
rownames(simple_diet2_macro) = simple_diet2_macro$Run
simple_diet2_macro$Run = NULL

#PCoA testing
#combine the two formatted sub-experiment dietary metadata tables
simple_diets_combined = as.data.frame(c(rownames(simple_diet1_micro), rownames(simple_diet2_micro)))
rownames(simple_diets_combined) = simple_diets_combined$`c(rownames(simple_diet1_micro), rownames(simple_diet2_micro))`
simple_diets_combined$`c(rownames(simple_diet1_micro), rownames(simple_diet2_micro))` = NULL
simple_diets_combined$Exp = NA
simple_diets_combined$Exp[1:79] = "exper1"
simple_diets_combined$Exp[80:165] = "exper2"
individual = c(simple_diet1_micro$individual,simple_diet2_micro$individual)
simple_diets_combined$individual = individual

#ensure the taxonomy table includes the matched samples found in the set of two simple diet subexperiments
spec = filttest[match(row.names(simple_diets_combined), row.names(filttest)), ]

#check to make sure the samples in the dataframes are identical
identical(rownames(spec),rownames(simple_diets_combined))

dist_matrix = vegdist(spec, method="bray",diag=TRUE)
pcoa = pcoa(dist_matrix)
pcoa_df = as.data.frame(pcoa$vectors[,1:2])

pdf(file = paste0(output_dir,"/simpleExperCombined.pdf"), width=6, height = 6)
ggplot(as.data.frame(pcoa_df), aes(x = Axis.1, y = Axis.2, color=simple_diets_combined$Exp)) + 
  geom_point(alpha = 4/10, size=3) +
  scale_color_manual(values = c('#E69F00','#56B4E9')) +
  labs(colour = "Sub-experiment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("PCoA - simple experiments") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(.9, .2)) +
  labs(x=paste0("PC1: ",round(pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(pcoa$values$Relative_eig[2],digits=3)*100,"%"))
dev.off()

adonis_obj_subexps = adonis2(spec ~ simple_diets_combined$Exp, data = simple_diets_combined, method = "bray",na.rm=TRUE, strata = simple_diets_combined$individual,by="margin")

#since there is no significant difference in species composition across the two refined diet sub-experiments, pooling the samples from the two sub-experiments 
diets_combined = rbind(simple_diet1_micro,simple_diet2_micro)

#ensuring the datatype for the dietary measurements is numeric
diets_combined[,5:17] <- lapply(diets_combined[,5:17], as.numeric)

#not enough variation in the below dietary variables, removing
diets_combined$`Maltodextrin Lo-Dex 10` = NULL
diets_combined$`Cellulose (Fiber)` = NULL
diets_combined$`Ethoxyquin (Liquid)` = NULL

#isolate diet variables for adonis testing
#diets_combined_adon = diets_combined[,5:14]

#### plot heatmap of simple combined diet variables #### - maybe split the 
#### map into macro vs. micro so that we can better tell the concentration
#### difference between the two, since the macro is masking the micro [ ]'s
mat_col = data.frame(individual = diets_combined$individual)
rownames(mat_col) = rownames(diets_combined)

mat_colors = RColorBrewer::brewer.pal(9,"RdPu")
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
z_diets_combined = as.data.frame(scale(diets_combined[,5:14]))
z_diets_combined = do.call(data.frame,lapply(z_diets_combined, function(x) replace(x, is.infinite(x),0)))
rownames(z_diets_combined) = rownames(diets_combined)
## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

cairo_pdf(file = paste0(output_dir,"/simple_diet_heatmap.pdf"), width=8, height = 7)
pheatmap(z_diets_combined, annotation_row = mat_col, border_color = NA,show_rownames = FALSE,annotation_names_row=FALSE, color = RColorBrewer::brewer.pal(9,"RdPu"))
dev.off()

#### Binning diets as opposed to running with raw values ###########
binned_diets = diets_combined
binned_diets$Casein[binned_diets$Casein<139] <- "low"
binned_diets$Casein[binned_diets$Casein>138 & binned_diets$Casein<600] <- "medium"
binned_diets$Casein[binned_diets$Casein!="low" & binned_diets$Casein!="medium"] <- "high"

binned_diets$`L-Cystine`[binned_diets$`L-Cystine`<2] <- "low"
binned_diets$`L-Cystine`[binned_diets$`L-Cystine`>2 & binned_diets$`L-Cystine`<40] <- "medium"
binned_diets$`L-Cystine`[binned_diets$`L-Cystine`!="low" & binned_diets$`L-Cystine`!="medium"] <- "high"

binned_diets$Sucrose[binned_diets$Sucrose<200] <- "low"
binned_diets$Sucrose[binned_diets$Sucrose>200 & binned_diets$Sucrose<500] <- "medium"
binned_diets$Sucrose[binned_diets$Sucrose!="low" & binned_diets$Sucrose!="medium"] <- "high"

binned_diets$`Corn Starch`[binned_diets$`Corn Starch`==0] <- "low"
binned_diets$`Corn Starch`[binned_diets$`Corn Starch`==100] <- "medium"
binned_diets$`Corn Starch`[binned_diets$`Corn Starch`==400] <- "high"

binned_diets$`Corn Oil`[binned_diets$`Corn Oil`<100] <- "low"
binned_diets$`Corn Oil`[binned_diets$`Corn Oil`!="low"] <- "high"

binned_diets$`79055 Mineral mix Ca-P Deficient`[binned_diets$`79055 Mineral mix Ca-P Deficient`<14] <- "low"
binned_diets$`79055 Mineral mix Ca-P Deficient`[binned_diets$`79055 Mineral mix Ca-P Deficient`==14.74] <- "medium"
binned_diets$`79055 Mineral mix Ca-P Deficient`[binned_diets$`79055 Mineral mix Ca-P Deficient`!="low" & binned_diets$`79055 Mineral mix Ca-P Deficient`!="medium"] <- "high"

binned_diets$`Calcium Carbonate`[binned_diets$`Calcium Carbonate`<6] <- "low"
binned_diets$`Calcium Carbonate`[binned_diets$`Calcium Carbonate`=="11.3"] <- "high"
binned_diets$`Calcium Carbonate`[binned_diets$`Calcium Carbonate`=="11.6"] <- "high"
binned_diets$`Calcium Carbonate`[binned_diets$`Calcium Carbonate`=="12.25"] <- "high"
binned_diets$`Calcium Carbonate`[binned_diets$`Calcium Carbonate`!="low" & binned_diets$`Calcium Carbonate`!= "high"] <- "medium"

binned_diets$`Calcium Phosphate`[binned_diets$`Calcium Phosphate`<4] <- "low"
binned_diets$`Calcium Phosphate`[binned_diets$`Calcium Phosphate`>4 & binned_diets$`Calcium Phosphate`<9] = "medium"
binned_diets$`Calcium Phosphate`[binned_diets$`Calcium Phosphate`!="low" & binned_diets$`Calcium Phosphate`!="medium"] <- "high"

binned_diets$`40077 Vitamin mix`[binned_diets$`40077 Vitamin mix`<16] <- "low"
binned_diets$`40077 Vitamin mix`[binned_diets$`40077 Vitamin mix`!="low"] <- "high"

binned_diets$`Choline Bitartrate`[binned_diets$`Choline Bitartrate`<2.4] <- "low"
binned_diets$`Choline Bitartrate`[binned_diets$`Choline Bitartrate`==2.40] <- "medium"
binned_diets$`Choline Bitartrate`[binned_diets$`Choline Bitartrate`==2.60] <- "high"

#subsetting metadata to include only diet variables for adonis2
binned_diets_adon = binned_diets[,5:14]

adonis_res_rsq_combined = vector()
adonis_res_pval_combined = vector()

#loop through diets (columns) and perform adonis2 test
for (col in seq_along(binned_diets_adon)){
  var_name = names(binned_diets_adon)[[col]]
  adonis_obj_subscombined = adonis2(spec ~ binned_diets_adon[,var_name], data = binned_diets_adon, method = "bray", strata = binned_diets$individual,na.rm=TRUE,by="margin")
  adonis_res_rsq_combined[col] = adonis_obj_subscombined[1,]$R2
  adonis_res_pval_combined[col] = adonis_obj_subscombined[1,]$`Pr(>F)`
  }

#extracting the R^2 and p-vals saved in the adonis_res_rsq_combined and adonis_res_pval_combined vectors and putting into dataframe
univar_res_rgnav_combined = rbind(adonis_res_pval_combined, adonis_res_rsq_combined)
univar_res_rgnav_combined = as.data.frame(t(univar_res_rgnav_combined))
names(univar_res_rgnav_combined) = c("P-Value", "R2")
univar_res_rgnav_combined$`P-Value` = as.numeric(univar_res_rgnav_combined$`P-Value`)
univar_res_rgnav_combined$R2 = as.numeric(univar_res_rgnav_combined$R2)
univar_res_rgnav_combined$p_adj = p.adjust(univar_res_rgnav_combined$`P-Value`, "fdr")
univar_res_rgnav_combined$stars = cut(univar_res_rgnav_combined$p_adj, c(0, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "-", ""))
univar_res_rgnav_combined$diet_variable = colnames(binned_diets_adon)
univar_res_rgnav_combined = univar_res_rgnav_combined[order(univar_res_rgnav_combined$R2),]

png(paste0(output_dir,"/simple_diets_combined_Binned_Analysis.png"))
ggplot(univar_res_rgnav_combined, aes(x=reorder(diet_variable,R2), y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() +
  geom_point(data = univar_res_rgnav_combined[univar_res_rgnav_combined$p_adj < 0.05, ],aes(x=reorder(diet_variable,R2), R2 + 0.01), shape = "*", size=4.233, color="black") +
  geom_point(data = univar_res_rgnav_combined[univar_res_rgnav_combined$p_adj < 0.1, ],aes(x=reorder(diet_variable,R2), R2 + 0.01), shape = "-", size=4.233, color="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  #ggtitle("Mouse Refined diets combined (binned)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(expression("R^2")) +
  xlab("Diet variable") +
  theme(plot.title = element_text(size=11))
dev.off()

simp_comb_adonis_binned = ggplot(univar_res_rgnav_combined, aes(x=reorder(diet_variable,R2), y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() +
  geom_point(data = univar_res_rgnav_combined[univar_res_rgnav_combined$p_adj < 0.1, ],aes(x=reorder(diet_variable,R2), R2 + 0.01), shape = "*", size=4.233, color="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  #ggtitle("Mouse Refined diets combined (binned)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(expression("R^2")) +
  xlab("") +
  theme(plot.title = element_text(size=11))

#formatting the binned diet tables for maaslin2
binned_diets_maaslin = binned_diets
binned_diets_maaslin$diet_type = NULL
binned_diets_maaslin$sampling_day = NULL
binned_diets_maaslin$diet_period = NULL
names(binned_diets_maaslin)<-make.names(names(binned_diets_maaslin),unique = TRUE)

spec_strict = spec
spec_strict$`Blautia coccoides` = NULL
spec_strict$`Lactococcus lactis` = NULL
spec_strict$`Desulfovibrio piger` = NULL

binned_diets_maaslin[,2:11] <- lapply(binned_diets_maaslin[,2:11], as.character)
#Maaslin2(input_data = spec_strict, input_metadata = binned_diets_maaslin[,2:11], output = paste0(output_dir,"/mouse_simple_diets_combined_character_binned"), random_effects = c("individual"), reference=c("Casein,low","L.Cystine,low","Sucrose,low","Corn.Starch,low","Corn.Oil,low","X79055.Mineral.mix.Ca.P.Deficient,low","Calcium.Carbonate,low","Calcium.Phosphate,low","40077.Vitamin.mix,low","Choline.Bitartrate,low"))

#put the two combine simple diet figures together for the Supplement
simple_combine_binned_EffectS <- as.data.frame(data.table::fread(paste0(output_dir,"/mouse_simple_diets_combined_character_binned/all_results.tsv"), header=T, check.names=T))
complex_results = as.data.frame(data.table::fread(paste0(output_dir,"/mouse_complexdiet_maaslin2_binned/all_results.tsv"), header=T, check.names=T))
simple_vec = simple_combine_binned_EffectS$coef
comp_vec = complex_results$coef
simple_forplot = data.frame(group = "Refined diet", value = simple_vec)
complex_forplot = data.frame(group = "Complex diet", value = comp_vec)
plot.data = rbind(simple_forplot,complex_forplot)
box = ggplot(plot.data, aes(x=group, y=abs(value))) + geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y=expression(beta~'coefficients'), x="")

t.test(abs(simple_forplot$value), abs(complex_forplot$value))
  
cairo_pdf(paste0(output_dir,"/mouse_refined_combined_supplement.pdf"),width = 10, height = 5)
ggarrange(simp_comb_adonis_binned, box, ncol = 2, nrow = 1,labels=c("A","B"))
dev.off()

#saving species table that includes samples from the complex diet + the two sub-experiment diets
length(union(rownames(spec),rownames(filttest_complex)))
filttest = filttest[rownames(filttest) %in% union(rownames(spec),rownames(filttest_complex)),]
write.table(filttest, paste0(output_dir,"/mouse_species.csv"), sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)


