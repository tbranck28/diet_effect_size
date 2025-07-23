#This script formats the taxonomic profile and metadata tables from the HMP2 cohort
#resulting number of samples after pairing taxonomic profiles with dietary metadata: 770
#resulting number of individuals: 63

#load packages
library(RNOmni)
library(ggplot2)
library(grid)
library(vegan)
library(tidyverse)
library(lme4)
library(factoextra)
library(ade4)
library(FactoMineR)
library(dplyr)
library(stringr)

#setting paths for input and output directories
input_dir = file.path("/home","tobynb","diet_test","input",fsep = "/")
output_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")

#load metadata and taxonomy tables
hmp2_metadata = data.table::fread(paste0(input_dir,"/HMP2/Initial/hmp2_metadata.csv"), header=T, check.names = F)
hmp2_extra_meta = data.table::fread(paste0(input_dir,"/HMP2/Initial/HMP2_MGX_metadata_09282018.csv"), header=T, check.names = F)
hmp2_species = data.table::fread(paste0(input_dir,"/HMP2/Initial/hmp2_species.txt"), header=T, check.names = F)
 
##subset metadata files
#grep for information that pertains to the metagenomic samples (there are several other sample types in this metadata)
hmp2_metadata = hmp2_metadata[hmp2_metadata$data_type == 'metagenomics',]

#subset to include consent age, abx, and state of dysbiosis
hmp2_extra_meta = hmp2_extra_meta[,c('External_ID','consent_age','antibiotic','diagnosis2')]

#subset to include diet info only (there is a non-diet column embedded here, so dropping that)
hmp2_meta = hmp2_metadata[,72:92] %>% select(-"1) In the past 2 weeks, have you received any of the following medications")
hmp2_meta$`Tea and coffee` = hmp2_metadata[,'Tea or coffee no sugar and no sugar replacement']
diet_vars = colnames(hmp2_meta)
write.csv(diet_vars, paste0(output_dir,"/hmp2_diet_vars_list.csv"))

disease = c(hmp2_metadata$diagnosis)

#coding diet normalized by serving's per day
old_values = c("No, I did not consume these products in the last 7 days",
               "Within the past 4 to 7 days",
               "Within the past 2 to 3 days",
               "Yesterday, 1 to 2 times",
               "Yesterday, 3 or more times")
new_values = c(0,.143,.286,1.5,3)
key = data.frame(old_values = old_values, new_values = new_values)
hmp2_meta[] <- as.data.frame(lapply(hmp2_meta, function(x) key$new_values[match(x, key$old_values)]))

#re-adding diet variable names to columns and adding ID columns
names(hmp2_meta) = diet_vars
hmp2_meta$IDs = hmp2_metadata$`External ID`
hmp2_meta$Participant_ID = hmp2_metadata$`Participant ID`

#making a column that bins the week numbers and a column for diagnosis.
hmp2_metadata = hmp2_metadata %>% mutate(bin = case_when(week_num >= 0 & week_num <= 7 ~ "(0, 7)",
                                                                      week_num >= 8 & week_num <= 15 ~ "(8, 15)",
                                                                      week_num >= 16 & week_num <= 24 ~ "(16, 24)",
                                                                      week_num >= 25 & week_num <= 33 ~ "(25, 33)",
                                                                      week_num >= 34 & week_num <= 42 ~ "(34, 42)",
                                                                      week_num >= 43 & week_num <= 57 ~ "(43, 57)"))

hmp2_meta$bin = hmp2_metadata$bin[match(hmp2_meta$IDs,hmp2_metadata$`External ID`)]
hmp2_meta$diagnosis = disease

hmp2_meta = hmp2_meta[complete.cases(hmp2_meta), ]
hmp2_meta = hmp2_meta %>% drop_na()

#sort by IDs
hmp2_meta = hmp2_meta[order(hmp2_meta$IDs), ]
#set rownames
rownames(hmp2_meta) = hmp2_meta$IDs
hmp2_meta$IDs = NULL

#incorporate additional metadata by merging two data frames
hmp2_meta$External_ID = rownames(hmp2_meta)
hmp2_meta = left_join(hmp2_meta, hmp2_extra_meta, by = "External_ID")

#selecting for adults
hmp2_meta = hmp2_meta[hmp2_meta$consent_age>=18,]

rownames(hmp2_meta) = hmp2_meta$External_ID
hmp2_meta$External_ID = NULL

#save metadata file (includes all phenotypes)
write.csv(hmp2_meta, paste0(output_dir,"/hmp2_meta_allPhenotypes.csv"))

#separating based on disease
hmp2_meta$Sample_ID = rownames(hmp2_meta)
df_CD = hmp2_meta[hmp2_meta$diagnosis == 'CD',]
df_UC = hmp2_meta[hmp2_meta$diagnosis == 'UC',]
df_nonIBD = hmp2_meta[hmp2_meta$diagnosis == 'nonIBD',]
df_list = list(df_CD, df_UC, df_nonIBD)
df_names = list('df_CD', 'df_UC', 'df_nonIBD')
hmp2_meta$Sample_ID = NULL
for (i in seq(df_list)) {
  print(i)
  df = df_list[[i]]
  df_name = df_names[i]
  rownames(df) = df$Sample_ID
  df$Sample_ID = NULL
  write.csv(df,paste0(output_dir,sprintf("/%s.csv",df_name)))
}

#format dataframe of taxonomic profiles
hmp2_species$NCBI_tax_id = NULL
names(hmp2_species) = gsub("_profile","",names(hmp2_species))
hmp2_species$clade_name = gsub(".*s__","",hmp2_species$clade_name)
hmp2_species$clade_name = gsub("_"," ",hmp2_species$clade_name)
rownames(hmp2_species) = hmp2_species$clade_name
bugnames = hmp2_species$clade_name

##filtering for prevalence and abundance: must be at least 10^-5 in at least 25% samples
species_filtered = hmp2_species[rowSums(hmp2_species[,2:dim(hmp2_species)[2]] > 0.00001) >= dim(hmp2_species)[2]*.25, ]
species_filtered = data.frame(species_filtered)
rownames(species_filtered) = species_filtered$clade_name
species_filtered$clade_name = NULL

#transposing species df; add person ID and average be person within bins
species_filtered = as.data.frame(t(species_filtered))
species_filtered = species_filtered[rownames(species_filtered) %in% rownames(hmp2_meta),]
hmp2_meta$samples = rownames(hmp2_meta)
hmp2_meta = hmp2_meta[rownames(hmp2_meta) %in% rownames(species_filtered),]
rownames(hmp2_meta) = hmp2_meta$samples
hmp2_meta$samples = NULL
species_filtered = species_filtered/100
species_filtered$bin = hmp2_meta$bin[match(rownames(species_filtered),rownames(hmp2_meta))]
species_filtered$Participant_ID = hmp2_meta$Participant_ID[match(rownames(species_filtered),rownames(hmp2_meta))]

write.csv(species_filtered, paste0(output_dir,"/hmp_species_df.csv"))
