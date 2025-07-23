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

#load metaphlan table, grep for species-level bacterial assignments, format column names and taxonomy names
mpa = data.table::fread(paste0(input_dir,"/HPFS/Initial/metaphlan_taxonomic_profiles.tsv"), header=T, check.names = F)
names(mpa) <- gsub("# ", "", names(mpa), fixed = TRUE)
mpa = mpa[- grep("Archaea", mpa$taxonomy),]
mpa = mpa[- grep("Eukaryota", mpa$taxonomy),]
mpa = mpa[grep("s__", mpa$taxonomy),]
names(mpa) = gsub(pattern = "_dna.*", replacement = "", x = names(mpa))
mpa$taxonomy <- gsub('^.*?s__', '', mpa$taxonomy)
mpa$taxonomy = gsub("_", " ", mpa$taxonomy)

#divide by 100 to put on 0-1 scale, preserving microbe names
bugnames = mpa$taxonomy
mpa = mpa[,2:930]/100
rownames(mpa) = bugnames

#Filter for prevalence / abundance; 10^-5 in at least 25% of the samples; creating a taxonomy column to preserve microbe ID during filtering
mpa$taxonomy = bugnames
mpa_filtered = mpa[rowSums(mpa[,-930] > 0.00001) >= dim(mpa)[2]*.25, ]
mpa_filtered = data.frame(mpa_filtered)
rownames(mpa_filtered) = mpa_filtered$taxonomy
mpa_filtered$taxonomy = NULL

#formatting sample names in metaphlan dataframe to be consistent with samples names in the diet metadata. SF05/s1 refers to the first sample (and so on) 
#of the 4 samples collected for each individual in the MLVS cohort. Recall that samples 1 and 2 were collected a week apart and then after 6 months,
#samples 3 and 4 were collected a week apart.
names(mpa_filtered) <- sub("^X", "", names(mpa_filtered))
column_names = colnames(mpa_filtered)
column_names = paste(substr(column_names, 1, 6), substr(column_names, 9, nchar(column_names)), sep='')
colnames(mpa_filtered) = column_names
colnames(mpa_filtered) <- gsub("SF05", "s1", colnames(mpa_filtered))
colnames(mpa_filtered) <- gsub("SF06", "s2", colnames(mpa_filtered))
colnames(mpa_filtered) <- gsub("SF07", "s3", colnames(mpa_filtered))
colnames(mpa_filtered) <- gsub("SF08", "s4", colnames(mpa_filtered))
mpa_filtered = as.data.frame(t(mpa_filtered))

#adding sample number (options: 1 through 4), week number (options: 1 or 2), and participant ID to mpa matrix
mpa_filtered = mpa_filtered %>% mutate(SampNum = case_when(str_detect(rownames(.), "_s1") ~ "1",
                                                   str_detect(rownames(.), "_s2") ~ "2",
                                                   str_detect(rownames(.), "_s3") ~ "3",
                                                   str_detect(rownames(.), "_s4") ~ "4"))
mpa_filtered = mpa_filtered %>% mutate(WeekNum = case_when(str_detect(rownames(.), "_s1") ~ "1",
                                                   str_detect(rownames(.), "_s2") ~ "1",
                                                   str_detect(rownames(.), "_s3") ~ "2",
                                                   str_detect(rownames(.), "_s4") ~ "2"))
mpa_filtered$Participant_ID = str_extract(rownames(mpa_filtered), "^.*(?=(_))")
write.table(mpa_filtered, paste0(output_dir,"/mlvs_species.csv"), sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

#stricter filter for prevalence / abundance; 10^-5 in at least 60% of the samples --> NOTE NOW USING SAME AS LESS STRICT FILTER (25%)
# mpa$taxonomy = bugnames
# mpa_filtered = mpa[rowSums(mpa[,-930] > 0.00001) >= dim(mpa)[2]*.25, ]
# mpa_filtered = data.frame(mpa_filtered)
# rownames(mpa_filtered) = mpa_filtered$taxonomy
# mpa_filtered$taxonomy = NULL
# names(mpa_filtered) <- sub("^X", "", names(mpa_filtered))
# column_names = colnames(mpa_filtered)
# column_names = paste(substr(column_names, 1, 6), substr(column_names, 9, nchar(column_names)), sep='')
# colnames(mpa_filtered) = column_names
# colnames(mpa_filtered) <- gsub("SF05", "s1", colnames(mpa_filtered))
# colnames(mpa_filtered) <- gsub("SF06", "s2", colnames(mpa_filtered))
# colnames(mpa_filtered) <- gsub("SF07", "s3", colnames(mpa_filtered))
# colnames(mpa_filtered) <- gsub("SF08", "s4", colnames(mpa_filtered))
# mpa_filtered = as.data.frame(t(mpa_filtered))
# samplist = rownames(mpa_filtered)
# samplist <- gsub(".*_s1.*", "1", samplist)
# samplist <- gsub(".*_s2.*", "2", samplist)
# samplist <- gsub(".*_s3.*", "3", samplist)
# samplist <- gsub(".*_s4.*", "4", samplist)
# weeklist <- gsub(".*2.*", "1", samplist)
# weeklist <- gsub(".*3.*", "2", weeklist)
# weeklist <- gsub(".*4.*", "2", weeklist)
# mpa_filtered$SampNum = samplist
# mpa_filtered$WeekNum = weeklist
# mpa_filtered$Participant_ID = str_extract(rownames(mpa_filtered), "^.*(?=(_))")
# write.table(mpa_filtered, "/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_species_strict_filtering.csv", sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)


#replication test
diet_initial = as.data.frame(data.table::fread(paste0(input_dir,"/HPFS/Initial/all_diet.txt"), header=T, check.names = F))
week1_diet_vars = grep("wk1avg", names(diet_initial), value = TRUE)
week2_diet_vars = grep("wk2avg", names(diet_initial), value = TRUE)
diet_week1 = diet_initial[,week1_diet_vars]
diet_week2 = diet_initial[,week2_diet_vars]
diet_week1$Participant_ID = diet_initial$aliasid
diet_week2$Participant_ID = diet_initial$aliasid
names(diet_week1) = gsub("wk1.*","",names(diet_week1))
diet_week1$WeekNum = 1
names(diet_week2) = gsub("wk2.*","",names(diet_week2))
diet_week2$WeekNum = 2
diet = rbind(diet_week1,diet_week2)
diet = diet[!rowSums(is.na(diet))>20, ]
diet = merge(mpa_filtered[,c("Participant_ID","SampNum","WeekNum")], diet, all.x=TRUE, by = c("Participant_ID","WeekNum"))
diet = diet[rowSums(is.na(diet[,4:dim(diet)[2]])) != ncol(diet[,4:dim(diet)[2]]), ]
diet$sampleID = paste0(diet$Participant_ID,"_s",diet$SampNum)

#replace column names with human-readable diet variable names
diet_names = data.table::fread(paste0(input_dir,"/HPFS/Initial/7ddr_codes_engadj.txt"), header=F, check.names = F)
diet_names_readable = data.table::fread(paste0(input_dir,"/HPFS/Initial/7ddr_diet_names_engadj.txt"), header=F, check.names = F)
diet_names$key = diet_names_readable$V1
diet_names$V1 = gsub("_","",diet_names$V1)
diet_names$V1 = gsub("w1d1","",diet_names$V1)
colnames(diet) <- dplyr::recode(
  colnames(diet), 
  !!!setNames(as.character(diet_names$key), diet_names$V1)
)

#read in file containing covariates and add covariates to diet info to make one combined metadata file
covs = data.table::fread(paste0(input_dir,"/HPFS/Initial/covariates7ddr.txt"), header=T, check.names = F)
diet$consent_age = covs$agemlvs[match(diet$sampleID,covs$sample)]
diet$antibiotic = covs$abx[match(diet$sampleID,covs$sample)]
rownames(diet) = diet$sampleID
diet$sampleID = NULL
diet = diet[order(row.names(diet)), ]
remove = setdiff(colnames(mlvs_diet_new),colnames(mlvs_diet))
diet = diet %>% select(-one_of(remove))
write.table(diet, paste0(output_dir,"/mlvs_diet_new.csv"), sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)



#importing and formatting diet + additional metadata; note we are using the 
diet = data.table::fread(paste0(input_dir,"/HPFS/Initial/7ddrEA_wkAvg_CovsInv.txt"), header=T, check.names = F)
diet = as.data.frame(diet)
rownames(diet) = diet$SampleID
diet$SampleID = NULL
diet = as.data.frame(t(diet))
diet$consent_age = as.character(diet$agemlvs)
diet$agemlvs = NULL
diet$antibiotic <- as.character(diet$abx)
diet$abx = NULL
diet$Participant_ID <- as.character(diet$PersonID)
diet$PersonID = NULL

#removes 2 samples that had all dietary records with all 0s
diet = diet[!rowSums(is.na(diet)) > 0, ]

#adding SampNum and WeekNum columns to the diet metadata as we did above in the metaphlan table
diet = diet %>% mutate(SampNum = case_when(str_detect(rownames(.), "_s1") ~ "1",
                                                           str_detect(rownames(.), "_s2") ~ "2",
                                                           str_detect(rownames(.), "_s3") ~ "3",
                                                           str_detect(rownames(.), "_s4") ~ "4"))
diet = diet %>% mutate(WeekNum = case_when(str_detect(rownames(.), "_s1") ~ "1",
                                                           str_detect(rownames(.), "_s2") ~ "1",
                                                           str_detect(rownames(.), "_s3") ~ "2",
                                                           str_detect(rownames(.), "_s4") ~ "2"))
write.table(diet, paste0(output_dir,"/mlvs_diet.csv"), sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)
