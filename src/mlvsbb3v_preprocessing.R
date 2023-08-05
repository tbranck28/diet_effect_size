setwd("/Users/tobynbranck/Documents/Effect_size/data/")
library("RNOmni")
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
library("grid", lib.loc="/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
library("vegan", lib.loc="/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
library(tidyverse)
library(lme4)
library(factoextra)
library(ade4)
library(FactoMineR)

test = data.table::fread("MLVS/metaphlan/metaphlan_taxonomic_profiles.tsv", header=T, check.names = F)
names(test) <- gsub("# ", "", names(test), fixed = TRUE)
test = test[- grep("Archaea", test$taxonomy),]
test = test[- grep("Eukaryota", test$taxonomy),]
test = test[grep("s__", test$taxonomy),]
names(test) = gsub(pattern = "_dna.*", replacement = "", x = names(test))
test$taxonomy <- gsub('^.*?s__', '', test$taxonomy)
test$taxonomy = gsub("_", " ", test$taxonomy)
#plot to observe distribution of sum of abundance by sample - for samples with < 100, how to handle? why are some so low (in 80-98 range)??
qplot(colSums(test[,-1]), geom="histogram")

#divide by 100 to put on 0-1 scale
rownames(test) = test$taxonomy
test$taxonomy = NULL
bugnames = rownames(test)
test = test/100
rownames(test) = bugnames

#less strict filter for prevalence / abundance; 10^-5 in at least 25% of the samples
test$taxonomy = bugnames
filttest = test[rowSums(test[,-930] > 0.00001) >= dim(test)[2]*.25, ]
filttest = data.frame(filttest)
rownames(filttest) = filttest$taxonomy
filttest$taxonomy = NULL
names(filttest) <- sub("^X", "", names(filttest))
column_names = colnames(filttest)
column_names = paste(substr(column_names, 1, 6), substr(column_names, 9, nchar(column_names)), sep='')
colnames(filttest) = column_names
colnames(filttest) <- gsub("SF05", "s1", colnames(filttest))
colnames(filttest) <- gsub("SF06", "s2", colnames(filttest))
colnames(filttest) <- gsub("SF07", "s3", colnames(filttest))
colnames(filttest) <- gsub("SF08", "s4", colnames(filttest))
filttest = as.data.frame(t(filttest))
samplist = rownames(filttest)
samplist <- gsub(".*_s1.*", "1", samplist)
samplist <- gsub(".*_s2.*", "2", samplist)
samplist <- gsub(".*_s3.*", "3", samplist)
samplist <- gsub(".*_s4.*", "4", samplist)
weeklist <- gsub(".*2.*", "1", samplist)
weeklist <- gsub(".*3.*", "2", weeklist)
weeklist <- gsub(".*4.*", "2", weeklist)
filttest$SampNum = samplist
filttest$WeekNum = weeklist
filttest$Participant_ID = str_extract(rownames(filttest), "^.*(?=(_))")
write.table(filttest, "/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_species.csv", sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

#stricter filter for prevalence / abundance; 10^-5 in at least 60% of the samples --> NOTE NOW USING SAME AS LESS STRICT FILTER (25%)
test$taxonomy = bugnames
filttest = test[rowSums(test[,-930] > 0.00001) >= dim(test)[2]*.25, ]
filttest = data.frame(filttest)
rownames(filttest) = filttest$taxonomy
filttest$taxonomy = NULL
names(filttest) <- sub("^X", "", names(filttest))
column_names = colnames(filttest)
column_names = paste(substr(column_names, 1, 6), substr(column_names, 9, nchar(column_names)), sep='')
colnames(filttest) = column_names
colnames(filttest) <- gsub("SF05", "s1", colnames(filttest))
colnames(filttest) <- gsub("SF06", "s2", colnames(filttest))
colnames(filttest) <- gsub("SF07", "s3", colnames(filttest))
colnames(filttest) <- gsub("SF08", "s4", colnames(filttest))
filttest = as.data.frame(t(filttest))
samplist = rownames(filttest)
samplist <- gsub(".*_s1.*", "1", samplist)
samplist <- gsub(".*_s2.*", "2", samplist)
samplist <- gsub(".*_s3.*", "3", samplist)
samplist <- gsub(".*_s4.*", "4", samplist)
weeklist <- gsub(".*2.*", "1", samplist)
weeklist <- gsub(".*3.*", "2", weeklist)
weeklist <- gsub(".*4.*", "2", weeklist)
filttest$SampNum = samplist
filttest$WeekNum = weeklist
filttest$Participant_ID = str_extract(rownames(filttest), "^.*(?=(_))")
write.table(filttest, "/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_species_strict_filtering.csv", sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)


#importing diet
diet_names = data.table::fread("/Users/tobynbranck/Documents/STARR_diet/maaslin2_EnergyAdj/7ddrEA_wkAvg_CovsInv.txt", header=T, check.names = F)
diet_names = as.data.frame(diet_names)
rownames(diet_names) = diet_names$SampleID
diet_names$SampleID = NULL
diet_names = t(diet_names)
diet_names = as.data.frame(diet_names)
diet_names$consent_age = as.character(diet_names$agemlvs)
diet_names$agemlvs = NULL
diet_names$antibiotic <- as.character(diet_names$abx)
diet_names$abx = NULL
diet_names$Participant_ID <- as.character(diet_names$PersonID)
diet_names$PersonID = NULL

#standardizing diet variables
A = function(x) scale(x)
##### or instead use inverse normal transformation ##########
##### A <- function(x) RankNorm(as.numeric(na.omit(x))) #####

#scaled.diet_names <- cbind(lapply(diet_names[,1:22], A), diet_names[,23:25])
scaled.diet_names = diet_names #delete this line if need to add back standardization
scaled.diet_names = scaled.diet_names[!rowSums(is.na(scaled.diet_names)) > 0, ]
scaled.diet_names = as.data.frame(scaled.diet_names)
samplist_diet = rownames(scaled.diet_names)
samplist_diet <- gsub(".*_s1.*", "1", samplist_diet)
samplist_diet <- gsub(".*_s2.*", "2", samplist_diet)
samplist_diet <- gsub(".*_s3.*", "3", samplist_diet)
samplist_diet <- gsub(".*_s4.*", "4", samplist_diet)
weeklist_diet <- gsub(".*2.*", "1", samplist_diet)
weeklist_diet <- gsub(".*3.*", "2", weeklist_diet)
weeklist_diet <- gsub(".*4.*", "2", weeklist_diet)
scaled.diet_names$SampNum = samplist_diet
scaled.diet_names$WeekNum = weeklist_diet
write.table(scaled.diet_names, "/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_diet.csv", sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)
