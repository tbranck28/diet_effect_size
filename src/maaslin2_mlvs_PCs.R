#This script performs feature-wise testing to calculate the effect size of individual diet features on individual
#microbiome features.
library(Maaslin2)
library(ggplot2)
library(grid)
library(vegan)
library(RNOmni)
library(tidyverse)
library(lme4)
library(factoextra)
library(ade4)
library(FactoMineR)
library(dplyr)
library(Rcpp)
library(reshape2)
library(gridExtra)
library(ggpubr)

#set directory paths and load taxonomy and dietary profiles, load diet_variables list for easy grepping when
#separating non-diet metadata from dietary metadata
input_dir = file.path("/home","tobynb","diet_test","input","HPFS","Initial",fsep = "/")
output_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")

bugs <- data.table::fread(paste0(output_dir,"/mlvs_species.csv"), header=T, check.names=F)
input_file <- data.table::fread(paste0(output_dir,"/mlvs_diet.csv"), header=T, check.names=F)

colnames(bugs)[1] = "Sample_ID"
colnames(input_file)[1] = "Sample_ID"  

#ensures both matrices include only the samples present in the other matrix
bugs = bugs[bugs$Sample_ID %in% input_file$Sample_ID,]
input_file = input_file[input_file$Sample_ID %in% bugs$Sample_ID,]

#dataframe formatting
bugs <- as.data.frame(bugs)
input_file <- as.data.frame(input_file)
row.names(bugs) <- bugs$Sample_ID
row.names(input_file) <- input_file$Sample_ID

#ensures samples in the taxonomy matrix have a corresponding dietary profile, and vice versa (and that the samples are in the same order)
setdiff(row.names(bugs), row.names(input_file))
setdiff(row.names(input_file), row.names(bugs))
identical(row.names(bugs), row.names(input_file))
input_file = input_file[match(row.names(bugs), row.names(input_file)), ]
identical(row.names(bugs), row.names(input_file))

#keep only diet variables
input_file_diet = input_file[,2:23]

#PCA to derive dietary patterns
ncomp=2
pca.western.prudent <- psych::principal(r=input_file_diet,
                                        residuals = TRUE,
                                        missing = TRUE,
                                        impute = 'median',
                                        rotate = 'varimax',
                                        nfactors = ncomp,
                                        scores = TRUE)
#put loadings into a table
loadings_table = as.table(pca.western.prudent$loadings)

#visualizes loadings and scree plot
plot(pca.western.prudent$loadings[,1:2],   # x and y data
     pch=21,              # point shape
     bg="black",          # point color
     cex=0,               # point size
     main="MLVS PCA Loadings")
text(pca.western.prudent$loadings[,1:2],             # sets position of labels
     labels=rownames(pca.western.prudent$loadings[,1:2]),cex=.7 )

scree(input_file_diet) #displays factors _and_ components

#grabs factor scores
factors = pca.western.prudent$scores

#removes non-diet metadata
bugs$Sample_ID = NULL
bugs$Participant_ID = NULL
bugs$WeekNum = NULL
bugs$SampNum = NULL
input_file$Sample_ID = NULL

#ensures samples in the factor matrix correspond to those in the dietary profiles, and vice versa (and that the samples are in the same order)
setdiff(row.names(factors), row.names(input_file))
setdiff(row.names(input_file), row.names(factors))
identical(row.names(input_file), row.names(factors))

#set up Maaslin2 input using factor scores + non-diet metadata. Below is code for running maaslin2 when
#the factor scores (dietary patterns) are coded as either a continuous or discrete variable. Note: manuscript reports the continuous version in main figure. Discrete version is in Supplement.
metadata_input = cbind(factors, input_file[, 23:27])
metadata_input_discretized = metadata_input
metadata_input_discretized$RC1 = arules::discretize(metadata_input_discretized$RC1,breaks=3)
metadata_input_discretized$RC2 = arules::discretize(metadata_input_discretized$RC2,breaks=3)

Maaslin2(input_data = bugs, input_metadata = metadata_input, output = paste0(output_dir,"/mlvs_maaslin2_PCs"), random_effects = c("Participant_ID"),fixed_effects = c('RC1','RC2','consent_age','antibiotic','SampNum','WeekNum'))
Maaslin2(input_data = bugs, input_metadata = metadata_input_discretized, output = paste0(output_dir,"/mlvs_maaslin2_PCs_discrete"), random_effects = c("Participant_ID"),fixed_effects = c('RC1','RC2','consent_age','antibiotic','SampNum','WeekNum'))

#reference = c('RC1,[-2.9,-0.516)','RC2,[-3.28,-0.463)'))

#saving a combined file of the factor scores, non-diet metadata, and bugs
mlvs_combined_df = cbind(metadata_input,bugs)
write.csv(mlvs_combined_df,paste0(output_dir,"/mlvs_combined_df.csv"))

############Extra code ###############################
#maaslin2 raw (inputting diet variables as their raw values, not in dietary pattern form)
input_file$SampNum = NULL
input_file$WeekNum = NULL
bugs$SampNum = NULL
bugs$WeekNum = NULL
bugs$Participant_ID = NULL
Maaslin2(input_data = bugs, input_metadata = input_file, output = paste0(output_dir,"/mlvs_maaslin2_raw"), random_effects = c("Participant_ID"))