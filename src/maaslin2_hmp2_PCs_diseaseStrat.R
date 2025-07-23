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

#set directory paths and load taxonomy and dietary profiles
input_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")
output_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")

bugs <- data.table::fread(paste0(input_dir,"/hmp_species_df.csv"), header=T, check.names=F)
input_file <- data.table::fread(paste0(input_dir,"/hmp2_meta_allPhenotypes.csv"), header=T, check.names=F)

#save week bin values as a variable
week_var = bugs[,"bin"]

#data frame formatting
colnames(bugs)[1] = "Sample_ID"  
colnames(input_file)[1] = "Sample_ID"
rownames(week_var) = bugs$Sample_ID

bugs = arrange(bugs, desc(Sample_ID))
input_file = arrange(input_file, desc(Sample_ID))

#ensures both matrices include only the samples present in the other matrix
bugs = bugs[bugs$Sample_ID %in% input_file$Sample_ID,]
input_file = input_file[input_file$Sample_ID %in% bugs$Sample_ID,]

#data frame formatting
bugs <- as.data.frame(bugs)
input_file <- as.data.frame(input_file)
row.names(bugs) <- bugs$Sample_ID
row.names(input_file) <- input_file$Sample_ID

#ensures samples in the taxonomy matrix have a corresponding dietary profile, and vice versa (and that the samples are in the same order)
bugs = bugs[row.names(bugs) %in% row.names(input_file), ]
identical(row.names(bugs), row.names(input_file))
input_file = input_file[match(row.names(bugs), row.names(input_file)), ]
identical(row.names(bugs), row.names(input_file))

#keep only diet variables
input_file_diet = input_file[, 2:20]

#manually (yuck) setting abbreviated column names
names(input_file_diet)[15] = "Red Meat"
names(input_file_diet)[16] = "White meat"
names(input_file_diet)[9] = "Vegetables"
names(input_file_diet)[12] = "Starch (bread, rice, etc.)"
names(input_file_diet)[19] = "Sweets"
names(input_file_diet)[8] = "Fruits (no juice)"
names(input_file_diet)[7] = "Dairy"
names(input_file_diet)[1] = "Soft drinks, tea or coffee with sugar"
names(input_file_diet)[6] = "Yogurt"
names(input_file_diet)[11] = "Whole grains"
names(input_file_diet)[4] = "Fruit juice"
names(input_file_diet)[18] = "Fish"
names(input_file_diet)[3] = "Diet drinks"
names(input_file_diet)[5] = "Alcohol"
names(input_file_diet)[10] = "Beans"
names(input_file_diet)[14] = "Processed meat"
names(input_file_diet)[17] = "Shellfish"

#PCA to derive dietary patterns
ncomp=2
pca.western.prudent <- psych::principal(r=input_file_diet,
                                        residuals = TRUE,
                                        missing = TRUE,
                                        impute = 'median',
                                        rotate = 'varimax',
                                        nfactors = ncomp,
                                        scores = TRUE)

#visualizes loadings and scree plot
plot(pca.western.prudent$loadings[,1:2],   # x and y data
     pch=21,              # point shape
     bg="black",          # point color
     cex=0,               # point size
     main="HMP2 PCA Loadings")
text(pca.western.prudent$loadings[,1:2],             # sets position of labels
     labels=rownames(pca.western.prudent$loadings[,1:2]),cex=.7 )
scree(input_file_diet) #displays factors _and_ components

#grabs factor scores
factors = pca.western.prudent$scores

factors_df = as.data.frame(factors)
bugs_PCs_plotting = bugs 
bugs_PCs_plotting$diagnosis = input_file$diagnosis[match(bugs_PCs_plotting$Sample_ID,input_file$Sample_ID)]
bugs_PCs_plotting$RC1 = factors_df$RC1[match(bugs_PCs_plotting$Sample_ID,rownames(factors_df))]
bugs_PCs_plotting$RC2 = factors_df$RC2[match(bugs_PCs_plotting$Sample_ID,rownames(factors_df))]
write.csv(bugs_PCs_plotting, paste0(output_dir,"/HMP2_bugs_PCs.csv"), row.names = TRUE)

#removes non-diet metadata and formats 
bugs$Sample_ID = NULL
bugs$Participant_ID = NULL
bugs$bin = NULL
input_file$Sample_ID = NULL
input_file$bin = NULL

#format for maaslin, specifically for creating reference variables
input_file$diagnosis = gsub('nonIBD', 'a_nonIBD', input_file$diagnosis)
input_file$diagnosis2 = gsub('nonIBD.non_dysbiosis', 'a_nonIBD.non_dysbiosis', input_file$diagnosis2)

#more formatting for maaslin, removing special characters in string
colnames(bugs) = sub(".*s__", "", colnames(bugs))

#double-checking that both matrices include only the samples present in the other matrix
setdiff(row.names(factors), row.names(input_file))
setdiff(row.names(input_file), row.names(factors))

#set up Maaslin2 input using factor scores + non-diet metadata. Below is code for running maaslin2 when
#the factor scores (dietary patterns) are coded as either a continuous or discrete variable. Note: manuscript reports the continuous version in main figure. Discrete version is in Supplement.
metadata_input = cbind(factors, input_file[, 22:26])
metadata_input$WeekBin = week_var$bin[match(rownames(metadata_input),rownames(week_var))]

#Split metadata_input data frame by diagnosis
X<-split(metadata_input, metadata_input$diagnosis)

#Run Maaslin2 for subset of samples for each diagnosis and when diet factors are continious and discrete
for (i in X) {
  i = as.data.frame(i)
  name = i$diagnosis[1]
  print(name)
  i$diagnosis = NULL
  i$diagnosis2 = gsub("^.*\\.", "", i$diagnosis2)
  i$WeekBin = factor(i$WeekBin)
  i_discrete = i
  i_discrete$RC1 = arules::discretize(i_discrete$RC1,breaks=3,labels = c("low","mid","high"))
  i_discrete$RC2 = arules::discretize(i_discrete$RC2,breaks=3,labels = c("low","mid","high"))
  Maaslin2(input_data = bugs, input_metadata = i, output = paste0(output_dir, sprintf("/hmp2_maaslin2_PCs_diseaseStrat_%s", name)), fixed_effects = c('RC1','RC2','consent_age','antibiotic','diagnosis2','WeekBin'), random_effects = c("Participant_ID"),max_significance = 1.0)
  Maaslin2(input_data = bugs, input_metadata = i_discrete, output = paste0(output_dir, sprintf("/hmp2_maaslin2_PCs_DISCRETE_diseaseStrat_%s", name)), fixed_effects = c('RC1','RC2','consent_age','antibiotic','diagnosis2','WeekBin'), random_effects = c("Participant_ID"),max_significance = 1.0)
  }

#### Run Maaslin2 for raw variables ###################################################################
#maaslin2
names(input_file) <- gsub(x = names(input_file), pattern = "\\(", replacement = "")  
names(input_file) <- gsub(x = names(input_file), pattern = "\\)", replacement = "")  
names(input_file) <- gsub(x = names(input_file), pattern = "\\.", replacement = "")  
names(input_file) <- gsub(x = names(input_file), pattern = "\\,", replacement = "")  
names(input_file) <- gsub(x = names(input_file), pattern = " ", replacement = "_")  
input_maaslin = input_file
input_maaslin$Participant_ID = NULL
Maaslin2(input_data = bugs, input_metadata = input_maaslin, output = paste0(output_dir, "/hmp2_maaslin2_raw"), random_effects = input_file$Participant_ID, reference = c("diagnosis,a_nonIBD","diagnosis2,a_nonIBD.non_dysbiosis"))
