setwd("/Users/tobynbranck/Documents/Effect_size/analysis/")
library("Maaslin2", lib.loc="/Library/Frameworks/R.framework/Versions/4.1/Resources/library")
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
library("grid", lib.loc="/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
library("vegan", lib.loc="/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
library("RNOmni", lib.loc="/Library/Frameworks/R.framework/Versions/4.1/Resources/library")
library(tidyverse)
library(lme4)
library(factoextra)
library(ade4)
library(FactoMineR)
library(dplyr)
library(imager)
library(Rcpp)
library(reshape2)
library(gridExtra)
library(ggpubr)

bugs <- data.table::fread("mlvs_species.csv", header=T, check.names=F)
colnames(bugs)[1] = "Sample_ID"  
arrange(bugs, desc(Sample_ID))

grep_list = read.csv("mlvs_diet_vars_list.csv", header=FALSE)$V1
input_file <- data.table::fread("mlvs_diet.csv", header=T, check.names=F)

colnames(input_file)[1] = "Sample_ID"  
arrange(input_file, desc(Sample_ID))

bugs = bugs[bugs$Sample_ID %in% input_file$Sample_ID,]
input_file = input_file[input_file$Sample_ID %in% bugs$Sample_ID,]

bugs <- as.data.frame(bugs)
input_file <- as.data.frame(input_file)
row.names(bugs) <- bugs$Sample_ID
row.names(input_file) <- input_file$Sample_ID

setdiff(row.names(bugs), row.names(input_file))
setdiff(row.names(input_file), row.names(bugs))

input_file <- na.omit(input_file)

bugs = bugs[row.names(bugs) %in% row.names(input_file), ]
identical(row.names(bugs), row.names(input_file))
input_file = input_file[match(row.names(bugs), row.names(input_file)), ]
identical(row.names(bugs), row.names(input_file))

str(input_file)
input_file_diet = input_file[, grep_list]

#PCA
ncomp=2
pca.western.prudent <- psych::principal(r=input_file_diet,
                                        residuals = TRUE,
                                        missing = TRUE,
                                        impute = 'median',
                                        rotate = 'varimax',
                                        nfactors = ncomp,
                                        scores = TRUE)
loadings_table = as.table(pca.western.prudent$loadings)

plot(pca.western.prudent$loadings[,1:2],   # x and y data
     pch=21,              # point shape
     bg="black",          # point color
     cex=0,               # point size
     main="MLVS PCA Loadings")
text(pca.western.prudent$loadings[,1:2],             # sets position of labels
     labels=rownames(pca.western.prudent$loadings[,1:2]),cex=.7 )

scree(input_file_diet) #displays factors _and_ components
factors = pca.western.prudent$scores

bugs$Sample_ID = NULL
bugs$Participant_ID = NULL
bugs$WeekNum = NULL
bugs$SampNum = NULL
input_file$Sample_ID = NULL

setdiff(row.names(factors), row.names(input_file))
setdiff(row.names(input_file), row.names(factors))


#maaslin2 factor scores
metadata_input = cbind(factors, input_file[, 23:27])
metadata_input_discretized = metadata_input
metadata_input_discretized$RC1 = discretize(metadata_input_discretized$RC1,breaks=3)
metadata_input_discretized$RC2 = discretize(metadata_input_discretized$RC2,breaks=3)

Maaslin2(input_data = bugs, input_metadata = metadata_input, output = "mlvs_maaslin2_PCs", random_effects = c("Participant_ID"))
Maaslin2(input_data = bugs, input_metadata = metadata_input_discretized, output = "mlvs_maaslin2_PCs_discrete", random_effects = c("Participant_ID"),reference = c('RC1,[-2.9,-0.516)','RC2,[-3.28,-0.463)'))

mlvs_combined_df = cbind(metadata_input,bugs)
write.csv(mlvs_combined_df,"mlvs_combined_df.csv")

############Extra code ###############################
#maaslin2 raw (no standardization of diet variables)
input_file$SampNum = NULL
input_file$WeekNum = NULL
bugs$SampNum = NULL
bugs$WeekNum = NULL
bugs$Participant_ID = NULL
Maaslin2(input_data = bugs, input_metadata = input_file, output = "mlvs_maaslin2_raw", random_effects = c("Participant_ID"))

#boxplot for MLVS maaslin effect sizes - factor scores
# mlvs_factor_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_maaslin2_PCs/all_results.tsv", header=T, check.names=T))
# 
# filttest = apply(mlvs_factor_results, 1, function(r) any(r %in% c("antibiotic", "consent_age")))
# test_covariates = mlvs_factor_results[filttest,]
# test_microbiome = mlvs_factor_results[!filttest,]
# 
# mlvs_factor = ggplot(mlvs_factor_results,aes(x=reorder(metadata,-abs(coef)), y=abs(coef))) +
#   geom_boxplot() +
#   labs(y=expression(beta~'coefficients'), x="") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#   theme(axis.text.x = element_text(size = 9)) +
#   theme(axis.text.y = element_text(size = 12)) + 
#   labs(title = "MLVS", subtitle = "Maaslin2 Effect Sizes - Factor Scores") +
#   scale_y_continuous(limits = c(0, 3))

#boxplot for mlvs maaslin effect sizes - raw diet
#mlvs_raw_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_maaslin2_raw/all_results.tsv", header=T, check.names=T))

# filttest = apply(mlvs_raw_results, 1, function(r) any(r %in% c("antibiotic", "consent_age")))
# test_covariates = mlvs_raw_results[filttest,]
# test_microbiome = mlvs_raw_results[!filttest,]
# 
# mlvs_raw = ggplot(mlvs_raw_results,aes(x=reorder(metadata,-abs(coef)), y=abs(coef))) +
#   geom_boxplot() +
#   labs(y=expression(beta~'coefficients'), x="") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#   theme(axis.text.x = element_text(size = 7)) +
#   theme(axis.text.y = element_text(size = 12)) + 
#   labs(title = "MLVS", subtitle = "Maaslin2 Effect Sizes - Raw Diet") +
#   scale_y_continuous(limits = c(0, 3))

###########################################################
#Supplemental: factor scores vs. raw mlvs and hmp2 maaslin#
###########################################################

# boxplot for hmp2 maaslin effect sizes - factor scores
# hmp2_factor_results_CD = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_diseaseStrat_CD/all_results.tsv", header=T, check.names=T))
# hmp2_factor_results_CD$metadata[hmp2_factor_results_CD$metadata == 'diagnosis2'] <- 'dysbiosis'
# 
# hmp2_factor_results_UC = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_diseaseStrat_UC/all_results.tsv", header=T, check.names=T))
# hmp2_factor_results_UC$metadata[hmp2_factor_results_UC$metadata == 'diagnosis2'] <- 'dysbiosis'
# 
# hmp2_factor_results_nonIBD = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_diseaseStrat_a_nonIBD/all_results.tsv", header=T, check.names=T))
# hmp2_factor_results_nonIBD$metadata[hmp2_factor_results_nonIBD$metadata == 'diagnosis2'] <- 'dysbiosis'
# 
# hmp2_factor_results = rbind(hmp2_factor_results_CD,hmp2_factor_results_UC,hmp2_factor_results_nonIBD)
# 
# filttest = apply(hmp2_factor_results, 1, function(r) any(r %in% c("antibiotic", "consent_age", "diagnosis", "dysbiosis")))
# test_covariates = hmp2_factor_results[filttest,]
# test_microbiome = hmp2_factor_results[!filttest,]
# 
# hmp2_factor = ggplot(hmp2_factor_results,aes(x=reorder(metadata,-abs(coef)), y=abs(coef))) +
#   geom_boxplot() +
#   labs(y=expression(beta~'coefficients'), x="") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#   theme(axis.text.x = element_text(size = 9)) +
#   theme(axis.text.y = element_text(size = 12)) + 
#   labs(title = "HMP2", subtitle = "Maaslin2 Effect Sizes - PC Scores")
# 
# #boxplot for hmp2 maaslin effect sizes - raw diet
# hmp2_raw_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_raw/all_results.tsv", header=T, check.names=T))
# hmp2_raw_results$metadata[hmp2_raw_results$metadata == 'diagnosis2'] <- 'dysbiosis'
# 
# filttest = apply(hmp2_raw_results, 1, function(r) any(r %in% c("antibiotic", "consent_age", "diagnosis", "diagnosis2")))
# test_covariates = hmp2_raw_results[filttest,]
# test_microbiome = hmp2_raw_results[!filttest,]
# 
# hmp2_raw = ggplot(hmp2_raw_results,aes(x=reorder(metadata,-abs(coef)), y=abs(coef))) +
#   geom_boxplot() +
#   labs(y=expression(beta~'coefficients'), x="") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#   theme(axis.text.x = element_text(size = 7)) +
#   theme(axis.text.y = element_text(size = 12)) + 
#   labs(title = "HMP2", subtitle = "Maaslin2 Effect Sizes - Raw Diet") +
#   scale_x_discrete(label=function(x) str_trunc(x, 15, "right"))
# 
# ############ combining plots ########################
# png("WesternAdult_PCVsRawBox.png")
# ggarrange(hmp2_factor,hmp2_raw,mlvs_factor,mlvs_raw, labels = c("A","B","C","D"),nrow = 2, ncol=2,align = "hv")
# dev.off()