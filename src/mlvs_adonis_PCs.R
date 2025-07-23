#Calculates PERMANOVA for the MLVS data; dietary patterns are derived from the dietary profiles
library(ggplot2)
library(grid)
library(vegan)
library(tidyverse)
library(lme4)
library(cowplot)
library(psych)

#set directories and load dietary profiles
input_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")
output_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")
input_files <- c(paste0(input_dir,"/mlvs_diet.csv"))

#initialize vectors which will store the adonis result
adonis_rsq = vector()
adonis_pval = vector()

#loops through input matrices, in the case of MLVS, we only have 1 for this cohort.
for (i in input_files) {
  #loads microbiome and diet matrices
  bugs <- data.table::fread(paste0(input_dir,"/mlvs_species.csv"), header=T, check.names=F)
  input_file <- data.table::fread(i, header=T, check.names=F)
  
  #formats both matrices consistently
  colnames(input_file)[1] = "Sample_ID"  
  input_file = arrange(input_file, Sample_ID)
  colnames(bugs)[1] = "Sample_ID"  
  bugs = arrange(bugs, Sample_ID)
  
  #ensures both matrices include only the samples present in the other matrix
  bugs = bugs[bugs$Sample_ID %in% input_file$Sample_ID,]
  input_file = input_file[input_file$Sample_ID %in% bugs$Sample_ID,]
  
  #more dataframe formatting
  bugs <- as.data.frame(bugs)
  input_file <- as.data.frame(input_file)
  row.names(bugs) <- bugs$Sample_ID
  row.names(input_file) <- input_file$Sample_ID
  
  #ensures samples in the taxonomy matrix have a corresponding dietary profile, and vice versa
  setdiff(row.names(bugs), row.names(input_file))
  setdiff(row.names(input_file), row.names(bugs))
  identical(row.names(bugs), row.names(input_file))
  
  bugs$Sample_ID = NULL
  bugs$Participant_ID = NULL
  bugs$bin = NULL
  bugs$SampNum = NULL
  bugs$WeekNum = NULL
  
  #### running adonis for individual diet features
  adonis_res_rsq_individual = vector()
  adonis_res_pval_individual = vector()
  for (col in colnames(input_file)[2:23]){
    print(col)
    diet_var = input_file[,col]
    adonis.univ = adonis2(bugs ~ antibiotic + consent_age + diet_var, data = input_file, method = "bray", strata = input_file$Participant_ID,by="margin")
    adonis_res_rsq_individual[col] = adonis.univ[3,]$R2
    adonis_res_pval_individual[col] = adonis.univ[3,]$`Pr(>F)`
    study = "MLVS"
  }
  adonis_df = as.data.frame(cbind(adonis_res_rsq_individual, adonis_res_pval_individual,study))
  write.csv(adonis_df, paste0(output_dir,"/mlvs_individual_diet_adonis_results.csv"), row.names = TRUE)
  
  #removing non-diet metadata
  input_file_diet = input_file[, 2:23]
  input_file$Sample_ID = NULL

  adonis_label = vector()
  ncomp=2
  #derives dietary patterns
  pca.western.prudent <- psych::principal(r=input_file_diet,
                                          residuals = TRUE,
                                          missing = TRUE,
                                          impute = 'median',
                                          rotate = 'varimax',
                                          nfactors = ncomp,
                                          scores = TRUE)
  
  #pulls dietary patterns from "psych::principal" output
  factors = as.data.frame(pca.western.prudent$scores)
  
  #performs permanova using dietary patterns + other covariates
  adonis.PCs = adonis2(bugs ~ antibiotic + consent_age + factors$RC1 + factors$RC2, data = input_file, method = "bray", strata = input_file$Participant_ID,by="margin")
  adonis_rsq = adonis.PCs[3:4,]$R2
  adonis_pval = adonis.PCs[3:4,]$`Pr(>F)`
  adonis_label = gsub('.*\\$', '', rownames(adonis.PCs[3:4,]))
}

#formats and saves adonis output which will be used in downstream visualizations that combine effect sizes
#from different cohorts
adonis_result = rbind(adonis_label, adonis_rsq, adonis_pval)
adonis_result = as.data.frame(t(adonis_result))
adonis_result$adonis_pval = as.numeric(adonis_result$adonis_pval)
adonis_result$adonis_rsq = as.numeric(adonis_result$adonis_rsq)
adonis_result$p_adj = p.adjust(adonis_result$adonis_pval, "fdr")
adonis_result$stars = cut(adonis_result$p_adj, c(0, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "-", ""))

write.csv(adonis_result, paste0(output_dir,"/mlvs_adonisPCs_results.csv"), row.names = FALSE)
loads = as.data.frame(unclass(pca.western.prudent$loadings))
write.csv(loads, paste0(output_dir,"/mlvs_combined_loadings.csv"), row.names = TRUE)

