setwd("/Users/tobynbranck/Documents/Effect_size/analysis/")
library("ggplot2")
library("grid")
library("vegan")
library(tidyverse)
library(lme4)
library("cowplot")
library(psych)
library(MicEco)

input_files <- c("mlvs_diet.csv")

names = vector()
adonis_rsq = vector()
adonis_pval = vector()
results = vector()
datafs = list()
n=0
for (i in input_files) {
  n=n+1
  name = gsub(".*[_]([^.]+)[.].*", "\\1", i)
  cohort = "hmp2"
  bugs <- data.table::fread("mlvs_species.csv", header=T, check.names=F)
  
  input_file <- data.table::fread(i, header=T, check.names=F)
  colnames(input_file)[1] = "Sample_ID"  
  arrange(input_file, desc(Sample_ID))
  
  colnames(bugs)[1] = "Sample_ID"  
  arrange(bugs, desc(Sample_ID))
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
  input_file_diet = input_file[, 2:23]
  adonis_output = vector()
  adonis_res_rsq = vector()
  adonis_res_pval = vector()
  adonis_label = vector()
  
  bugs$Sample_ID = NULL
  bugs$Participant_ID = NULL
  bugs$bin = NULL
  
  bugs$SampNum = NULL
  bugs$WeekNum = NULL
  input_file$Sample_ID = NULL
  ncomp=2
  pca.western.prudent <- psych::principal(r=input_file_diet,
                                          residuals = TRUE,
                                          missing = TRUE,
                                          impute = 'median',
                                          rotate = 'varimax',
                                          nfactors = ncomp,
                                          scores = TRUE)
  
  scree(input_file_diet) #displays factors _and_ components
  factors = as.data.frame(pca.western.prudent$scores)
  
  adonis.PCs = adonis(bugs ~ antibiotic + consent_age + factors$RC1 + factors$RC2, data = input_file, method = "bray", strata = input_file$Participant_ID)
  print (adonis.PCs)
  adonis_rsq = adonis.PCs$aov.tab[3:4,]$R2
  adonis_pval = adonis.PCs$aov.tab[3:4,]$`Pr(>F)`
  adonis_label = gsub('.*\\$', '', rownames(adonis.PCs$aov.tab[3:4,]))
}
univar_res_rgnav = rbind(adonis_label, adonis_rsq, adonis_pval)
univar_res_rgnav = as.data.frame(t(univar_res_rgnav))
univar_res_rgnav$adonis_pval = as.numeric(univar_res_rgnav$adonis_pval)
univar_res_rgnav$adonis_rsq = as.numeric(univar_res_rgnav$adonis_rsq)
univar_res_rgnav$p_adj = p.adjust(univar_res_rgnav$adonis_pval, "fdr")
univar_res_rgnav$stars = cut(univar_res_rgnav$p_adj, c(0, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "-", ""))

adonis_OmegaSq(adonis.PCs,partial = TRUE)

write.csv(univar_res_rgnav, "/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_adonisPCs_results.csv", row.names = FALSE)
loads = as.data.frame(unclass(pca.western.prudent$loadings))
write.csv(loads, "/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_combined_loadings.csv", row.names = TRUE)