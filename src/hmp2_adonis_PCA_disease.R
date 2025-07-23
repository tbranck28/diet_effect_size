#Calculates PERMANOVA for the HMP2 data; dietary patterns are derived from the dietary profiles
library(ggplot2)
library(grid)
library(vegan)
library(tidyverse)
library(lme4)
library(cowplot)
library(psych)


#set directories and load dietary profiles
input_dir = file.path("/home","tobynb","diet_test","input","HMP2","Intermediate",fsep = "/")
output_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")
files = c("df_CD.csv","df_nonIBD.csv","df_UC.csv")
allPts <- data.table::rbindlist(lapply(files, function(f) data.table::fread(file.path(output_dir, f), header = TRUE, check.names = FALSE)))
cohort = "hmp2"

#load microbiome data
bugs <- data.table::fread(paste0(output_dir,"/hmp_species_df.csv"), header=T, check.names=F)

#initialize vectors which will store the adonis result
adonis_rsq = vector()
adonis_pval = vector()
adonis_label = vector()

#formats both matrices consistently
colnames(allPts)[1] = "Sample_ID"  
allPts = arrange(allPts, desc(Sample_ID))
colnames(bugs)[1] = "Sample_ID"  
bugs = arrange(bugs, desc(Sample_ID))

#ensures both matrices include only the samples present in the other matrix
bugs = bugs[bugs$Sample_ID %in% allPts$Sample_ID,]
allPts = allPts[allPts$Sample_ID %in% bugs$Sample_ID,]

#more dataframe formatting
bugs <- as.data.frame(bugs)
allPts <- as.data.frame(allPts)
row.names(bugs) <- bugs$Sample_ID
row.names(allPts) <- allPts$Sample_ID

#ensures samples in the taxonomy matrix have a corresponding dietary profile, and vice versa
setdiff(row.names(bugs), row.names(allPts))
setdiff(row.names(allPts), row.names(bugs))
identical(row.names(bugs), row.names(allPts))

#removing non-diet metadata
input_file_diet = allPts[, 2:22]
bugs$Sample_ID = NULL
bugs$Participant_ID = NULL
bugs$bin = NULL

# asks whether diet is different across disease phenotypes
diet_adonis = adonis2(input_file_diet ~ diagnosis, data = allPts, method = "manhattan", strata = allPts$Participant_ID)

#derives dietary patterns
ncomp <- 2
pca.western.prudent <- psych::principal(r=input_file_diet,
                                        residuals = TRUE,
                                        missing = TRUE,
                                        impute = 'median',
                                        rotate = 'varimax',
                                        nfactors = ncomp,
                                        scores = TRUE)

scree(input_file_diet) #displays factors _and_ components

#pulls dietary patterns from "psych::principal" output
factors = as.data.frame(pca.western.prudent$scores)

#splits larger matrix containing diet vectors into diagnosis-specific dataframes (UC, CD, and nonIBD)
X<-split(allPts, allPts$diagnosis)

#initialize iteration number
n=0

#creates empty list which we'll add adonis output to sequentially as we loop through diagnosis
datalist = list()
for (i in X) {
  n=n+1
  
  #filter diet patterns table to include only the samples contained in diet input pertaining a particular diagnosis
  factors_df = factors %>% filter(row.names(factors) %in% rownames(i))
  
  #filter taxonomy table to include only the samples contained in diet input pertaining a particular diagnosis
  bugs_df = bugs %>% filter(row.names(bugs) %in% rownames(i))
  
  print(identical(row.names(bugs_df), row.names(factors_df)))
  
  #performs permanova using dietary patterns + other covariates
  adonis.PCs = adonis2(bugs_df ~ antibiotic + consent_age + diagnosis2 + factors_df$RC1 + factors_df$RC2, data = i, method = "bray", strata = i$Participant_ID,by="margin")
  adonis_rsq = adonis.PCs[4:5,]$R2
  adonis_pval = adonis.PCs[4:5,]$`Pr(>F)`
  adonis_label = gsub('.*\\$', '', rownames(adonis.PCs[4:5,]))
  
  univar_res_rgnav = rbind(adonis_label, adonis_rsq, adonis_pval)
  univar_res_rgnav = as.data.frame(t(univar_res_rgnav))
  univar_res_rgnav$adonis_pval = as.numeric(univar_res_rgnav$adonis_pval)
  univar_res_rgnav$adonis_rsq = as.numeric(univar_res_rgnav$adonis_rsq)
  univar_res_rgnav$p_adj = p.adjust(univar_res_rgnav$adonis_pval, "fdr")
  univar_res_rgnav$stars = cut(univar_res_rgnav$p_adj, c(0, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "-", ""))
  datalist[[n]] = univar_res_rgnav
}

df_comb = rbind(datalist[[1]],datalist[[2]],datalist[[3]])

df_comb$state = c("CD","CD","nonIBD","nonIBD","UC","UC")
df_comb$cat = paste(df_comb$state,df_comb$adonis_label)
df_comb$cat = factor(df_comb$cat, levels=c("nonIBD RC2","nonIBD RC1","CD RC2","CD RC1","UC RC2","UC RC1"))
write.csv(df_comb, paste0(output_dir,"/HMP2_adonisPCs_byDisease_results.csv"), row.names = FALSE)
loads = as.data.frame(unclass(pca.western.prudent$loadings))
write.csv(loads, paste0(output_dir,"/hmp2_combined_loadings.csv"), row.names = TRUE)

### running adonis for individual diet features
#initialize empty list that will store each phenotypes adonis results df
datafs = list()
#initialize iteration number
n=0
for (i in files) {
  n=n+1
  bugs_df = bugs
  input_file <- data.table::fread(paste0(output_dir, "/",i), header=T, check.names=F)
  colnames(input_file)[1] = "Sample_ID"  
  
  input_file = input_file[gtools::mixedorder(input_file$Sample_ID), ]
  bugs_df = bugs_df[gtools::mixedorder(rownames(bugs_df)), ]
  
  bugs_subset = bugs_df[rownames(bugs_df) %in% input_file$Sample_ID,]
  input_file = input_file[input_file$Sample_ID %in% rownames(bugs_subset),]
  
  bugs_subset <- as.data.frame(bugs_subset)
  input_file <- as.data.frame(input_file)
  row.names(input_file) <- input_file$Sample_ID
  
  #check that data frames match-up
  setdiff(row.names(bugs_subset), row.names(input_file))
  setdiff(row.names(input_file), row.names(bugs_subset))
  identical(row.names(bugs_subset), row.names(input_file))

  adonis_res_rsq_individual = vector()
  adonis_res_pval_individual = vector()
  for (col in colnames(input_file)[2:22]){
    print(col)
    diet_var = input_file[,col]
    adonis.univ = adonis2(bugs_subset ~ antibiotic + consent_age + diagnosis2 + diet_var, data = input_file, method = "bray", strata = input_file$Participant_ID,by="margin")
    adonis_res_rsq_individual[col] = adonis.univ[4,]$R2
    adonis_res_pval_individual[col] = adonis.univ[4,]$`Pr(>F)`
    diagnosis = input_file$diagnosis[1]
  }
  adonis_df = as.data.frame(cbind(adonis_res_rsq_individual, adonis_res_pval_individual,diagnosis))
  datafs[[n]] = adonis_df
  }

comb = do.call(rbind,datafs)
write.csv(comb, paste0(output_dir,"/hmp2_individual_diet_adonis_results.csv"), row.names = TRUE)

