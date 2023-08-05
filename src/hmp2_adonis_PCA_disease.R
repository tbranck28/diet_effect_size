setwd("/Users/tobynbranck/Documents/Effect_size/analysis/")
library("ggplot2")
library("grid")
library("vegan")
library(tidyverse)
library(lme4)
library("cowplot")
library(psy)
library(psych)
library(stats)
library(MicEco)

inputfile1 <- data.table::fread("df_CD.csv", header=T, check.names=F)
inputfile2 <- data.table::fread("df_nonIBD.csv", header=T, check.names=F)
inputfile3 <- data.table::fread("df_UC.csv", header=T, check.names=F)
allPts = rbind(inputfile1,inputfile2,inputfile3)

adonis_rsq = vector()
adonis_pval = vector()
results = vector()
datafs = list()
cohort = "hmp2"
bugs <- data.table::fread("hmp_species_df.csv", header=T, check.names=F)

input_file = allPts
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
input_file_diet = input_file[, 2:20]
adonis_label = vector()
######is diet different across disease phenotypes?#######
diet_adonis = adonis(input_file_diet ~ diagnosis, data = input_file, method = "manhattan", strata = input_file$Participant_ID)
#########################################################

bugs$Sample_ID = NULL
bugs$Participant_ID = NULL
bugs$bin = NULL
print (head(bugs))

ncomp <- 2
pca.western.prudent <- psych::principal(r=input_file_diet,
                                        residuals = TRUE,
                                        missing = TRUE,
                                        impute = 'median',
                                        rotate = 'varimax',
                                        nfactors = ncomp,
                                        scores = TRUE)

scree(input_file_diet) #displays factors _and_ components
factors = as.data.frame(pca.western.prudent$scores)
X<-split(input_file, input_file$diagnosis)
n=0
datalist = list()
for (i in X) {
  n=n+1
  input_file_df = i
  factors_df = factors %>% filter(row.names(factors) %in% rownames(input_file_df))
  bugs_df = bugs %>% filter(row.names(bugs) %in% rownames(input_file_df))
  print(dim(input_file_df))
  print(dim(factors_df))
  print(dim(bugs_df))
  print(input_file_df$diagnosis)
  print(summary(as.factor(input_file_df$antibiotic)))
  adonis.PCs = adonis(bugs_df ~ antibiotic + consent_age + diagnosis2 + factors_df$RC1 + factors_df$RC2, data = input_file_df, method = "bray", strata = input_file_df$Participant_ID)
  print (adonis.PCs)
  adonis.PCs_PartialOmega = adonis_OmegaSq(adonis.PCs, partial = TRUE)
  print (adonis.PCs_PartialOmega)
  adonis_rsq = adonis.PCs$aov.tab[4:5,]$R2
  adonis_pval = adonis.PCs$aov.tab[4:5,]$`Pr(>F)`
  adonis_label = rownames(adonis.PCs$aov.tab[4:5,])
  adonis_label = gsub('.*\\$', '', rownames(adonis.PCs$aov.tab[4:5,]))
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
ggplot(data=df_comb,aes(x=cat,y=df_comb$adonis_rsq, fill=state)) + geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values=c("lemonchiffon3", "azure3", "indianred3")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  ylab(~R^2) + xlab("") +
  labs(fill = "Diagnosis") 
write.csv(df_comb,"/Users/tobynbranck/Documents/Effect_size/analysis/HMP2_adonisPCs_byDisease_results.csv", row.names = FALSE)
loads = as.data.frame(unclass(pca.western.prudent$loadings))
write.csv(loads, "/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_combined_loadings.csv", row.names = TRUE)
