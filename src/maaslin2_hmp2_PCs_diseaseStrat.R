############ all phenotypes ##################
setwd("/Users/tobynbranck/Documents/Effect_size/analysis/")
library("Maaslin2")
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/4.1/Resources/library")
library("grid", lib.loc="/Library/Frameworks/R.framework/Versions/4.1/Resources/library")
library(permute)
library(vegan)
library(tidyverse)
library(lme4)
library(factoextra)
library(ade4)
library(FactoMineR)
library(psych)
library(stringr)
library(reshape2)

bugs <- data.table::fread("hmp_species_df.csv", header=T, check.names=F)

colnames(bugs)[1] = "Sample_ID"  
arrange(bugs, desc(Sample_ID))

week_var = bugs[,"bin"]
rownames(week_var) = bugs$Sample_ID

input_file <- data.table::fread("hmp2_meta_allPhenotypes.csv", header=T, check.names=F)
colnames(input_file)[1] = "Sample_ID"  
arrange(input_file, desc(Sample_ID))

bugs = bugs[bugs$Sample_ID %in% input_file$Sample_ID,]
input_file = input_file[input_file$Sample_ID %in% bugs$Sample_ID,]

bugs <- as.data.frame(bugs)
input_file <- as.data.frame(input_file)
row.names(bugs) <- bugs$Sample_ID
row.names(input_file) <- input_file$Sample_ID

input_file <- na.omit(input_file)

bugs = bugs[row.names(bugs) %in% row.names(input_file), ]
identical(row.names(bugs), row.names(input_file))
input_file = input_file[match(row.names(bugs), row.names(input_file)), ]
identical(row.names(bugs), row.names(input_file))

str(input_file)
input_file_diet = input_file[, 2:20]

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

# PCA
ncomp=2
pca.western.prudent <- psych::principal(r=input_file_diet,
                                        residuals = TRUE,
                                        missing = TRUE,
                                        impute = 'median',
                                        rotate = 'varimax',
                                        nfactors = ncomp,
                                        scores = TRUE)
plot(pca.western.prudent$loadings[,1:2],   # x and y data
     pch=21,              # point shape
     bg="black",          # point color
     cex=0,               # point size
     main="HMP2 PCA Loadings")
text(pca.western.prudent$loadings[,1:2],             # sets position of labels
     labels=rownames(pca.western.prudent$loadings[,1:2]),cex=.7 )
scree(input_file_diet) #displays factors _and_ components
factors = pca.western.prudent$scores

bugs$Sample_ID = NULL
bugs$Participant_ID = NULL
bugs$bin = NULL
input_file$Sample_ID = NULL
input_file$bin = NULL
input_file$diagnosis = gsub('nonIBD', 'a_nonIBD', input_file$diagnosis)
input_file$diagnosis2 = gsub('nonIBD.non_dysbiosis', 'a_nonIBD.non_dysbiosis', input_file$diagnosis2)

### remove special characters in string
colnames(bugs) = sub(".*s__", "", colnames(bugs))

setdiff(row.names(factors), row.names(input_file))
setdiff(row.names(input_file), row.names(factors))
input_file = input_file[rownames(input_file) %in% rownames(factors),]
metadata_input = cbind(factors, input_file[, 20:24])
metadata_input$WeekBin = week_var$bin[match(rownames(metadata_input),rownames(week_var))]
### Split metadata_input data frame by diagnosis ###
X<-split(metadata_input, metadata_input$diagnosis)

#maaslin2
for (i in X) {
  i = as.data.frame(i)
  name = i$diagnosis[1]
  print(name)
  i$diagnosis = NULL
  i$diagnosis2 = gsub("^.*\\.", "", i$diagnosis2)
  i_discrete = i
  i_discrete$RC1 = discretize(i_discrete$RC1,breaks=3)
  i_discrete$RC2 = discretize(i_discrete$RC2,breaks=3)
  #Maaslin2(input_data = bugs, input_metadata = i, output = sprintf("hmp2_maaslin2_PCs_diseaseStrat_%s", name), random_effects = c("Participant_ID"),reference = c("WeekBin,(0, 7)"),max_significance = 1.0)
  Maaslin2(input_data = bugs, input_metadata = i_discrete, output = sprintf("hmp2_maaslin2_PCs_DISCRETE_diseaseStrat_%s", name), random_effects = c("Participant_ID"),reference = c("WeekBin,(0, 7)",'RC1,[-2.14,-0.139)','RC2,[-1.77,-0.47)'),max_significance = 1.0)
  }

hmp2_cd = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_diseaseStrat_CD/all_results.tsv", header=T, check.names=T))
hmp2_uc = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_diseaseStrat_UC/all_results.tsv", header=T, check.names=T))
hmp2_non = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_diseaseStrat_a_nonIBD/all_results.tsv", header=T, check.names=T))

keep = c("RC1","RC2")
hmp2_cd = hmp2_cd[hmp2_cd$metadata %in% keep, ]
hmp2_uc = hmp2_uc[hmp2_uc$metadata %in% keep, ]
hmp2_non = hmp2_non[hmp2_non$metadata %in% keep, ]

combined = as.data.frame(cbind(hmp2_cd$coef,hmp2_uc$coef,hmp2_non$coef))
names(combined)[1] = "CD_coeffs"
names(combined)[2] = "UC_coeffs"
names(combined)[3] = "NonIBD_coeffs"

ggplot(data = melt(combined), aes(x=variable, y=abs(value))) + geom_boxplot(aes(fill=variable))

######### scatter plots ##########
r = cbind(bugs, metadata_input)
write.csv(r,"HMP2_PCA_df2.csv",row.names = TRUE)
#https://www-nature-com.ezp-prod1.hul.harvard.edu/articles/s41586-019-1065-y

ggplot(r, aes(x=r$RC2, y = r$`Anaerostipes hadrus`, color = diagnosis)) +
  geom_point(alpha=0.35, size=2) +
  geom_smooth(aes(col=diagnosis), method="lm", se=F, fullrange=TRUE) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_colour_manual("HMP2 Diagnosis", 
                      breaks = c("a_nonIBD", "CD", "UC"),
                      values = c("azure4", "burlywood3", "indianred3")) +
  xlab("Dietary pattern RC2") +
  ylab("Anaerostipes hadrus")

ggplot(r, aes(x=r$RC1, y = r$`Anaerostipes hadrus`, color = diagnosis)) +
  geom_point(alpha=0.35, size=2) +
  geom_smooth(aes(col=diagnosis), method="lm", se=F, fullrange=TRUE) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_colour_manual("HMP2 Diagnosis", 
                      breaks = c("a_nonIBD", "CD", "UC"),
                      values = c("azure4", "burlywood3", "indianred3")) +
  xlab("Dietary pattern RC1") +
  ylab("Anaerostipes hadrus")

splitdf_uc = as.data.frame(split(r,r$diagnosis)[3])
splitdf_cd = as.data.frame(split(r,r$diagnosis)[2])
splitdf_nonibd = as.data.frame(split(r,r$diagnosis)[1])

ggplot(splitdf_uc, aes(x=splitdf_uc$UC.RC1, y = splitdf_uc$UC.Anaerostipes.hadrus)) +
  geom_point(alpha=0.35, size=2) +
  geom_smooth(method="lm", se=F, fullrange=TRUE)

ggplot(splitdf_cd, aes(x=splitdf_cd$CD.RC1, y = splitdf_cd$CD.Anaerostipes.hadrus)) +
  geom_point(alpha=0.35, size=2) +
  geom_smooth(method="lm", se=F, fullrange=TRUE)

ggplot(splitdf_nonibd, aes(x=splitdf_nonibd$a_nonIBD.RC1, y = splitdf_nonibd$a_nonIBD.Anaerostipes.hadrus)) +
  geom_point(alpha=0.35, size=2) +
  geom_smooth(method="lm", se=F, fullrange=TRUE)

#### Extra code for raw variables ###################################################################
#maaslin2
names(input_file) <- gsub(x = names(input_file), pattern = "\\(", replacement = "")  
names(input_file) <- gsub(x = names(input_file), pattern = "\\)", replacement = "")  
names(input_file) <- gsub(x = names(input_file), pattern = "\\.", replacement = "")  
names(input_file) <- gsub(x = names(input_file), pattern = "\\,", replacement = "")  
names(input_file) <- gsub(x = names(input_file), pattern = " ", replacement = "_")  
Maaslin2(input_data = bugs, input_metadata = input_file, output = "hmp2_maaslin2_raw", random_effects = c("Participant_ID"), reference = c("diagnosis,a_nonIBD","diagnosis2,a_nonIBD.non_dysbiosis"))