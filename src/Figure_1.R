setwd("/Users/tobynbranck/Documents/Effect_size/analysis/")
library("ggplot2")
library("grid")
library("vegan")
library("dbplyr")
library("plyr")
library(tidyverse)
library("phyloseq")
library(gtools)
library("viridis")
library("colorspace")
library("cowplot")
library("ggpubr", lib.loc="/Library/Frameworks/R.framework/Versions/4.1/Resources/library")
library("png", lib.loc="/Library/Frameworks/R.framework/Versions/4.1/Resources/library")
library("wesanderson", lib.loc="/Library/Frameworks/R.framework/Versions/4.1/Resources/library")
library("ape")
library(psych)

hmp2 <- data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp_species_df.csv", header=T, check.names=F)
mlvs <- data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_species.csv", header=T, check.names=F)
mouse <- data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mouse_species.csv", header=T, check.names=F)
hadza <- data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hadza_species.csv", header=T, check.names=F)

mouse_complex_species = data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mouse_complex_species.csv", header=T, check.names=F)
mouse_complex_meta = data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mouse_complex_metadata.csv", header=T, check.names=F)
backhed_species = data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/backhed_species_for_joint_ord.csv", header=T, check.names=F)
backhed_metadata = data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/backhed_meta_joint_ord.csv", header=T, check.names=F)
murphy_species = data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/murphy_species_for_joint_ord.csv", header=T, check.names=F)
murphy_metadata = data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/murphy_meta_joint_ord.csv", header=T, check.names=F)

droplist = c('SRR7351797','SRR7411394')
murphy_species = murphy_species[ ! murphy_species$V1 %in% droplist, ]
murphy_metadata = murphy_metadata[ ! murphy_metadata$V1 %in% droplist, ]

backhed_species = backhed_species[gtools::mixedorder(backhed_species$V1), ]
backhed_metadata = backhed_metadata[gtools::mixedorder(backhed_metadata$V1), ]
murphy_species = murphy_species[gtools::mixedorder(murphy_species$V1), ]
murphy_metadata = murphy_metadata[gtools::mixedorder(murphy_metadata$V1), ]

rownames(backhed_species) = backhed_species$V1
backhed_species$V1 = NULL

names(backhed_species)<-gsub("\\_"," ",names(backhed_species))
rownames(backhed_metadata) = backhed_metadata$V1
backhed_metadata$V1 = NULL
rownames(murphy_species) = murphy_species$V1
murphy_species$V1 = NULL
rownames(murphy_metadata) = murphy_metadata$V1
murphy_metadata$V1 = NULL

rownames(hmp2) = hmp2$V1
hmp2$V1 = NULL
hmp2$bin = NULL
hmp2$Participant_ID = NULL

rownames(mlvs) = mlvs$V1
mlvs$V1 = NULL
mlvs$SampNum = NULL
mlvs$WeekNum = NULL
mlvs$Participant_ID = NULL

rownames(mouse) = mouse$V1
mouse$V1 = NULL
rownames(hadza) = hadza$V1
hadza$V1 = NULL

hadza$`Collinsella massiliensis` = hadza$`.Collinsella. massiliensis`
hadza$`.Collinsella. massiliensis` = NULL

backhed_group_var = c(backhed_species$Group)
murphy_age_var = c(murphy_species$Host_Age)

backhed_species$Group = NULL
murphy_species$Host_Age = NULL
combined_df = rbind.fill(hmp2,mlvs,backhed_species,murphy_species,mouse)
all_rownames = c(rownames(hmp2),rownames(mlvs),rownames(backhed_species),rownames(murphy_species),rownames(mouse))
rownames(combined_df) = all_rownames
combined_df[is.na(combined_df)] <- 0
dist_matrix = vegdist(combined_df, method="bray",diag=TRUE)

metadata_df <- data.frame(matrix(ncol = 1, nrow = 2610))
rownames(metadata_df) = all_rownames
metadata_df$matrix.ncol...1..nrow...2610. = NULL
metadata_df$Study = "NA"
metadata_df$Study[1:768] = 'HMP2'
metadata_df$Study[769:1697] = 'MLVS'
metadata_df$Study[1698:1960] = 'Infants_Backhed'
metadata_df$Study[1961:2209] = 'Infants_Murphy'
metadata_df$Study[2210:2610] = 'Mouse'

pcoa = pcoa(dist_matrix)
pcoa_df = as.data.frame(pcoa$vectors[,1:2])

cairo_pdf(file = "figure1/jointplot_fig1.pdf", width = 4, height = 4)
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color=metadata_df$Study)) + 
  geom_point(alpha = 6/10, size=3) +
  scale_color_manual(values = c('navy','#CC0099','#FFCCCC','skyblue','lightsalmon4')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(.9, .2)) +
  labs(x=paste0("PC1: ",round(pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
  theme(legend.key.size = unit(.25, 'cm')) +
  #theme(legend.position="top", legend.direction="horizontal") +
  theme(axis.text.x = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=8)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = c(0.84, 0.17), legend.box.background = element_rect(color="black", size=.75)) +
  theme(legend.title=element_blank(), legend.margin = ggplot2::margin(4,4,4,4))
dev.off()  

####### INDIVIDUAL STUDY ORDINATIONS #############
### for hmp2 and mlvs and mouse, use Factor scores
### categorical coloring for NHP, hadza

### putting diet information in metadata_df
metadata_df$Diet = NA
metadata_df$Diet_color = NA

#### NHP ####
NHP_high_fat = c("SRR2866538", "SRR2866539", "SRR2866540", "SRR2866541", "SRR2866542", "SRR2866543", "SRR2866544",
                        "SRR2866545", "SRR2866546", "SRR2866547", "SRR2866548", "SRR2866549", "SRR2866550", "SRR2866551",
                        "SRR2866552", "SRR2866553", "SRR2866554", "SRR2866555", "SRR2866556", "SRR2866557", "SRR2866558",
                        "SRR2866559", "SRR2866560", "SRR2866561", "SRR2866562", "SRR2866563", "SRR2866564", "SRR2866565",
                        "SRR2866566", "SRR2866567", "SRR2866568", "SRR2866569", "SRR2866570", "SRR2866571", "SRR2866572",
                        "SRR2866573", "SRR2866574", "SRR2866575", "SRR2866576", "SRR2866577", "SRR2866578", "SRR2866579",
                        "SRR2866580", "SRR2866581", "SRR2866582", "SRR2866583", "SRR2866594", "SRR2866624", "SRR2866699",
                        "SRR2866700")
NHP_herbivore = c("SRR2866584", "SRR2866585", "SRR2866586", "SRR2866587", "SRR2866588", "SRR2866589", "SRR2866590",
                  "SRR2866591", "SRR2866592", "SRR2866593", "SRR2866595", "SRR2866596", "SRR2866597", "SRR2866598",
                  "SRR2866599", "SRR2866600", "SRR2866601", "SRR2866602", "SRR2866603", "SRR2866604", "SRR2866645",
                  "SRR2866674", "SRR2866681", "SRR2866697", "SRR2866698")

nhp_metadata <- data.frame(matrix(ncol = 1, nrow = 75))
rownames(nhp_metadata) = c(NHP_high_fat,NHP_herbivore)
nhp_metadata$Diet = NA
nhp_metadata$Diet_color = NA
nhp_metadata$matrix.ncol...1..nrow...75. = NULL
nhp_metadata[NHP_high_fat,"Diet"] = 'High fat'
nhp_metadata[NHP_herbivore,"Diet"] = 'Herbivore'
nhp_metadata[NHP_high_fat,"Diet_color"] = '#999999'
nhp_metadata[NHP_herbivore,"Diet_color"] = '#E69F00'

nhp <- data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/dada2_test/nonHumanPrim_Phenotp_.001abun.25prev_table.csv", header=T, check.names=F)
rownames(nhp) = nhp$V1
nhp$V1 = NULL
nhp_metadata = subset(nhp_metadata, row.names(nhp_metadata) %in% row.names(nhp))

nhp[is.na(nhp)] <- 0
nhp_dist_matrix = vegdist(nhp, method="bray",diag=TRUE)
nhp_pcoa = pcoa(nhp_dist_matrix)
nhp_pcoa_df = as.data.frame(nhp_pcoa$vectors[,1:2])

cairo_pdf(file = "figure1/NHP_ord_fig1.pdf", width=4, height = 4)
ggplot(nhp_pcoa_df, aes(x = Axis.1, y = Axis.2, color=nhp_metadata$Diet)) + 
  geom_point(alpha = 6/10, size=3) +
  scale_color_manual(values = c('lightskyblue4', 'maroon2')) +
  labs(colour = "Diet") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(.9, .1)) +
  labs(x=paste0("PC1: ",round(nhp_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(nhp_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
  theme(legend.title = element_text( size=10), legend.text=element_text(size=10)) +
  theme(legend.key.size = unit(.35, 'cm')) +
  theme(legend.position="top", legend.direction="horizontal") +
  theme(axis.text.x = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=8)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position=c(.3,1),legend.justification=c(1,1),
        legend.direction="vertical",
        legend.box="horizontal",
        legend.box.just = c("top"), 
        legend.background = element_rect(fill=alpha('white', 0)))
dev.off()
# hmp2 importing factor analysis, 1st factor
bugs_hmp2 <- data.table::fread("hmp_species_df.csv", header=T, check.names=F)
colnames(bugs_hmp2)[1] = "Sample_ID"  
arrange(bugs_hmp2, desc(Sample_ID))

grep_list_hmp2 = read.csv("hmp2_diet_vars_list.csv", header=FALSE)$V1
input_file_hmp2 <- data.table::fread("hmp2_meta_allPhenotypes.csv", header=T, check.names=F)
colnames(input_file_hmp2)[1] = "Sample_ID"  
arrange(input_file_hmp2, desc(Sample_ID))
bugs_hmp2 = bugs_hmp2[bugs_hmp2$Sample_ID %in% input_file_hmp2$Sample_ID,]
input_file_hmp2 = input_file_hmp2[input_file_hmp2$Sample_ID %in% bugs_hmp2$Sample_ID,]
bugs_hmp2 <- as.data.frame(bugs_hmp2)
input_file_hmp2 <- as.data.frame(input_file_hmp2)
row.names(bugs_hmp2) <- bugs_hmp2$Sample_ID
row.names(input_file_hmp2) <- input_file_hmp2$Sample_ID
input_file_hmp2 <- na.omit(input_file_hmp2)
bugs_hmp2 = bugs_hmp2[row.names(bugs_hmp2) %in% row.names(input_file_hmp2), ]
identical(row.names(bugs_hmp2), row.names(input_file_hmp2))
input_file_hmp2 = input_file_hmp2[match(row.names(bugs_hmp2), row.names(input_file_hmp2)), ]
identical(row.names(bugs_hmp2), row.names(input_file_hmp2))

input_file_diet_hmp2 = input_file_hmp2[,grep_list_hmp2]
bugs_hmp2$Sample_ID = NULL
names(input_file_diet_hmp2) = substring(names(input_file_diet_hmp2),1,18)

#food.fa.raw.hmp2 <- factanal(input_file_diet_hmp2, factors = 3, scores='regression')
#hmp2_fact_scores = as.data.frame(food.fa.raw.hmp2$scores)
ncomp = 2
food.fa.raw.hmp2 <- psych::principal(r=input_file_diet_hmp2,
                                     residuals = TRUE,
                                     missing = TRUE,
                                     impute = 'median',
                                     rotate = 'varimax',
                                     nfactors = ncomp,
                                     scores = TRUE)
hmp2_fact_scores = as.data.frame(food.fa.raw.hmp2$scores)

metadata_df$Diet[1:768] = hmp2_fact_scores$RC1[match(rownames(metadata_df)[1:768],rownames(hmp2_fact_scores))]
hmp2_meta = metadata_df[1:768,]

hmp2[is.na(hmp2)] <- 0
hmp2_dist_matrix = vegdist(hmp2, method="bray",diag=TRUE)
hmp2_pcoa = pcoa(hmp2_dist_matrix)
hmp2_pcoa_df = as.data.frame(hmp2_pcoa$vectors[,1:2])

hmp2_pcoa_df$Diet = hmp2_meta$Diet
hmp2_pcoa_df = hmp2_pcoa_df[!is.na(hmp2_pcoa_df$Diet),]

cairo_pdf(file = "figure1/HMP2_ord_fig1.pdf", width=4, height =4)
ggplot(hmp2_pcoa_df, aes(x = Axis.1, y = Axis.2, color=as.numeric(Diet))) + 
  geom_point(alpha = 6/10,size=3) +
  scale_color_gradient2(midpoint = 0, mid = "gray84", low = "#E69F00", high = "navy", space = "Lab", limits=c(-4,4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(.95, .85),legend.key = element_rect(fill="transparent",colour = "transparent")) +
  labs(colour = "Dietary pattern (PC1)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=paste0("PC1: ",round(hmp2_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(hmp2_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
  theme(legend.title = element_text( size=10), legend.text=element_text(size=10)) +
  theme(legend.key.size = unit(.45, 'cm')) +
  #theme(legend.position="top", legend.direction="horizontal") +
  theme(axis.text.x = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=8)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  #theme(legend.position = c(0.17, 0.84), legend.key = element_rect(fill = alpha('white',1)))
  theme(legend.position=c(.3,1),legend.justification=c(1,1),
      legend.direction="vertical",
      legend.box="horizontal",
      legend.box.just = c("top"), 
      legend.background = element_rect(fill=alpha('white', 0)))
dev.off()

### mlvs
bugs_mlvs <- data.table::fread("mlvs_species.csv", header=T, check.names=F)
colnames(bugs_mlvs)[1] = "Sample_ID"  
arrange(bugs_mlvs, desc(Sample_ID))

grep_list_mlvs = read.csv("mlvs_diet_vars_list.csv", header=FALSE)$V1
input_file_mlvs <- data.table::fread("mlvs_diet.csv", header=T, check.names=F)
colnames(input_file_mlvs)[1] = "Sample_ID"  
arrange(input_file_mlvs, desc(Sample_ID))

bugs_mlvs = bugs_mlvs[bugs_mlvs$Sample_ID %in% input_file_mlvs$Sample_ID,]
input_file_mlvs = input_file_mlvs[input_file_mlvs$Sample_ID %in% bugs_mlvs$Sample_ID,]
bugs_mlvs <- as.data.frame(bugs_mlvs)
input_file_mlvs <- as.data.frame(input_file_mlvs)
row.names(bugs_mlvs) <- bugs_mlvs$Sample_ID
row.names(input_file_mlvs) <- input_file_mlvs$Sample_ID

input_file_diet_mlvs = input_file_mlvs[, grep_list_mlvs]
bugs_mlvs$Sample_ID = NULL

ncomp = 2
food.fa.raw.mlvs <- psych::principal(r=input_file_diet_mlvs,
                                     residuals = TRUE,
                                     missing = TRUE,
                                     impute = 'median',
                                     rotate = 'varimax',
                                     nfactors = ncomp,
                                     scores = TRUE)
mlvs_fact_scores = as.data.frame(food.fa.raw.mlvs$scores)
metadata_df$Diet[769:1697] = mlvs_fact_scores$RC1[match(rownames(metadata_df)[769:1697],rownames(mlvs_fact_scores))]

mlvs_meta = metadata_df[769:1697,]

mlvs[is.na(mlvs)] <- 0
mlvs_dist_matrix = vegdist(mlvs, method="bray",diag=TRUE)
mlvs_pcoa = pcoa(mlvs_dist_matrix)
mlvs_pcoa_df = as.data.frame(mlvs_pcoa$vectors[,1:2])

mlvs_pcoa_df$Diet = mlvs_meta$Diet
mlvs_pcoa_df = mlvs_pcoa_df[!is.na(mlvs_pcoa_df$Diet),]

cairo_pdf(file = "figure1/MLVS_ord_fig1.pdf", width=4, height = 4)
ggplot(mlvs_pcoa_df, aes(x = Axis.1, y = Axis.2, color=as.numeric(Diet))) + 
  geom_point(alpha = 6/10,size=3) +
  scale_color_gradient2(midpoint = 0, mid = "gray84", low = "#E69F00", high = "#56B4E9",space = "Lab", limits=c(-4,5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(.95, .89),legend.key = element_rect(fill="transparent",colour = "transparent")) +
  labs(colour = "Dietary pattern (PC1)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=paste0("PC1: ",round(mlvs_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(mlvs_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
  theme(legend.title = element_text(size=10), legend.text=element_text(size=10)) +
  theme(legend.key.size = unit(.45, 'cm')) +
  theme(legend.position="top", legend.direction="horizontal") +
  theme(axis.text.x = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=8)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position=c(.3,1),legend.justification=c(1,1),
        legend.direction="vertical",
        legend.box="horizontal",
        legend.box.just = c("top"), 
        legend.background = element_rect(fill=alpha('white', 0)))
dev.off()
# mouse importing factor analysis, 1st and 2nd factors
rownames(mouse_complex_meta) = mouse_complex_meta$V1
mouse_complex_meta$V1 = NULL
rownames(mouse_complex_species) = mouse_complex_species$V1
mouse_complex_species$V1 = NULL

#factor analysis for mice
#food.fa.raw.mouse <- factanal(mouse_complex_meta[,4:11], factors = 3, scores='regression',lower = 0.1)
ncomp = 2
food.fa.raw.mouse <- psych::principal(r=mouse_complex_meta[,4:11],
                                     residuals = TRUE,
                                     missing = TRUE,
                                     impute = 'median',
                                     rotate = 'varimax',
                                     nfactors = ncomp,
                                     scores = TRUE)

mouse_fact_scores = as.data.frame(food.fa.raw.mouse$scores)
rownames(mouse_fact_scores) = rownames(mouse_complex_meta)

mouse_complex_species[is.na(mouse_complex_species)] <- 0
mouse_complex_species_dist_matrix = vegdist(mouse_complex_species, method="bray",diag=TRUE)
mouse_pcoa = pcoa(mouse_complex_species_dist_matrix)
mouse_pcoa_df = as.data.frame(mouse_pcoa$vectors[,1:2])

cairo_pdf(file = "figure1/MOUSE_ord_fig1.pdf", width=4, height = 4)
ggplot(mouse_pcoa_df, aes(x = Axis.1, y = Axis.2, color=mouse_fact_scores$RC1)) + 
  geom_point(alpha = 7/10,size=3) +
  scale_color_gradient2(midpoint = 0, mid = "snow2", low = "gray19", high = "red4",space = "Lab", limits=c(-4,4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(.95, .89),legend.key = element_rect(fill="transparent",colour = "transparent")) +
  labs(colour = "Mouse Diet Factor 1") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=paste0("PC1: ",round(mouse_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(mouse_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
  theme(legend.title = element_text( size=10), legend.text=element_text(size=10)) +
  theme(legend.key.size = unit(.45, 'cm')) +
  theme(legend.position="top", legend.direction="horizontal") +
  theme(axis.text.x = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=8)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position=c(.3,1),legend.justification=c(1,1),
        legend.direction="vertical",
        legend.box="horizontal",
        legend.box.just = c("top"), 
        legend.background = element_rect(fill=alpha('white', 0)))
dev.off()
#pcoa for backhed
backhed_dist_matrix = vegdist(backhed_species, method="bray",diag=TRUE)
backhed_pcoa = pcoa(backhed_dist_matrix)
backhed_pcoa_df = as.data.frame(backhed_pcoa$vectors[,1:2])

backhed_pcoa_df$feeding = backhed_metadata$feeding
backhed_pcoa_df = backhed_pcoa_df[!is.na(backhed_pcoa_df$feeding),]
backhed_pcoa_df$age = backhed_metadata$age
backhed_metadata$group = sub(".*_", "", backhed_metadata$Individual)
backhed_pcoa_df$group = backhed_metadata$group

backhed_plot = ggplot(backhed_pcoa_df, aes(x = Axis.1, y = Axis.2, color=feeding)) + 
  geom_point(alpha = 6/10,size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(.95, .89),legend.key = element_rect(fill="transparent",colour = "transparent")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=paste0("PC1: ",round(backhed_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(backhed_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
  theme(legend.title = element_text(size=10), legend.text=element_text(size=10)) +
  theme(legend.key.size = unit(.45, 'cm')) +
  theme(legend.position="top", legend.direction="horizontal") +
  theme(axis.text.x = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=7)) +
  theme(legend.title=element_text(size=7)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  +
  theme(legend.direction = "vertical",
        legend.box = "horizontal",
        legend.position = "right",
        legend.justification = c(0, 1)) +
  scale_color_manual(values=c('#999999','#E69F00','#56B4E9','violetred4','slateblue4')) +
  theme(legend.position=c(1.05,.33),legend.justification=c(1,1),
        legend.direction="vertical",
        legend.box="horizontal",
        legend.box.just = c("top"))
  theme(legend.title= element_blank())

cairo_pdf(file = "figure1/BackhedDiet_ord_fig1.pdf", width=4, height = 4)
backhed_plot
dev.off() 

backhed_pcoa_df$group <- factor(backhed_pcoa_df$group, levels = c("B", "4M", "12M"))

backhed_plot_age = ggplot(backhed_pcoa_df, aes(x = Axis.1, y = Axis.2, color=group)) + 
  geom_point(alpha = 6/10,size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=paste0("PC1: ",round(backhed_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(backhed_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
  theme(axis.text.x = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.title=element_text(size=9)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.direction = "vertical",
        legend.box = "horizontal",
        legend.position = c(.95, .3)) +
  scale_color_manual(values=c('#999999','#E69F00','#56B4E9'))

cairo_pdf(file = "figure1/BackhedAge_ord_fig1.pdf", width=4, height = 4)
backhed_plot_age
dev.off()

pdf(file = "figure1/backhed_age_feeding_ordinations.pdf", width=14, height = 5)
ggarrange(backhed_plot_age,backhed_plot, ncol = 2, nrow = 1)
dev.off()
#pcoa for murphy

#backhed_species$Group = backhed_species
#backhed_species_B = backhed_species[backhed_species$Group == 'B',]

murphy_dist_matrix = vegdist(murphy_species, method="bray",diag=TRUE)
murphy_pcoa = pcoa(murphy_dist_matrix)
murphy_pcoa_df = as.data.frame(murphy_pcoa$vectors[,1:2])

murphy_pcoa_df$feeding = murphy_metadata$feeding
murphy_pcoa_df$age = murphy_metadata$age
murphy_pcoa_df = murphy_pcoa_df[!is.na(murphy_pcoa_df$feeding),]

murphy_plot = ggplot(murphy_pcoa_df, aes(x = Axis.1, y = Axis.2, color=feeding)) + 
  geom_point(alpha = 6/10,size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(.95, .89),legend.key = element_rect(fill="transparent",colour = "transparent")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=paste0("PC1: ",round(murphy_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(murphy_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
  theme(legend.title = element_text(size=10), legend.text=element_text(size=10)) +
  theme(legend.key.size = unit(.45, 'cm')) +
  theme(legend.position="top", legend.direction="horizontal") +
  theme(axis.text.x = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.title=element_text(size=9)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.direction = "vertical",
        legend.box = "horizontal",
        legend.position = "right",
        legend.justification = c(0, 1)) +
  scale_color_manual(values=c('#E69F00','#56B4E9','violetred4','slateblue4')) +
  theme(legend.position=c(1,1),legend.justification=c(1,1),
        legend.direction="vertical",
        legend.box="horizontal",
        legend.box.just = c("top"))

cairo_pdf(file = "figure1/MurphyDiet_ord_fig1.pdf", width=4, height = 4)
murphy_plot
dev.off()

murphy_plot_age = ggplot(murphy_pcoa_df, aes(x = Axis.1, y = Axis.2, color=age)) + 
  geom_point(alpha = 6/10,size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(.95, .89),legend.key = element_rect(fill="transparent",colour = "transparent")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=paste0("PC1: ",round(murphy_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(murphy_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
  theme(legend.title = element_text(size=10), legend.text=element_text(size=10)) +
  theme(legend.key.size = unit(.45, 'cm')) +
  theme(legend.position="top", legend.direction="horizontal") +
  theme(axis.text.x = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.title=element_text(size=9)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  +
  theme(legend.direction = "vertical",
        legend.box = "horizontal",
        legend.position = c(.80, .95),
        legend.justification = c(0, 1)) +
  scale_color_manual(values=c('#999999','#E69F00','#56B4E9'))

cairo_pdf(file = "figure1/MurphyAge_ord_fig1.pdf", width=4, height = 4)
murphy_plot_age
dev.off()

cairo_pdf(file = "figure1/murphy_age_feeding_ordinations.pdf", width=14, height = 5)
ggarrange(murphy_plot_age,murphy_plot, ncol = 2, nrow = 1)
dev.off()
#hadza seasonal
hadza_seas_sp = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hadza_seasonal_species.csv", header=T, check.names=F))
rownames(hadza_seas_sp) = hadza_seas_sp$V1
hadza_seas_sp$V1 = NULL
hadza_seas_metadata = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hadza_seasonal_metadata_processed.csv", header=T, check.names=F))
rownames(hadza_seas_metadata) = hadza_seas_metadata$V1
hadza_seas_metadata$V1 = NULL
hadza_seas_metadata = subset(hadza_seas_metadata, row.names(hadza_seas_metadata) %in% row.names(hadza_seas_sp))

hadza_seas_sp[is.na(hadza_seas_sp)] <- 0
hadza_seas_sp_dist_matrix = vegdist(hadza_seas_sp, method="bray",diag=TRUE)
hadza_seas_pcoa = pcoa(hadza_seas_sp_dist_matrix)
hadza_seas_pcoa_df = as.data.frame(hadza_seas_pcoa$vectors[,1:2])
hadza_seas_metadata$Season_binary = hadza_seas_metadata$Season
hadza_seas_metadata$Season_binary[hadza_seas_metadata$Season_binary == "2013-LD"] = 'Dry'
hadza_seas_metadata$Season_binary[hadza_seas_metadata$Season_binary == "2014-ED"] = 'Dry'
hadza_seas_metadata$Season_binary[hadza_seas_metadata$Season_binary == "2014-LD"] = 'Dry'
hadza_seas_metadata$Season_binary[hadza_seas_metadata$Season_binary == "2014-EW"] = 'Wet'
hadza_seas_metadata$Season_binary[hadza_seas_metadata$Season_binary == "2014-LW"] = 'Wet'


cairo_pdf(file = "figure1/Hadza_ord_fig1.pdf", width=4, height = 4)
ggplot(hadza_seas_pcoa_df, aes(x = Axis.1, y = Axis.2, color=hadza_seas_metadata$Season_binary)) + 
  geom_point(alpha = 6/10, size=3) +
  scale_color_manual(values = c("aquamarine4","darkgoldenrod1")) +
  #scale_fill_manual(values=c("aquamarine4","bisque4")) +
  labs(colour = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  #theme(legend.position = c(.9, .1)) +
  theme(legend.position="top", legend.direction="horizontal") +
  theme(plot.title = element_blank()) +
  theme(legend.title = element_text(size=6), legend.text=element_text(size=5.5)) +
  theme(legend.key.size = unit(.23, 'cm')) +
  labs(x=paste0("PC1: ",round(hadza_seas_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(hadza_seas_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
  theme(legend.title = element_text(size=8), legend.text=element_text(size=8)) +
  theme(legend.key.size = unit(.2, 'cm')) +
  theme(legend.position="top", legend.direction="horizontal") +
  theme(axis.text.x = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=10)) +
  theme(legend.title=element_text(size=10)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position=c(1,1.05),legend.justification=c(1,1),
        legend.direction="vertical",
        legend.box="horizontal",
        legend.box.just = c("top"), 
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.box.background = element_rect(colour = "black"))
dev.off()  
##### combining into one plot
# row_one = ggarrange(hmp2_plot + rremove("xlab"), mlvs_plot + rremove("xlab") + rremove("ylab"))
# row_two = ggarrange(hadza_seas_plot,nhp_plot + rremove("ylab"))
# 
# pdf(file = "fig1_indiv_ordinations_vertical.pdf", width=10, height = 10)
# ggarrange(row_one,row_two, ncol = 1, nrow = 2)
# dev.off()
# 
# pdf(file = "combined_ordinations.pdf", width=20, height = 4)
# ggarrange(joint_plot, hmp2_plot,mlvs_plot,hadza_seas_plot,nhp_plot, ncol = 5, nrow = 1)
# dev.off()
# 
# ######## overview fig + ordinations ##########
# img1 <- readPNG("/Users/tobynbranck/Desktop/overview6_21.png")
# im_A <- ggplot() + 
#   background_image(img1) +
#   theme(
#     panel.background = element_rect(fill = "transparent"), # bg of the panel
#     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#     panel.grid.major = element_blank(), # get rid of major grid
#     panel.grid.minor = element_blank(), # get rid of minor grid
#     legend.background = element_rect(fill = "transparent"),# get rid of legend bg
#     legend.box.background = element_rect(fill = "transparent") + # get rid of legend panel bg
#   # This ensures that the image leaves some space at the edges
#   theme(plot.margin = ggplot2::margin(t=1, l=1, r=1, b=1, unit = "cm"))) +
#   theme(panel.border = element_blank()) 

# metadata histograms for Figure 1
hmp2_meta = as.data.frame(data.table::fread('/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_meta_allPhenotypes.csv', header=T, check.names=F))
mlvs_meta = as.data.frame(data.table::fread('/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_diet.csv', header=T, check.names=F))
nhp_meta = as.data.frame(data.table::fread('/Users/tobynbranck/Documents/Effect_size/data/nonhumanPrimate_metadata.txt', header=T, check.names=F))

#hadza_seasonal_meta$`Library Name` = as.character(hadza_seasonal_meta$`Library Name`)
#hadza_seasonal_meta = hadza_seasonal_meta[hadza_seasonal_meta$`Library Name` %in% hadza_seas_samps_keep$V1,]

# hmp2 diagnosis plot
cairo_pdf(file = "figure1/diagnosis_fig1SUPP.pdf", width=4, height = 4)
ggplot(hmp2_meta[order(hmp2_meta$diagnosis2, decreasing = T),], aes(x = factor(diagnosis, levels=c("CD","UC","nonIBD")), fill=factor(diagnosis2, levels=c("CD.non_dysbiosis","CD.dysbiosis","UC.non_dysbiosis","UC.dysbiosis","nonIBD.non_dysbiosis", "nonIBD.dysbiosis" ))))+
  geom_bar(position = 'stack') +
  xlab("HMP2 Disease Phenotypes") + ylab("Number Samples") +
  scale_fill_manual(values = c("#CC0099","hotpink4","lemonchiffon3","lemonchiffon4","azure3","azure4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.key.size = unit(0.23, "cm")) +
  theme(legend.position = c(0.8, 0.9)) +
  theme(legend.background = element_rect(fill = NA)) +
  #theme(legend.title=element_blank()) +
  labs(fill = "Diagnosis") +
  theme(axis.text.x = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=8)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

# age plot (hmp2, mlvs, eventually diabimmune + backhed)
b_meta = backhed_metadata[,c("Individual","age")]
b_meta$"Age (years)" = b_meta$age / 365
b_meta$age = NULL
b_meta$Individual = as.character(b_meta$Individual)
colnames(b_meta) = c("Participant_ID", "consent_age")

murphy_meta = murphy_metadata[,c("Individual","age")]
murphy_meta$age[murphy_meta$age == 'baseline'] = 0
murphy_meta$age[murphy_meta$age == 'three months'] = 0.246
murphy_meta$age[murphy_meta$age == 'twelve months'] = 1
murphy_meta$Individual = as.character(murphy_meta$Individual)
colnames(murphy_meta) = c("Participant_ID", "consent_age")
murphy_meta$consent_age = as.numeric(murphy_meta$consent_age)

h_meta = hmp2_meta[,c('Participant_ID' , 'consent_age')]
m_meta = mlvs_meta[,c('Participant_ID' , 'consent_age')]
h_meta$consent_age = as.numeric(h_meta$consent_age)
m_meta$consent_age = as.numeric(m_meta$consent_age)
h_meta$Participant_ID = as.character(h_meta$Participant_ID)
m_meta$Participant_ID = as.character(m_meta$Participant_ID)
haz_meta = hadza_seas_metadata[,c('host_subject_id','Host_Age')]
colnames(haz_meta) <- c("Participant_ID", "consent_age")
haz_meta$Participant_ID = as.character(haz_meta$Participant_ID)

agedf = dplyr::bind_rows(h_meta, m_meta,haz_meta,b_meta,murphy_meta)
agedf$Study = NA
agedf$Study[1:768] = "HMP2"
agedf$Study[769:1691] = "MLVS"
agedf$Study[1692:1766] = "Hadza seasonal"
agedf$Study[1767:2029] = "Infants_Backhed"
agedf$Study[2030:2278] = "Infants_Murphy"

agedf = na.omit(agedf)
agedf_groupedbyParticipant = agedf %>%
  group_by(Participant_ID, Study) %>% 
  summarise_each(funs(mean))

haz_meta = haz_meta[haz_meta$consent_age >= 18,]
agedf_groupedbyParticipant$Study[agedf_groupedbyParticipant$Study == "Hadza seasonal"] = 'Hadza'

cairo_pdf(file = "figure1/age_histogram_fig1.pdf", width=4, height = 4)
ggplot(agedf_groupedbyParticipant, aes(x = consent_age, fill = factor(Study, levels = c("HMP2","MLVS","Hadza","Infants_Backhed","Infants_Murphy")))) +                       # Draw overlaying histogram
  geom_histogram(position = "identity", alpha = 0.6, bins = 20) + 
  xlab("Age (years)") + ylab("# Individuals") +
  labs(fill = " ") +
  scale_fill_manual(values = c('navy','skyblue','springgreen3','#CC0099','#FFCCCC')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.key.size = unit(0.35, "cm")) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=8)) +
  theme(legend.position = c(0.35, 0.8)) +
  theme(legend.background = element_rect(fill = NA)) +
  theme(axis.text.x = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=8)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

#### longitudinal data - number of samples per individual
rownames(nhp_meta) = nhp_meta$Run
nhp_meta = nhp_meta[rownames(nhp_meta) %in% rownames(nhp),]
nhp_meta$Participant_ID = nhp_meta$Run
nhp_meta$Run = NULL
n_meta = as.data.frame(nhp_meta[,c('Participant_ID')])
n_meta$Participant_ID = n_meta$`nhp_meta[, c("Participant_ID")]`
n_meta$`nhp_meta[, c("Participant_ID")]` = NULL
agedf = dplyr::bind_rows(agedf, n_meta)
agedf$Study[2279:2312] = "Nonhuman Primates"
agedf[nrow(agedf)+62,] <- NA
agedf$Study[2313:2374] = "Mouse"
agedf$Participant_ID[2313:2374] = mouse_complex_meta$Individual

num_samples = agedf %>% dplyr::count(Participant_ID, Study)
num_samples$Study[num_samples$Study=='Hadza seasonal'] = 'Hadza'

cairo_pdf(file = "figure1/longitudinality_fig1SUPP.pdf", width=4, height = 4)
ggplot(num_samples, aes(x=factor(Study, levels = c("HMP2", "MLVS", "Hadza", "Infants_Backhed","Infants_Murphy","Nonhuman Primates","Mouse")), y=n, color=Study)) + 
  geom_boxplot() +
  xlab("Study") + 
  ylab("Time points / subject") +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_color_manual(values=c('springgreen3','navy','#CC0099','#FFCCCC','skyblue','lightsalmon4','orange2')) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete(labels= c("HMP2", "MLVS", "Hadza", "Infants_Backhed","Infants_Murphy","NHP","Mouse")) +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(size = 12,angle = 45, vjust = 1, hjust = 1)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=8)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()
# hadza seasonal information
cairo_pdf(file = "figure1/hazda_seasons_fig1SUPP.pdf", width=4, height = 4)
ggplot(hadza_seas_metadata) + 
  geom_bar(aes(x = factor(Season, levels=c("2013-LD","2014-EW",'2014-LW','2014-ED','2014-LD')),fill=factor(Season, levels=c("2013-LD","2014-EW",'2014-LW','2014-ED','2014-LD')))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Subseason") + 
  ylab("Number Individuals") +
  scale_fill_viridis_d() +
  scale_x_discrete(labels = c("2013-LD" = "2013\nLD", "2014-EW" = "2014\nEW", "2014-LW" = "2014\nLW", "2014-ED" = "2014\nED","2014-LD" = "2014\nLD")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=8)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()
######## Figure 1 ##########  
# ords = ggarrange(joint_plot,hmp2_plot,mlvs_plot,hadza_seas_plot,backhed_plot,
#           nhp_plot,mouse_plot, ncol = 3, nrow = 3)
# hist = ggarrange(samples, age, diagnosis, seasons, labels = c("B","C","D","E"),ncol = 2, nrow = 2)
# ords2 = ggarrange(joint_plot,hmp2_plot,mlvs_plot,hadza_seas_plot,backhed_plot,
#                          nhp_plot,mouse_plot, ncol = 2, nrow = 4)
# hist2 = ggarrange(samples, age, diagnosis, seasons, ncol = 2, nrow = 2)
# test = ggarrange(im_A,hist2,ncol=1,nrow=2,labels = c("A", "C"),heights = c(5, 7))
# 
# pdf("/Users/tobynbranck/Desktop/ordinations.pdf",width=7,height = 12)
# ggarrange(joint_plot,hmp2_plot,mlvs_plot,hadza_seas_plot,backhed_plot,
#                  nhp_plot,mouse_plot, labels = c("A","B","C","D","E","F","G",""), ncol = 2, nrow = 4)
# dev.off()
# 
# pdf("/Users/tobynbranck/Desktop/histograms.pdf",width=12,height = 3)
# ggarrange(samples, age, diagnosis, seasons, labels = c("A","B","C","D"),ncol = 4, nrow = 1)
# dev.off()
# 
# # ggarrange(test, ords2,
# #           labels = c("", "B"),
# #           ncol = 2, nrow = 1)
# pdf("/Users/tobynbranck/Desktop/overview_test3.pdf", width=16, height = 12)
# ggarrange(test, ords2,
#           labels = c("", "B"),
#           ncol = 2, nrow = 1)
# dev.off()
# 
# ords3 = ggarrange(joint_plot,hmp2_plot,mlvs_plot,hadza_seas_plot,backhed_plot,
#                   nhp_plot,mouse_plot, labels = c("F","G","H","I","J","K","L"),ncol = 5, nrow = 2)
# ovw_ords = ggarrange(im_A,hist,ncol=2,nrow=1,labels = c("A", ""),widths = c(9, 7))
# pdf("/Users/tobynbranck/Desktop/overview_test5.pdf", width=15, height = 12)
# ggarrange(ovw_ords, ords3,
#           ncol = 1, nrow = 2, heights = c(4,5))
# dev.off()
# 
# hist_revised = ggarrange(samples, age, labels = c("B","C"),ncol = 1, nrow = 2,heights = c(6,4))
# ords_revised = ggarrange(joint_plot,hmp2_plot,mlvs_plot,hadza_seas_plot, labels = c("D","E","F","G"),ncol = 4, nrow = 1)
# ords_revised2 = ggarrange(backhed_plot,nhp_plot,mouse_plot, labels = c("H","I","J"),ncol = 3, nrow = 1, widths=c(6,4,4))
# ovw_hist = ggarrange(im_A,hist_revised,ncol=2,nrow=1,labels = c("A", ""),widths = c(10, 4))
# pdf("/Users/tobynbranck/Documents/Effect_size/analysis/Figure_1_Jun21.pdf", width=15, height = 13)
# ggarrange(ovw_hist, ords_revised, ords_revised2,
#           ncol = 1, nrow = 3, heights = c(6,4,4))
# dev.off()

#### overview figure panel A quantitative #diet features vs. #samples/individuals
#adding a small value (.1) to one of the infant datasets to prevent overlap on scatter
inf_murphy_ids = data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/Infants_murphy/ids.txt", header=T, check.names=F)
inf_murphy_ids = as.data.frame(inf_murphy_ids)
avg_murphy = mean(table(inf_murphy_ids$Patient.ID)) ##avg number of samples / individual for murphy dataset

study <- c('HMP2','MLVS','Hadza_seasonal','infants_Backhed','infants_Murphy','Non-human primates','Mouse')
num_diet_features <- c(21, 22, 5, 2.1, 2, 2, 8)
aggregate(n ~ Study, num_samples, mean) #average number of samples / individual / study
avg_num_samples_indiv = c(12.19,3,1,2,avg_murphy,1,7.75)
NumSamps_DietRes <- data.frame(study, num_diet_features, avg_num_samples_indiv)
NumSamps_DietRes$num_diet_features = log(NumSamps_DietRes$num_diet_features)
NumSamps_DietRes$avg_num_samples_indiv = log(NumSamps_DietRes$avg_num_samples_indiv)

pdf("overview_scatter.pdf", width=10, height = 10)
ggplot(NumSamps_DietRes, aes(x=avg_num_samples_indiv, y=num_diet_features,color=study)) + 
  geom_point(alpha = 0.5, size = 10,show.legend = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Level of time series data (Log average number samples/individual)") + 
  ylab("Resolution of dietary data (Log number of dietary variables)") +
  geom_text(label=NumSamps_DietRes$study,size=10,vjust=-.85,hjust=.005,show.legend = FALSE) +
  scale_fill_brewer() +
  theme(axis.text=element_text(size=22),axis.title=element_text(size=22))
dev.off()

######## EXTRA CODE (older): FOR BACKHED PRIOR TO BY-AGE FILTERING #########
#pcoa for backhed
# rownames(backhed_metadata) = backhed_metadata$Individual
# backhed_metadata$V1 = NULL
# backhed_metadata$samps = NULL
# rownames(backhed_metadata) = mixedsort(rownames(backhed_metadata))
# 
# backhed_species[is.na(backhed_species)] <- 0
# rownames(backhed_species) = backhed_species$Individual
# rownames(backhed_species) = mixedsort(rownames(backhed_species))
# backhed_species = backhed_species %>%
#   filter(rownames(backhed_species) %in% rownames(backhed_metadata))
# backhed_species$Individual = NULL
# 
# backhed_species_dist_matrix = vegdist(backhed_species, method="bray",diag=TRUE)
# backhed_pcoa = pcoa(backhed_species_dist_matrix)
# backhed_pcoa_df = as.data.frame(backhed_pcoa$vectors[,1:2])
# 
# backhed_plot = ggplot(backhed_pcoa_df, aes(x = Axis.1, y = Axis.2, color=backhed_metadata$feeding_strategy, shape = backhed_metadata$Time_point)) + 
#   geom_point(alpha = 8/10, size=3) +
#   #labs(colour = "Diet switch") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   theme(legend.position = c(.9, .1)) +
#   labs(x=paste0("PC1: ",round(backhed_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
#        y=paste0("PC2: ",round(backhed_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
#   theme(legend.title = element_text( size=10), legend.text=element_text(size=10)) +
#   theme(legend.key.size = unit(.35, 'cm')) +
#   #theme(legend.position="top", legend.direction="horizontal") +
#   theme(legend.position="right", legend.direction = "vertical") +
#   scale_color_manual(values = c("exclusively breastfeeding" = "blue4","mixed feeding" = "lightgoldenrod",
#                                 "exclusively formula feeding" ="dark#CC0099","no breastfeeding"="darksalmon","any breastfeeding"="deepskyblue")) +
#   labs(color = "Feeding type", shape = "Time") +
#   theme(axis.text.x = element_text(size = 12)) + 
#   theme(axis.text.y = element_text(size = 12)) +
#   theme(axis.title.x = element_text(size = 12)) + 
#   theme(axis.title.y = element_text(size = 12)) +
#   theme(legend.text=element_text(size=8)) +
#   theme(legend.title=element_text(size=8)) +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) 

