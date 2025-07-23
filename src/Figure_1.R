### Figure 1 and accompanying supplemental plots ###
library(ggplot2)
library(grid)
library(vegan)
library(dbplyr)
library(plyr)
library(tidyverse)
library(phyloseq)
library(gtools)
library(viridis)
library(colorspace)
library(cowplot)
library(ggpubr)
library(png)
library(wesanderson)
library(ape)
library(psych)

#set up directory path
input_dir = file.path("/home","tobynb","diet_test","input",fsep = "/")
output_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")

#import taxonomic profiles for cohorts/datasets with wmgx data (Figure 4A does not include Hazda or nonhuman primate data sets - 16S)
hmp2 <- data.table::fread(paste0(output_dir,"/hmp_species_df.csv"), header=T, check.names=F) #post filtering; final # samples
mlvs <- data.table::fread(paste0(output_dir,"/mlvs_species.csv"), header=T, check.names=F) #post filtering; final # samples
mouse <- data.table::fread(paste0(output_dir,"/mouse_species.csv"), header=T, check.names=F) #post filtering; final # samples
#hadza <- data.table::fread(paste0(output_dir,"/hadza_species.csv"), header=T, check.names=F)
backhed_species = data.table::fread(paste0(output_dir,"/backhed_species_for_joint_ord.csv"), header=T, check.names=F) #post filtering by age group; final # samples
murphy_species = data.table::fread(paste0(output_dir,"/murphy_species_for_joint_ord.csv"), header=T, check.names=F) #post filtering by age group; final # samples

#importing metadata for hmp2 and MLVS
input_file_hmp2 <- data.table::fread(paste0(output_dir,"/hmp2_meta_allPhenotypes.csv"), header=T, check.names=F)
names(input_file_hmp2)[1] = "Sample_ID"

input_file_mlvs <- data.table::fread(paste0(output_dir,"/mlvs_diet.csv"), header=T, check.names=F)
names(input_file_mlvs)[1] = "Sample_ID"

#importing metadata for infant cohorts
backhed_metadata = data.table::fread(paste0(output_dir,"/backhed_meta_joint_ord.csv"), header=T, check.names=F)
murphy_metadata = data.table::fread(paste0(output_dir,"/murphy_meta_joint_ord.csv"), header=T, check.names=F)

#ordering samples in infant taxonomic profiles and metadata
backhed_species = backhed_species[gtools::mixedorder(backhed_species$V1), ]
backhed_metadata = backhed_metadata[gtools::mixedorder(backhed_metadata$V1), ]
murphy_species = murphy_species[gtools::mixedorder(murphy_species$V1), ]
murphy_metadata = murphy_metadata[gtools::mixedorder(murphy_metadata$V1), ]

#set rownames in infant datasets
rownames(backhed_species) = backhed_species$V1
backhed_species$V1 = NULL
names(backhed_species)<-gsub("\\_"," ",names(backhed_species))
rownames(backhed_metadata) = backhed_metadata$V1
backhed_metadata$V1 = NULL

rownames(murphy_species) = murphy_species$V1
murphy_species$V1 = NULL
rownames(murphy_metadata) = murphy_metadata$V1
murphy_metadata$V1 = NULL

#isolate age group variables in infant datasets
backhed_group_var = c(backhed_species$Group)
murphy_age_var = c(murphy_species$Host_Age)
backhed_species$Group = NULL
murphy_species$Host_Age = NULL

#import Hadza seasonal data sets
hadza_seas_sp = as.data.frame(data.table::fread(paste0(output_dir,"/hadza_seasonal_species.csv"), header=T, check.names=F))
rownames(hadza_seas_sp) = hadza_seas_sp$V1
hadza_seas_sp$V1 = NULL
hadza_seas_metadata = as.data.frame(data.table::fread(paste0(output_dir,"/hadza_seasonal_metadata_processed.csv"), header=T, check.names=F))
rownames(hadza_seas_metadata) = hadza_seas_metadata$V1
hadza_seas_metadata$V1 = NULL
hadza_seas_sp = subset(hadza_seas_sp, row.names(hadza_seas_sp) %in% row.names(hadza_seas_metadata))
hadza_seas_sp = hadza_seas_sp[gtools::mixedorder(rownames(hadza_seas_sp)), ]
hadza_seas_metadata = hadza_seas_metadata[gtools::mixedorder(rownames(hadza_seas_metadata)), ]

#importing taxonomic profiles and metadata for complex diet mouse samples
mouse_complex_species = data.table::fread(paste0(output_dir,"/mouse_complex_species.csv"), header=T, check.names=F)
mouse_complex_meta = data.table::fread(paste0(output_dir,"/mouse_complex_metadata.csv"), header=T, check.names=F)

#importing metadata that includes all mouse samples (even additional samples not used in this study). subsetting to include all 227 samples in the complex and refined diets
mouse_all_metadata = as.data.frame(data.table::fread(paste0(output_dir,"/mouse_all_metadata.csv"), header=T, check.names=F))
mouse_all_metadata = subset(mouse_all_metadata, row.names(mouse_all_metadata) %in% row.names(mouse))

#import NHP taxonomic profiles and metadata
nhp_metadata <- data.table::fread(paste0(output_dir,"/nonHumanPrim_processedmetadata.csv"), header=T, check.names=F)
rownames(nhp_metadata) = nhp_metadata$V1
nhp_metadata$V1 = NULL

nhp <- data.table::fread(paste0(output_dir,"/nonHumanPrim_Phenotp_.001abun.25prev_table.csv"), header=T, check.names=F)
rownames(nhp) = nhp$V1
nhp$V1 = NULL

#format HMP2, MLVS, and mouse data to include only species columns for joint ordination
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

##############################################################################
### FIGURE 1A: ordination of taxonomic profiles for wgmx cohorts/data sets ### 
##############################################################################
#combining taxonomic profiles for the wmgx datasets: hmp2, mlvs, backhed, murphy and mouse
combined_df = rbind.fill(hmp2,mlvs,backhed_species,murphy_species,mouse)
all_rownames = c(rownames(hmp2),rownames(mlvs),rownames(backhed_species),rownames(murphy_species),rownames(mouse))
rownames(combined_df) = all_rownames
#where a species didn't exist in the data frame before, set as 0
combined_df[is.na(combined_df)] <- 0
dist_matrix = vegdist(combined_df, method="bray",diag=TRUE)

#add column to each data frame that links sample to study
hmp2$Study = 'HMP2'
mlvs$Study = 'MLVS'
backhed_species$Study = 'Infants_Backhed'
murphy_species$Study = 'Infants_Murphy'
mouse$Study = 'Mouse'
all_studies = c(hmp2$Study,mlvs$Study,backhed_species$Study,murphy_species$Study,mouse$Study)

#use the sample names and linked study names to build a simple metadata data frame for plotting the joint ordination
metadata_df = data.frame(all_rownames, all_studies)
names(metadata_df)[2] = 'Study'
rownames(metadata_df) = metadata_df$all_rownames
metadata_df$all_rownames = NULL

pcoa = pcoa(dist_matrix)
pcoa_df = as.data.frame(pcoa$vectors[,1:2])

cairo_pdf(file = paste0(output_dir,"/jointplot_fig1.pdf"), width = 4, height = 4)
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
  theme(axis.text.x = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(text=element_text(family="Arial")) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=8)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = c(0.90, 0.22), legend.box.background = element_rect(color="black", size=.75)) +
  theme(legend.title=element_blank(), legend.margin = ggplot2::margin(4,4,4,4))
dev.off()  

#########################################################################################################################################
####### FIGURE 1C-D and SUPPLEMENTAL FIGURE 1A-C: Individual study ordinations of taxonomic profiles overlayed by diet info #############
#########################################################################################################################################
### for hmp2 and mlvs and mouse, use Factor scores
### categorical coloring for NHP, hadza

### adding diet information to the metadata_df
metadata_df$Diet = NA
metadata_df$Diet_color = NA

#### NHP ordination ####
nhp_dist_matrix = vegdist(nhp, method="bray",diag=TRUE)
nhp_pcoa = pcoa(nhp_dist_matrix)
nhp_pcoa_df = as.data.frame(nhp_pcoa$vectors[,1:2])

cairo_pdf(file = paste0(output_dir,"/NHP_ord_fig1.pdf"), width=4, height = 4)
ggplot(nhp_pcoa_df, aes(x = Axis.1, y = Axis.2, color=nhp_metadata$Host_diet)) + 
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
        legend.background = element_rect(fill=scales::alpha('white', 0))) +
  theme(text=element_text(family="Arial")) +
dev.off()

#### HMP2 ordination ####
#remove Study column from the hmp2 taxonomic matrix; table formatting
bugs_hmp2 = select(hmp2, -Study)
hmp2_samp_names = rownames(bugs_hmp2)
bugs_hmp2 = as.data.frame(bugs_hmp2)
rownames(bugs_hmp2) = hmp2_samp_names

#hmp2 metadata table formatting; check to make sure samples line up in metadata and taxonomic profiles
input_file_hmp2 <- as.data.frame(input_file_hmp2)
rownames(input_file_hmp2) <- input_file_hmp2$Sample_ID
identical(rownames(bugs_hmp2), rownames(input_file_hmp2))

#samples didn't line up, so applying mixedorder function
bugs_hmp2 = bugs_hmp2[gtools::mixedorder(rownames(bugs_hmp2)), ]
input_file_hmp2 = input_file_hmp2[gtools::mixedorder(rownames(input_file_hmp2)), ]
identical(rownames(bugs_hmp2), rownames(input_file_hmp2))

#isolate the diet variables in the metadata
grep_list_hmp2 = read.csv(paste0(output_dir,"/hmp2_diet_vars_list.csv"), header=FALSE)$V2[2:22]
input_file_diet_hmp2 = input_file_hmp2[,colnames(input_file_hmp2) %in% grep_list_hmp2]

#truncates diet names for ease during plotting
names(input_file_diet_hmp2) = substring(names(input_file_diet_hmp2),1,18)

#derive factor scores from the HMP2 diet matrix
ncomp = 2
food.fa.raw.hmp2 <- psych::principal(r=input_file_diet_hmp2,
                                     residuals = TRUE,
                                     missing = TRUE,
                                     impute = 'median',
                                     rotate = 'varimax',
                                     nfactors = ncomp,
                                     scores = TRUE)
hmp2_fact_scores = as.data.frame(food.fa.raw.hmp2$scores)

#calculate Bray-Curtis distances and perform ordinatio
hmp2_dist_matrix = vegdist(bugs_hmp2, method="bray",diag=TRUE)
hmp2_pcoa = pcoa(hmp2_dist_matrix)
hmp2_pcoa_df = as.data.frame(hmp2_pcoa$vectors[,1:2])

#attach first diet factor score to PCs derived from taxonomy profiles for plotting
hmp2_pcoa_df$Diet_Factor_1 = hmp2_fact_scores$RC1[match(rownames(hmp2_fact_scores),rownames(hmp2_pcoa_df))]

cairo_pdf(file = paste0(output_dir,"/HMP2_ord_fig1.pdf"), width=4, height =4)
ggplot(hmp2_pcoa_df, aes(x = Axis.1, y = Axis.2, color=as.numeric(Diet_Factor_1))) + 
  geom_point(alpha = 6/10,size=3) +
  scale_color_gradient2(midpoint = 0, mid = "gray84", low = "#E69F00", high = "navy", space = "Lab", limits=c(-4,4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(.95, .85),legend.key = element_rect(fill="transparent",colour = "transparent")) +
  labs(colour = "HMP2 dietary pattern (PC1)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=paste0("PC1: ",round(hmp2_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(hmp2_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
  theme(legend.title = element_text( size=10), legend.text=element_text(size=10)) +
  theme(legend.key.size = unit(.45, 'cm')) +
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
  theme(legend.position=c(.97,.25),legend.justification=c(1,1),
      legend.direction="vertical",
      legend.box="horizontal",
      legend.box.just = c("top"), 
      legend.background = element_rect(fill=scales::alpha('white', 0))) +
  theme(text=element_text(family="Arial"))
dev.off()

#### MLVS ordination ####
#remove Study column from the MLVS taxonomic matrix; table formatting
bugs_mlvs = select(mlvs, -Study)
mlvs_samp_names = rownames(bugs_mlvs)
bugs_mlvs = as.data.frame(bugs_mlvs)
rownames(bugs_mlvs) = mlvs_samp_names

#MLVS metadata table formatting; check to make sure samples line up in metadata and taxonomic profiles
input_file_mlvs <- as.data.frame(input_file_mlvs)
rownames(input_file_mlvs) <- input_file_mlvs$Sample_ID
identical(rownames(bugs_mlvs), rownames(input_file_mlvs))

#samples didn't line up, so applying mixedorder function
bugs_mlvs = bugs_mlvs[gtools::mixedorder(rownames(bugs_mlvs)), ]
input_file_mlvs = input_file_mlvs[gtools::mixedorder(rownames(input_file_mlvs)), ]
identical(rownames(bugs_mlvs), rownames(input_file_mlvs))

#isolate the diet variables in the metadata
input_file_diet_mlvs = input_file_mlvs %>% select(-c(Sample_ID,consent_age,antibiotic,Participant_ID,SampNum,WeekNum))

#derive factor scores from the MLVS diet matrix
food.fa.raw.mlvs <- psych::principal(r=input_file_diet_mlvs,
                                     residuals = TRUE,
                                     missing = TRUE,
                                     impute = 'median',
                                     rotate = 'varimax',
                                     nfactors = ncomp,
                                     scores = TRUE)
mlvs_fact_scores = as.data.frame(food.fa.raw.mlvs$scores)

#calculate Bray-Curtis distances and perform ordination
mlvs_dist_matrix = vegdist(bugs_mlvs, method="bray",diag=TRUE)
mlvs_pcoa = pcoa(mlvs_dist_matrix)
mlvs_pcoa_df = as.data.frame(mlvs_pcoa$vectors[,1:2])

#attach first diet factor score to PCs derived from taxonomy profiles for plotting
mlvs_pcoa_df$Diet_Factor_1 = mlvs_fact_scores$RC1[match(rownames(mlvs_fact_scores),rownames(mlvs_pcoa_df))]

cairo_pdf(file = paste0(output_dir,"/MLVS_ord_fig1.pdf"), width=4, height = 4)
ggplot(mlvs_pcoa_df, aes(x = Axis.1, y = Axis.2, color=as.numeric(Diet_Factor_1))) + 
  geom_point(alpha = 6/10,size=3) +
  scale_color_gradient2(midpoint = 0, mid = "gray84", low = "#E69F00", high = "#56B4E9",space = "Lab", limits=c(-4,5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(.95, .89),legend.key = element_rect(fill="transparent",colour = "transparent")) +
  labs(colour = "MLVS dietary pattern (PC1)") +
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
        legend.background = element_rect(fill=scales::alpha('white', 0))) +
  theme(text=element_text(family="Arial"))
dev.off()

#### Mouse ordination ####
mouse_complex_meta = as.data.frame(mouse_complex_meta)
rownames(mouse_complex_meta) = mouse_complex_meta$V1
mouse_complex_meta$V1 = NULL
mouse_complex_species = as.data.frame(mouse_complex_species)
rownames(mouse_complex_species) = mouse_complex_species$V1
mouse_complex_species$V1 = NULL

input_file_diet_mouse = mouse_complex_meta %>% select(-c(individual,sampling_day,diet_period))

#factor analysis for mice
food.fa.raw.mouse <- psych::principal(r=input_file_diet_mouse,
                                     residuals = TRUE,
                                     missing = TRUE,
                                     impute = 'median',
                                     rotate = 'varimax',
                                     nfactors = ncomp,
                                     scores = TRUE)
mouse_fact_scores = as.data.frame(food.fa.raw.mouse$scores)

#calculate Bray-Curtis distances and perform ordination
mouse_complex_species_dist_matrix = vegdist(mouse_complex_species, method="bray",diag=TRUE)
mouse_pcoa = pcoa(mouse_complex_species_dist_matrix)
mouse_pcoa_df = as.data.frame(mouse_pcoa$vectors[,1:2])

#attach first diet factor score to PCs derived from taxonomy profiles for plotting
mouse_pcoa_df$Diet_Factor_1 = mouse_fact_scores$RC1[match(rownames(mouse_fact_scores),rownames(mouse_pcoa_df))]

cairo_pdf(file = paste0(output_dir,"/MOUSE_ord_fig1.pdf"), width=4, height = 4)
ggplot(mouse_pcoa_df, aes(x = Axis.1, y = Axis.2, color=Diet_Factor_1)) + 
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
        legend.background = element_rect(fill=scales::alpha('white', 0))) +
  theme(text=element_text(family="Arial"))
dev.off()

#### Infants (Backhed et al.) ordination ####
identical(rownames(backhed_metadata),rownames(backhed_species))
backhed_species$Study = NULL
backhed_species_samples = rownames(backhed_species)
backhed_dist_matrix = vegdist(backhed_species, method="bray",diag=TRUE)
backhed_pcoa = pcoa(backhed_dist_matrix)
backhed_metadata_samples = rownames(backhed_metadata)
backhed_pcoa_df = as.data.frame(backhed_pcoa$vectors[,1:2])

backhed_pcoa_df$feeding = backhed_metadata$feeding
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
  theme(legend.title= element_blank()) +
  theme(text=element_text(family="Arial"))
    

cairo_pdf(file = paste0(output_dir,"/BackhedDiet_ord_fig1.pdf"), width=4, height = 4)
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
  scale_color_manual(values=c('#999999','#E69F00','#56B4E9')) +
  theme(text=element_text(family="Arial"))

cairo_pdf(file = paste0(output_dir,"/BackhedAge_ord_fig1.pdf"), width=4, height = 4)
backhed_plot_age
dev.off()

pdf(file = paste0(output_dir,"/backhed_age_feeding_ordinations.pdf"), width=14, height = 5)
ggarrange(backhed_plot_age,backhed_plot, ncol = 2, nrow = 1)
dev.off()

#### Infants (Murphy et al.) ordination ####
identical(rownames(murphy_metadata),rownames(murphy_species))
murphy_species$Study = NULL
murphy_species_samples = rownames(murphy_species)
murphy_dist_matrix = vegdist(murphy_species, method="bray",diag=TRUE)
murphy_pcoa = pcoa(murphy_dist_matrix)
murphy_metadata_samples = rownames(murphy_metadata)
murphy_pcoa_df = as.data.frame(murphy_pcoa$vectors[,1:2])

murphy_pcoa_df$feeding = murphy_metadata$feeding
murphy_pcoa_df$age = murphy_metadata$age

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
        legend.box.just = c("top")) +
  theme(text=element_text(family="Arial"))

cairo_pdf(file = paste0(output_dir,"/MurphyDiet_ord_fig1.pdf"), width=4, height = 4)
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
  scale_color_manual(values=c('#999999','#E69F00','#56B4E9')) +
  theme(text=element_text(family="Arial"))

cairo_pdf(file = paste0(output_dir,"/MurphyAge_ord_fig1.pdf"), width=4, height = 4)
murphy_plot_age
dev.off()

cairo_pdf(file = paste0(output_dir,"/murphy_age_feeding_ordinations.pdf"), width=14, height = 5)
ggarrange(murphy_plot_age,murphy_plot, ncol = 2, nrow = 1)
dev.off()

#### Hadza ordination ####
identical(rownames(hadza_seas_metadata),rownames(hadza_seas_sp))

hadza_seas_sp_dist_matrix = vegdist(hadza_seas_sp, method="bray",diag=TRUE)
hadza_seas_pcoa = pcoa(hadza_seas_sp_dist_matrix)
hadza_seas_pcoa_df = as.data.frame(hadza_seas_pcoa$vectors[,1:2])

#recoding the season variable to be a binary variable
hadza_seas_metadata <- hadza_seas_metadata %>%
  mutate(Season_binary = recode(Season,
                                "2013-LD" = "Dry",
                                "2014-ED" = "Dry",
                                "2014-LD" = "Dry",
                                "2014-EW" = "Wet",
                                "2014-LW" = "Wet"))

cairo_pdf(file = paste0(output_dir,"/Hadza_ord_fig1.pdf"), width=4, height = 4)
ggplot(hadza_seas_pcoa_df, aes(x = Axis.1, y = Axis.2, color=hadza_seas_metadata$Season_binary)) + 
  geom_point(alpha = 6/10, size=3) +
  scale_color_manual(values = c("aquamarine4","darkgoldenrod1")) +
  labs(colour = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
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
        legend.background = element_rect(fill=scales::alpha('white', 0)),
        legend.box.background = element_rect(colour = "black")) +
  theme(text=element_text(family="Arial"))
dev.off() 

#########################################################################################################################################
####### FIGURE 1F and SUPPLEMENTAL FIGURE 1D-E: metadata features #######################################################################
#########################################################################################################################################
#Supplemental Figure 1E: hmp2 diagnosis plot
cairo_pdf(file = paste0(output_dir,"/diagnosis_fig1SUPP.pdf"), width=4, height = 4)
ggplot(input_file_hmp2[order(input_file_hmp2$diagnosis2, decreasing = T),], aes(x = factor(diagnosis, levels=c("CD","UC","nonIBD")), fill=factor(diagnosis2, levels=c("CD.non_dysbiosis","CD.dysbiosis","UC.non_dysbiosis","UC.dysbiosis","nonIBD.non_dysbiosis", "nonIBD.dysbiosis" ))))+
  geom_bar(position = 'stack') +
  xlab("HMP2 Disease Phenotypes") + ylab("Number Samples") +
  scale_fill_manual(values = c("#CC0099","hotpink4","lemonchiffon3","lemonchiffon4","azure3","azure4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.key.size = unit(0.23, "cm")) +
  theme(legend.position = c(0.8, 0.9)) +
  theme(legend.background = element_rect(fill = NA)) +
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
        panel.background = element_blank()) +
  theme(text=element_text(family="Arial"))
dev.off()

#Figure 1G: age plot (hmp2, mlvs, hadza, backhed, murphy)
#isolating the age and individual participant pieces of the metadata
b_meta = backhed_metadata[,c("Individual","age")]
b_meta$"Age (years)" = b_meta$age / 365
b_meta$age = NULL
b_meta$Individual = as.character(b_meta$Individual)
b_meta$Individual = as.character(sub("_.*", "", b_meta$Individual))
colnames(b_meta) = c("Participant_ID", "consent_age")
b_meta$Study = "Infants_Backhed"

murphy_meta = murphy_metadata[,c("Individual","age")]
murphy_meta$age[murphy_meta$age == 'baseline'] = 0
murphy_meta$age[murphy_meta$age == 'three months'] = 0.246
murphy_meta$age[murphy_meta$age == 'twelve months'] = 1
murphy_meta$Individual = as.character(murphy_meta$Individual)
colnames(murphy_meta) = c("Participant_ID", "consent_age")
murphy_meta$consent_age = as.numeric(murphy_meta$consent_age)
murphy_meta$Study = "Infants_Murphy"

h_meta = input_file_hmp2[,c('Participant_ID' , 'consent_age')]
m_meta = input_file_mlvs[,c('Participant_ID' , 'consent_age')]
h_meta$consent_age = as.numeric(h_meta$consent_age)
m_meta$consent_age = as.numeric(m_meta$consent_age)
h_meta$Participant_ID = as.character(h_meta$Participant_ID)
m_meta$Participant_ID = as.character(m_meta$Participant_ID)
h_meta$Study = "HMP2"
m_meta$Study = "MLVS"

haz_meta = hadza_seas_metadata[,c('host_subject_id','Host_Age')]
colnames(haz_meta) <- c("Participant_ID", "consent_age")
haz_meta$Participant_ID = as.character(haz_meta$Participant_ID)
haz_meta$Study = "Hadza seasonal"

agedf = dplyr::bind_rows(h_meta, m_meta,haz_meta,b_meta,murphy_meta)

agedf_groupedbyParticipant = agedf %>%
  group_by(Participant_ID, Study) %>% 
  summarise_each(funs(mean))

agedf_groupedbyParticipant$Study[agedf_groupedbyParticipant$Study == "Hadza seasonal"] = 'Hadza'

cairo_pdf(file = paste0(output_dir,"/age_histogram_fig1.pdf"), width=4, height = 4)
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
        panel.background = element_blank()) +
  theme(text=element_text(family="Arial"))
dev.off()

#Supplemental Figure 1D: longitudinal data - number of samples per individual
identical(rownames(nhp_metadata), rownames(nhp))
n_meta = nhp_metadata
n_meta$Participant_ID = rownames(n_meta)
n_meta$consent_age = NA
n_meta$Study = "Nonhuman Primates"
n_meta$Host_diet = NULL

agedf = dplyr::bind_rows(agedf, n_meta)

mouse_meta = as.data.frame(mouse_all_metadata[,c('individual')])
names(mouse_meta)[1] = "Participant_ID"
mouse_meta$consent_age = NA
mouse_meta$Study = "Mouse"
agedf = dplyr::bind_rows(agedf, mouse_meta)

num_samples = agedf %>% dplyr::count(Participant_ID, Study)
num_samples$Study[num_samples$Study=='Hadza seasonal'] = 'Hadza'

cairo_pdf(file = paste0(output_dir,"/longitudinality_fig1SUPP.pdf"), width=4, height = 4)
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
        panel.background = element_blank()) +
  theme(text=element_text(family="Arial"))
dev.off()

#hadza seasonal information (not included in manuscript)
cairo_pdf(file = paste0(output_dir,"/hazda_seasons_fig1SUPP.pdf"), width=4, height = 4)
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
        panel.background = element_blank()) +
  theme(text=element_text(family="Arial"))
dev.off()

#### overview figure panel A quantitative #diet features vs. #samples/individuals
#adding a small value (.1) to one of the infant datasets to prevent overlap on scatter
inf_murphy_ids = data.table::fread(paste0(input_dir,"/Murphy_infants/Initial/ids.txt"), header=T, check.names=F)
inf_murphy_ids = as.data.frame(inf_murphy_ids)
avg_murphy = mean(table(inf_murphy_ids$Patient.ID)) ##avg number of samples / individual for murphy dataset

study <- c('HMP2','MLVS','Hadza_seasonal','infants_Backhed','infants_Murphy','Non-human primates','Mouse')
num_diet_features <- c(21, 22, 5, 2.1, 2, 2, 8)
aggregate(n ~ Study, num_samples, mean) #average number of samples / individual / study
avg_num_samples_indiv = c(12.19,3,1,2,avg_murphy,1,7.75)
NumSamps_DietRes <- data.frame(study, num_diet_features, avg_num_samples_indiv)
NumSamps_DietRes$num_diet_features = log(NumSamps_DietRes$num_diet_features)
NumSamps_DietRes$avg_num_samples_indiv = log(NumSamps_DietRes$avg_num_samples_indiv)

pdf(paste0(output_dir,"/overview_scatter.pdf"), width=10, height = 10)
ggplot(NumSamps_DietRes, aes(x=avg_num_samples_indiv, y=num_diet_features,color=study)) + 
  geom_point(alpha = 0.5, size = 10,show.legend = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Level of time series data (Log average number samples/individual)") + 
  ylab("Resolution of dietary data (Log number of dietary variables)") +
  geom_text(label=NumSamps_DietRes$study,size=10,vjust=-.85,hjust=.005,show.legend = FALSE) +
  scale_fill_brewer() +
  theme(axis.text=element_text(size=22),axis.title=element_text(size=22)) +
  theme(text=element_text(family="Arial"))
dev.off()
