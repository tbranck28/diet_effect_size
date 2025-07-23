#This script formats and analyzes the microbiome/diet group information for the Hadza population
#resulting number of samples after pairing taxonomic profiles with dietary metadata: 75 (60 dry season; 15 wet season)
#resulting number of individuals: 75

#load packages
library(Maaslin2)
library(ggplot2)
library(grid)
library(vegan)
library(dplyr)
library(tidyverse)
library(gtools)
library(MicEco)

#set directory paths and load taxonomy profiles
input_dir = file.path("/home","tobynb","diet_test","input","hadza","Initial",fsep = "/")
output_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")
tax = as.data.frame(data.table::fread(paste0(input_dir,"/all_samples_taxonomy_closed_reference.tsv"), sep = "\t", header=T, check.names=F))

#formatting taxonomy table and aggregating counts by taxonomy
tax$taxonomy = gsub("; s__.*", "", tax$taxonomy)
rownames(tax) = tax$V1
asvs = tax$V1
tax$V1 = NULL
#tax$taxonomy = gsub("\\s*\\([^\\)]+\\)","",as.character(tax$taxonomy))
#tax = tax[- grep("k__Archaea", tax$taxonomy),]
#tax = tax %>% group_by(taxonomy) %>% summarise_all(list(sum))
tax = tax %>% group_by(taxonomy) %>% summarise_all(list(sum))
phy_dat = tax

#formatting taxonomic names
tax_levels = strsplit(as.character(tax$taxonomy), split = "; ", fixed = TRUE)
tax_levels = data.frame(t(sapply(tax_levels, function(x) x)))
row.names(tax_levels) = tax$taxonomy
names(tax_levels) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
tax_levels = data.frame(tax_levels)
tax_levels = subset(tax_levels, Domain != "k__ " &  Domain != "k__Eukaryota") #drops an ASV that mapped to Eukaryota
tax_levels$Domain = factor(tax_levels$Domain)
tax_levels$Domain = droplevels(tax_levels$Domain)
tax_levels$Phylum = ifelse(tax_levels$Phylum == "p__ ", "p__Unclassified", as.character(tax_levels$Phylum)) #replaces "p__" with "p__Unclassified"
phy_levels = tax_levels
tax_levels$cleaned_names = ifelse(tax_levels$Genus == "g__ ", as.character(tax_levels$Family), as.character(tax_levels$Genus))
tax_levels$cleaned_names = ifelse(tax_levels$cleaned_names == "f__ ", as.character(tax_levels$Order), as.character(tax_levels$cleaned_names))
tax_levels$cleaned_names = ifelse(tax_levels$cleaned_names == "o__ ", as.character(tax_levels$Class), as.character(tax_levels$cleaned_names))
tax_levels$cleaned_names = ifelse(tax_levels$cleaned_names == "c__ ", as.character(tax_levels$Phylum), as.character(tax_levels$cleaned_names))

#removes the Eukaryota ASV from the taxonomy table
tax = tax[tax$taxonomy %in% row.names(tax_levels), ]
tax = data.frame(tax)

#sets rownames of the taxonomy table to be the cleaned taxonomic names from the taxa levels table
row.names(tax) = tax_levels[match(tax$taxonomy, row.names(tax_levels)), "cleaned_names"]
tax$taxonomy = NULL

#calculates relative abundance
tax = data.frame(t(tax))
tax = tax/rowSums(tax)
rowSums(tax)
head(tax)

#note hadza_seasonal_metadata table came from SRA website metadata; 11358_prep_3753_20210706-110244.txt
#came from qiita (https://qiita.ucsd.edu/study/description/11358#, file downloaded after clicking on "Prep Info"). 
#Using this file to map the sample IDs (in the form of "11358.109" for example) to the SRR accession numbers ("SRRXXXX", found under
#the srr_accession column in the qiita file and is what is used in the all_samples_taxonomy_closed_reference.tsv and hadza_seasonal_metadata). 
#This mapping is needed because the SRR metadata was not clear as to which samples belonged to which metadata
otu_df = tax
metadata_df = as.data.frame(data.table::fread(paste0(input_dir,"/hadza_seasonal_metadata.txt"), header=T, check.names=F))
qiita_for_ID = as.data.frame(data.table::fread(paste0(input_dir,"/11358_prep_3753_20210706-110244.txt"), header=T, check.names=F, colClasses=c('character')))
qiita_for_ID = qiita_for_ID[,c("sample_name","srr_accession")]
sample_name = c(qiita_for_ID$sample_name)
metadata_df$`Library Name` = as.character(metadata_df$`Library Name`)
#subsets metadata file to include only the records (samples) we care about
metadata_df = metadata_df[metadata_df$`Library Name` %in% sample_name,]
metadata_df$sample_name = metadata_df$`Library Name`
metadata_df$`Library Name` = NULL
#adds srr_accession to metadata file by merging on sample name (sample name is in the form of "11358.109", for example)
metadata_df = merge(metadata_df, qiita_for_ID[, c("sample_name", "srr_accession")], by="sample_name")

#including adults only and "seasonality_unique_subjects == 1" which corresponds to individuals who are not repeated across populations (see original publication)
metadata_df = metadata_df[metadata_df$host_life_stage == 'adult',]
metadata_df = metadata_df[metadata_df$seasonality_unique_subjects == 1,]

#setting rownames to be sample names and dropping unneeded column
rownames(metadata_df) = metadata_df$srr_accession
metadata_df$srr_accession = NULL

#number of samples per season
factor_snps <- factor(metadata_df$Season)
summary(factor_snps)

#noting which codes refer to wet and dry seasons
wet_seasons = c("2014-EW","2014-LW")
dry_seasons = c("2013-LD","2014-ED","2014-LD")

#subset otu_df samples = to metadata samples
otu_df = otu_df[match(row.names(metadata_df), row.names(otu_df)), ]

#filtering by diet/season groups
#separate metadata tables by wet vs dry season
wet_season_df = metadata_df[metadata_df$Season %in% wet_seasons,]
dry_season_df = metadata_df[metadata_df$Season %in% dry_seasons,]

#separate otu table to match separated metadata-by-season tables
otu_df_wet = subset(otu_df, row.names(otu_df) %in% row.names(wet_season_df))
otu_df_dry = subset(otu_df, row.names(otu_df) %in% row.names(dry_season_df))
otu_df_wet = as.data.frame(t(otu_df_wet))
otu_df_dry = as.data.frame(t(otu_df_dry))

#filtering ASVs from wet season samples (less strict)
filtered_wet = otu_df_wet[rowSums(otu_df_wet > 0.001) >= dim(otu_df_wet)[2]*.25, ]
filtered_wet = as.data.frame(t(filtered_wet))
#filtering ASVs from otu dry season samples (less strict)
filtered_dry = otu_df_dry[rowSums(otu_df_dry > 0.001) >= dim(otu_df_dry)[2]*.25, ]
filtered_dry = as.data.frame(t(filtered_dry))

#filtering otu-wet samples (more strict)
filtered_strict_wet = otu_df_wet[rowSums(otu_df_wet > 0.001) >= dim(otu_df_wet)[2]*.25, ]
filtered_strict_wet = data.frame(filtered_strict_wet)
filtered_strict_wet = as.data.frame(t(filtered_strict_wet))
#filtering otu-dry samples (more strict)
filtered_strict_dry = otu_df_dry[rowSums(otu_df_dry > 0.001) >= dim(otu_df_dry)[2]*.25, ]
filtered_strict_dry = data.frame(filtered_strict_dry)
filtered_strict_dry = as.data.frame(t(filtered_strict_dry))

#find union of ASVs after apply prevalence/abundance filtering to the subsetted-by-season samples -
#two groups - high fat and herbivore; merge filtered_wet and filtered_dry into filtered
wet_asvs = colnames(filtered_wet)
dry_asvs = colnames(filtered_dry)
union_filt = union(wet_asvs,dry_asvs)
filtered = otu_df[,union_filt]
filtered = filtered[,mixedsort(colnames(filtered))]

write.csv(filtered, paste0(output_dir,"/hadza_seasonal_species.csv"), row.names = TRUE)

# if going to use 2 seasons for adonis, make new column for either wet or dry season
metadata_df$Binary_season = NA
metadata_df$Binary_season = ifelse(metadata_df$Season %in% wet_seasons, "Wet Season", "Dry Season")

write.csv(metadata_df, paste0(output_dir,"/hadza_seasonal_metadata_processed.csv"), row.names = TRUE)

#formatting season column
metadata_df$Season = gsub("-","",metadata_df$Season)

######## adonis section ##########
#checking to make sure the samples are in the same order for both the taxonomic profiles and metadata
identical(row.names(filtered), row.names(metadata_df))
#filtered = filtered[rownames(filtered) %in% rownames(metadata_df),]
metadata_df_adonis = metadata_df[rownames(metadata_df) %in% rownames(filtered),]
summary(factor(metadata_df_adonis$host_subject_id))

#note that in the manuscript, we report the R^2 for the binary season variable rather than the sub-seasons (since one of the sub-seasons has only 4 samples)
adonis_obj_binary_season = adonis2(filtered ~ metadata_df_adonis$Host_Age + metadata_df_adonis$Binary_season, data = metadata_df_adonis, method = "bray",na.rm=TRUE,by="margin")
adonis_obj_subseasons = adonis2(filtered ~ metadata_df_adonis$Host_Age + metadata_df_adonis$Season, data = metadata_df_adonis, method = "bray",na.rm=TRUE,by="margin")

adjust_bin = adonis_OmegaSq(adonis_obj_binary_season,partial = TRUE)
adjust_subseason = adonis_OmegaSq(adonis_obj_subseasons,partial = TRUE)

######## Maaslin2 section ########
#subsetting metadata for ease when inputting into maaslin2
name_list = c('Binary_season','Host_Age','host_subject_id')
metadata_df_maaslin = metadata_df[,name_list]
filtered = filtered[rownames(filtered) %in% rownames(metadata_df_maaslin),]
metadata_df_maaslin = metadata_df_maaslin[rownames(metadata_df_maaslin) %in% rownames(filtered),]
metadata_df_maaslin = as.data.frame(metadata_df_maaslin)
metadata_df_maaslin$Binary_season = gsub(" ","",metadata_df_maaslin$Binary_season)
Maaslin2(input_data = filtered, input_metadata = metadata_df_maaslin[,1:2], output = paste0(output_dir,"/hadza_seasonal_maaslin2_output_binary"),reference = c("Binary_season,DrySeason"),max_significance = 1)
