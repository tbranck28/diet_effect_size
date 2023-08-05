setwd("/Users/tobynbranck/Documents/Effect_size/analysis/")
library("Maaslin2", lib.loc="/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
library("grid", lib.loc="/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
library("vegan", lib.loc="/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
library("dbplyr", lib.loc="/Library/Frameworks/R.framework/Versions/4.1/Resources/library")
library(tidyverse)
library(gtools)
library(MicEco)
tax = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/Hadza_seasonal/all_samples_taxonomy_closed_reference.tsv", sep = "\t", header=T, check.names=F))
#aggragate counts by taxonomy
tax$taxonomy = gsub("; s__.*", "", tax$taxonomy)
rownames(tax) = tax$V1
asvs = tax$V1
tax$V1 = NULL
#tax$taxonomy = gsub("\\s*\\([^\\)]+\\)","",as.character(tax$taxonomy))
#tax = tax[- grep("k__Archaea", tax$taxonomy),]
#tax = tax %>% group_by(taxonomy) %>% summarise_all(list(sum))

tax = tax %>% group_by(taxonomy) %>% summarise_all(list(sum))
phy_dat = tax
tax_levels = strsplit(as.character(tax$taxonomy), split = "; ", fixed = TRUE)
tax_levels = data.frame(t(sapply(tax_levels, function(x) x)))
row.names(tax_levels) = tax$taxonomy
names(tax_levels) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
tax_levels = data.frame(tax_levels)
tax_levels = subset(tax_levels, Domain != "k__ " &  Domain != "k__Eukaryota")
tax_levels$Domain = droplevels(tax_levels$Domain)
tax_levels$Phylum = ifelse(tax_levels$Phylum == "p__ ", "p__Unclassified", as.character(tax_levels$Phylum))
phy_levels = tax_levels
tax_levels$cleaned_names = ifelse(tax_levels$Genus == "g__ ", as.character(tax_levels$Family), as.character(tax_levels$Genus))
tax_levels$cleaned_names = ifelse(tax_levels$cleaned_names == "f__ ", as.character(tax_levels$Order), as.character(tax_levels$cleaned_names))
tax_levels$cleaned_names = ifelse(tax_levels$cleaned_names == "o__ ", as.character(tax_levels$Class), as.character(tax_levels$cleaned_names))
tax_levels$cleaned_names = ifelse(tax_levels$cleaned_names == "c__ ", as.character(tax_levels$Phylum), as.character(tax_levels$cleaned_names))
tax = tax[tax$taxonomy %in% row.names(tax_levels), ]
tax = data.frame(tax)
row.names(tax) = tax_levels[match(tax$taxonomy, row.names(tax_levels)), "cleaned_names"]
tax$taxonomy = NULL
tax = data.frame(t(tax))
tax = tax/rowSums(tax)
rowSums(tax)
head(tax)
#note hadza_seasonal_metadata table came from SRA website metadata; 11358_prep_3753_20210706-110244.txt
#came from qiita (https://qiita.ucsd.edu/study/description/11358#, file downloaded after clicking on "Prep Info") 
#Using this file to map the sample ID (in the form of "11358.109" for example) to the SRR accession number ("SRRXXXX", found under
#the srr_accession column in the qiita file and is what is used in the all_samples_taxonomy_closed_reference.tsv and hadza_seasonal_metadata). 
#This mapping is needed because the SRR metadata was not clear as to which samples
#belonged to which metadata
# otu_df = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/Hadza_seasonal/all_samples_taxonomy_closed_reference.tsv", sep = "\t", header=T, check.names=F))
# with_tax = otu_df
# otu_df$taxonomy = NULL
# otu_df = setNames(data.frame(t(otu_df[,-1])), otu_df[,1])
otu_df = tax
metadata_df = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/Hadza_seasonal/hadza_seasonal_metadata.txt", header=T, check.names=F))
qiita_for_ID = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/Hadza_seasonal/11358_prep_3753_20210706-110244.txt", header=T, check.names=F, colClasses=c('character')))

qiita_for_ID = qiita_for_ID[,c("sample_name","srr_accession")]
sample_name = c(qiita_for_ID$sample_name)
metadata_df$`Library Name` = as.character(metadata_df$`Library Name`)
metadata_df = metadata_df[metadata_df$`Library Name` %in% sample_name,]
metadata_df$sample_name = metadata_df$`Library Name`
metadata_df$`Library Name` = NULL
metadata_df = merge(metadata_df, qiita_for_ID[, c("sample_name", "srr_accession")], by="sample_name")
# use adults only and "seasonality_unique_subjects == 1" which corresponds to individuals who are not repeated across populations
metadata_df = metadata_df[metadata_df$host_life_stage == 'adult',]
metadata_df = metadata_df[metadata_df$seasonality_unique_subjects == 1,]

rownames(metadata_df) = metadata_df$srr_accession
metadata_df$srr_accession = NULL

#number per season
factor_snps <- factor(metadata_df$Season)
summary(factor_snps)

wet_seasons = c("2014-EW","2014-LW")
dry_seasons = c("2013-LD","2014-ED","2014-LD")

# setting to proportional data and 0-1 scale
otu_df = otu_df / rowSums(otu_df)
# set otu_df samples = to metadata samples
otu_df = otu_df[match(row.names(metadata_df), row.names(otu_df)), ]


## filter by phenotype
#separate metadata tables by wet vs dry season
wet_season_df = metadata_df[metadata_df$Season %in% wet_seasons,]
dry_season_df = metadata_df[metadata_df$Season %in% dry_seasons,]

#separate otu table to match separated metadata-by-season tables
otu_df_wet = subset(otu_df, row.names(otu_df) %in% row.names(wet_season_df))
otu_df_dry = subset(otu_df, row.names(otu_df) %in% row.names(dry_season_df))
otu_df_wet = as.data.frame(t(otu_df_wet))
otu_df_dry = as.data.frame(t(otu_df_dry))

eubdry = as.data.frame(t(otu_df_dry))
eubwet = as.data.frame(t(otu_df_wet))
eubcomb = rbind(eubdry,eubwet)
eubcomb$season = "NA"
eubcomb$season[1:60] = "dry"
eubcomb$season[61:75] = "wet"
hadza_box = ggplot(eubcomb, aes(x = season, y = eubcomb$g__Flavonifractor)) + 
  geom_boxplot() +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  xlab("Hadza seasons") +
  ylab("g__Flavonifractor") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(axis.text.x = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14))

#filtering otu-wet samples (less strict)
filtered_wet = otu_df_wet[rowSums(otu_df_wet > 0.001) >= dim(otu_df_wet)[2]*.25, ]
filtered_wet = data.frame(filtered_wet)
filtered_wet = as.data.frame(t(filtered_wet))
#filtering otu-dry samples (less strict)
filtered_dry = otu_df_dry[rowSums(otu_df_dry > 0.001) >= dim(otu_df_dry)[2]*.25, ]
filtered_dry = data.frame(filtered_dry)
filtered_dry = as.data.frame(t(filtered_dry))

#filtering otu-wet samples (more strict)
filtered_strict_wet = otu_df_wet[rowSums(otu_df_wet > 0.001) >= dim(otu_df_wet)[2]*.25, ]
filtered_strict_wet = data.frame(filtered_strict_wet)
filtered_strict_wet = as.data.frame(t(filtered_strict_wet))
#filtering otu-dry samples (more strict)
filtered_strict_dry = otu_df_dry[rowSums(otu_df_dry > 0.001) >= dim(otu_df_dry)[2]*.25, ]
filtered_strict_dry = data.frame(filtered_strict_dry)
filtered_strict_dry = as.data.frame(t(filtered_strict_dry))

# find union of ASVs after apply prevalence/abundance filtering to the subsetted-by-phenotype samples -
# two groups - high fat and herbivore; merge filtered_wet and filtered_dry into filtered
wet_asvs = colnames(filtered_wet)
dry_asvs = colnames(filtered_dry)
union_filt = union(wet_asvs,dry_asvs)
filtered = otu_df[,union_filt]
filtered = filtered[,mixedsort(colnames(filtered))]
# same for the strict filtering; merge filtered_strict_wet and filtered_strict_dry into filtered_strict
wet_asvs_strict = colnames(filtered_strict_wet)
dry_asvs_strict = colnames(filtered_strict_dry)
union_filt_strict = union(wet_asvs_strict,dry_asvs_strict)
filtered_strict = otu_df[,union_filt_strict]
filtered_strict = filtered_strict[,mixedsort(colnames(filtered_strict))]

write.csv(filtered, "hadza_seasonal_species.csv", row.names = TRUE)
write.csv(filtered_strict, "hadza_seasonal_species_strict.csv", row.names = TRUE)

# if going to use 2 seasons for adonis, make new column for either wet or dry season
metadata_df$Binary_season = NA
metadata_df$Binary_season = ifelse(metadata_df$Season %in% wet_seasons, "Wet Season", "Dry Season")
factor_binary_season <- factor(metadata_df$Binary_season)
summary(factor_binary_season)
summary(factor(metadata_df$Season))

write.csv(metadata_df, "hadza_seasonal_metadata_processed.csv", row.names = TRUE)

metadata_df$Season = gsub("-","",metadata_df$Season)
######## adonis section ##########
identical(row.names(filtered), row.names(metadata_df))
filtered = filtered[rownames(filtered) %in% rownames(metadata_df),]
metadata_df_adonis = metadata_df[rownames(metadata_df) %in% rownames(filtered),]
summary(factor(metadata_df_adonis$host_subject_id))


adonis_obj_binary_season = adonis(filtered ~ metadata_df_adonis$Host_Age + metadata_df_adonis$Binary_season, data = metadata_df_adonis, method = "bray",na.rm=TRUE, strata = metadata_df_adonis$host_subject_id)
adonis_obj_binary_season = adonis(filtered ~ metadata_df_adonis$Host_Age + metadata_df_adonis$Binary_season, data = metadata_df_adonis, method = "bray",na.rm=TRUE)

adonis_obj_subseasons = adonis(filtered ~ metadata_df_adonis$Host_Age + metadata_df_adonis$Season, data = metadata_df_adonis, method = "bray",na.rm=TRUE)

adjust_bin = adonis_OmegaSq(adonis_obj_binary_season,partial = TRUE)
adjust_subseason = adonis_OmegaSq(adonis_obj_subseasons,partial = TRUE)
######## maaslin2 section ########
name_list = c('Season','Host_Age')
metadata_df_maaslin = metadata_df[,name_list]
filtered_strict = filtered_strict[rownames(filtered_strict) %in% rownames(metadata_df_maaslin),]
metadata_df_maaslin = metadata_df_maaslin[rownames(metadata_df_maaslin) %in% rownames(filtered_strict),]
metadata_df_maaslin = as.data.frame(metadata_df_maaslin)

#Maaslin2(input_data = filtered_strict, input_metadata = metadata_df_maaslin, output = "hadza_seasonal_maaslin2_output",reference = c("Season,2013LD"))

name_list = c('Binary_season','Host_Age','host_subject_id')
metadata_df_maaslin = metadata_df[,name_list]
filtered_strict = filtered_strict[rownames(filtered_strict) %in% rownames(metadata_df_maaslin),]
metadata_df_maaslin = metadata_df_maaslin[rownames(metadata_df_maaslin) %in% rownames(filtered_strict),]
metadata_df_maaslin = as.data.frame(metadata_df_maaslin)
metadata_df_maaslin$Binary_season = gsub(" ","",metadata_df_maaslin$Binary_season)
Maaslin2(input_data = filtered_strict, input_metadata = metadata_df_maaslin[,1:2], output = "hadza_seasonal_maaslin2_output_binary",reference = c("Binary_season,DrySeason"),max_significance = 1)
