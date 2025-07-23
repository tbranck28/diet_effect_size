#This script analyzes one of infant microbiome datasets (Backhed et al.) to
#understand the effects of diet (breastfeeding, intro of solid food) on microbiome
#resulting number of samples after pairing taxonomic profiles with dietary metadata: 263 total (baseline: 80; four: 91; twelve: 92)
#resulting number of individuals: 98

library(RNOmni)
library(ggplot2)
library(grid)
library(vegan)
library(tidyverse)
library(lme4)
library(factoextra)
library(ade4)
library(FactoMineR)
library(dbplyr)
library(plyr)
library(phyloseq)
library(gtools)
library(viridis)
library(colorspace)
library(cowplot)
library(reshape2)
library(Maaslin2)
library(randomForest)
require(caTools)
library(ape)
library(rstatix)
library(gridExtra)
library(ggthemes)
library(ggpubr)
library(pheatmap)

#set directory paths and load taxonomy and dietary profiles
input_dir = file.path("/home","tobynb","diet_test","input","Backhed_infants","Initial",fsep = "/")
output_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")

#load metadata
meta = data.table::fread(paste0(input_dir,"/mmc2.csv"), header=T, check.names = F) #acquired from manuscript supplement
meta_id = data.table::fread(paste0(input_dir,"/studyID_person_association.csv"), header=T, check.names = F) #person ID derived from "Sample alias" found in ENA accession metadata
age_table = as.data.frame(data.table::fread(paste0(input_dir,"/filereport_read_run_PRJEB6456_tsv.txt"), header=T, check.names = F)) #acquired from ENA
age_table = age_table[,c("run_accession","sample_alias")]

#load microbiome profiles and format the table column headers
tax = data.table::fread(paste0(input_dir,"/metaphlan_taxonomic_profiles.tsv"), header=T, check.names = F)
names(tax) <- gsub("Fixed_", "", names(tax), fixed = TRUE)
names(tax) <- gsub("_taxonomic_profile", "", names(tax), fixed = TRUE)
names(tax) <- gsub("# ", "", names(tax), fixed = TRUE)

#grep for species-level taxonomy and format taxonomic names
tax = tax[grep("s__", tax$taxonomy),]
names(tax) = gsub(pattern = "_.*", replacement = "", x = names(tax))
tax$taxonomy <- gsub('^.*?s__', '', tax$taxonomy)
tax$taxonomy = gsub("_", " ", tax$taxonomy)

#divide by 100 to put on 0-1 scale; add taxonomy to rownames
tax <- tax %>%
  mutate(across(-taxonomy, ~ . / 100)) 

#table formatting for downstream operations
bugnames = tax$taxonomy
filttest = tax #storing taxonomy profiles so that unfilterted table can be used later
filttest = data.frame(filttest)
rownames(filttest) = filttest$taxonomy
filttest$taxonomy = NULL
filttest = as.data.frame(t(filttest)) #transposing so that samples are rows

#attaches individual ID to metadata table (note: metadata table did not have the same type of individual ID
#as what is corresponding to sample accession IDs in ENA, so this matches up the metadata's version of individual
#ID which is 'STUDY ID' to the ENA version of the individual ID, found in the 'Sample alias' column).
meta_test = cbind(Individual = meta_id$Individual[ match(meta$`STUDY ID`, meta_id$`STUDY ID`) ],meta)

#renames feeding practice by age columns
names(meta_test)[names(meta_test) == 'feeding practice first week'] = "_B"
names(meta_test)[names(meta_test) == 'feeding practice 4M'] = "_4M"
names(meta_test)[names(meta_test) == 'Any breastfeeding 12 M'] = "_12M"

#recodes antibiotic use to binarize
meta_test$`Antibiotic treatment to infant 0-4M times` = ifelse(meta_test$`Antibiotic treatment to infant 0-4M times`==0, "no", "yes")
meta_test$`Antibiotic treatment to infant 4-12M times` = ifelse(meta_test$`Antibiotic treatment to infant 4-12M times`==0, "no", "yes")

meta_test$feeding_history = paste(meta_test$'_B',meta_test$'_4M')

#remove microbiome samples where sample contains _M in age_table (maternal info)
age_table = age_table[!grepl("_M", age_table$sample_alias),]
accessions_keep = age_table$run_accession
filttest = filttest[rownames(filttest) %in% accessions_keep,]

#adding sample alias to taxonomic profiles table
filttest$Individual <- age_table$sample_alias[match(rownames(filttest), age_table$run_accession)]

#formatting column names in taxonomic profiles dataframe
names(filttest) <- gsub(x = names(filttest), pattern = "\\[", replacement = "_")  
names(filttest) <- gsub(x = names(filttest), pattern = "\\]", replacement = "_")
names(filttest) <- gsub(x = names(filttest), pattern = " ", replacement = "_")
write.table(filttest, paste0(output_dir,"/backhed_species.csv"), sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

#Setting up taxonomy profiles table for prevalence / abundance filtering by age group; 10^-5 in at least 25% of the samples
filttest_baseline = filttest[grep("_B", filttest$Individual), ]
filttest_four = filttest[grep("_4M", filttest$Individual), ]
filttest_twelve = filttest[grep("_12M", filttest$Individual), ]

#isolating sample alias IDs
filttest_baseline_IDs = filttest_baseline$Individual
filttest_four_IDs = filttest_four$Individual
filttest_twelve_IDs = filttest_twelve$Individual

#setting age-specific dataframe rownames as sample IDs
rownames(filttest_baseline) = filttest_baseline_IDs
rownames(filttest_four) = filttest_four_IDs
rownames(filttest_twelve) = filttest_twelve_IDs

###ADDING IN THIS LINE - MAKE SURE IT DOESNT MESS SOMETHING UP DOWNSTREAM
filttest_baseline$Individual = NULL
filttest_four$Individual = NULL
filttest_twelve$Individual = NULL

#transposing for filtering microbes and setting all column data types to be numeric 
filttest_baseline = as.data.frame(t(filttest_baseline))
filttest_baseline <- mutate_all(filttest_baseline, function(x) as.numeric(as.character(x)))
filttest_four = as.data.frame(t(filttest_four))
filttest_four <- mutate_all(filttest_four, function(x) as.numeric(as.character(x)))
filttest_twelve = as.data.frame(t(filttest_twelve))
filttest_twelve <- mutate_all(filttest_twelve, function(x) as.numeric(as.character(x)))

#apply filtering criteria to each age-specific table of taxonomic profiles
filttest_baseline = filttest_baseline[rowSums(filttest_baseline > 0.00001) >= dim(filttest_baseline)[2]*.25, ]
filttest_four = filttest_four[rowSums(filttest_four > 0.00001) >= dim(filttest_four)[2]*.25, ]
filttest_twelve = filttest_twelve[rowSums(filttest_twelve > 0.00001) >= dim(filttest_twelve)[2]*.25, ]

#transposes 
filttest_baseline = as.data.frame(t(filttest_baseline))
filttest_four = as.data.frame(t(filttest_four))
filttest_twelve = as.data.frame(t(filttest_twelve))

#creating data frames of metadata, one for each time point
baseline_meta = meta_test[,c("Individual","_B","GENDER","Delivery mode", "Age at sample Newborn (days)")]
baseline_meta$Individual <- paste0(baseline_meta$Individual, "_B")
baseline_meta = as.data.frame(baseline_meta)

four_meta = meta_test[,c("Individual","_4M","GENDER","Delivery mode","Age at Sampling \n4 M","Antibiotic treatment to infant 0-4M times")]
four_meta$Individual <- paste0(four_meta$Individual, "_4M")
four_meta = as.data.frame(four_meta)

twelve_meta = meta_test[,c("Individual","_12M","GENDER","Delivery mode","Age at Sampling 12M","feeding_history","Antibiotic treatment to infant 4-12M times")]
twelve_meta$Individual <- paste0(twelve_meta$Individual, "_12M")
twelve_meta = as.data.frame(twelve_meta)

#ordering rows by bug name
filttest_baseline = filttest_baseline[order(rownames(filttest_baseline)), ]
filttest_four = filttest_four[order(rownames(filttest_four)), ]
filttest_twelve = filttest_twelve[order(rownames(filttest_twelve)), ]

#dropping rows with missing data
baseline_meta = baseline_meta %>% drop_na()
four_meta = four_meta %>% drop_na()
twelve_meta = twelve_meta %>% drop_na()

#ordering by sample ID
baseline_meta = baseline_meta[order(baseline_meta$Individual), ]
four_meta = four_meta[order(four_meta$Individual), ]
twelve_meta = twelve_meta[order(twelve_meta$Individual), ]

#ensuring only samples with feeding data are kept; ensuring that age-specific taxonomy tables include only the
#samples in the corresponding age-specific metadata
baseline_meta = baseline_meta[!(is.na(baseline_meta$'_B') | baseline_meta$'_B'==""), ]
baseline_meta = baseline_meta[!baseline_meta$'_B' == "Formula feeding", ] #only 1 case
baseline_intersect <- intersect(rownames(filttest_baseline), baseline_meta$Individual)
filttest_baseline <- subset(filttest_baseline, rownames(filttest_baseline) %in% baseline_intersect)
baseline_meta <- subset(baseline_meta, Individual %in% baseline_intersect)

four_meta = four_meta[!(is.na(four_meta$`_4M`) | four_meta$`_4M`==""), ]
four_intersect <- intersect(rownames(filttest_four), four_meta$Individual)
filttest_four <- subset(filttest_four, rownames(filttest_four) %in% four_intersect)
four_meta <- subset(four_meta, Individual %in% four_intersect)

twelve_intersect <- intersect(rownames(filttest_twelve), twelve_meta$Individual)
filttest_twelve <- subset(filttest_twelve, rownames(filttest_twelve) %in% twelve_intersect)
twelve_meta <- subset(twelve_meta, Individual %in% twelve_intersect)

###### heatmap for time ######
#get union of top abundant species for each time point
baseline_vector = c(colnames(filttest_baseline))
four_vector = c(colnames(filttest_four))
twelve_vector = c(colnames(filttest_twelve))
bug_union = union(baseline_vector,union(four_vector,twelve_vector))

#subset original abundance table based on the union of bugs
filttest_for_heat = filttest
rownames(filttest_for_heat) = filttest_for_heat$Individual
names.use <- names(filttest_for_heat)[(names(filttest_for_heat) %in% bug_union)]
filttest_for_heat <- filttest_for_heat[, names.use]
filttest_for_heat <- mutate_all(filttest_for_heat, function(x) as.numeric(x))
#adds age column for grouping by age
filttest_for_heat$Group <- sub("^[^_]*_", "", rownames(filttest_for_heat))

#sets color based on age group
rc <- ifelse(filttest_for_heat$Group == "B", "#440154FF",
             ifelse(filttest_for_heat$Group == "4M", "#31688EFF",
                    ifelse(filttest_for_heat$Group == "12M",
                           "#35B779FF", "#FDE725FF")))
filttest_for_heat$Group <- ordered(filttest_for_heat$Group, levels = c("B", "4M", "12M"))

#formatting and saving taxonomy and metadata profiles 
filttest_for_heat2 = filttest_for_heat
base_comb = baseline_meta
four_comb = four_meta
twelve_comb = twelve_meta
names(base_comb)[2] = 'feeding'
names(four_comb)[2] = 'feeding'
names(twelve_comb)[2] = 'feeding'
names(base_comb)[5] = 'age'
names(four_comb)[5] = 'age'
names(twelve_comb)[5] = 'age'
base_comb$GENDER = NULL
four_comb$GENDER = NULL
twelve_comb$GENDER = NULL
base_comb$`Delivery mode` = NULL
four_comb$`Delivery mode` = NULL
twelve_comb$`Delivery mode` = NULL
four_comb$`Antibiotic treatment to infant 0-4M times`=NULL
twelve_comb$`Antibiotic treatment to infant 4-12M times` = NULL
twelve_comb$feeding_history = NULL
combined_meta = rbind(base_comb, four_comb, twelve_comb)
rownames(combined_meta) = combined_meta$Individual
combined_meta = subset(combined_meta, rownames(combined_meta) %in% rownames(filttest_for_heat2))
filttest_for_heat2 = subset(filttest_for_heat2, rownames(filttest_for_heat2) %in% rownames(combined_meta))
write.csv(filttest_for_heat2, paste0(output_dir,"/backhed_species_for_joint_ord.csv"))
write.csv(combined_meta, paste0(output_dir,"/backhed_meta_joint_ord.csv"))

combined_meta = combined_meta[gtools::mixedorder(rownames(combined_meta)), ]
filttest_for_heat2 = filttest_for_heat2[gtools::mixedorder(rownames(filttest_for_heat2)), ]
identical(rownames(combined_meta),rownames(filttest_for_heat2))

#testing all samples (ages) together - these two tests show that age is a confounder
all_ages = adonis2(filttest_for_heat2[,1:101] ~ combined_meta$feeding, data = combined_meta, method = "bray",na.rm=TRUE,by="margin")
all_ages = adonis2(filttest_for_heat2[,1:101] ~ combined_meta$age + combined_meta$feeding, data = combined_meta, method = "bray",na.rm=TRUE,by="margin")

#creates heatmap (Supplemental Figure 7)
par(lend = 1)
annotation_row = subset(filttest_for_heat, select = c("Group"))
pheatmap(t(filttest_for_heat[1:(length(filttest_for_heat)-1)]), cluster_cols = T, cluster_rows = T, annotation_col = annotation_row, fontsize_col = 6,show_rownames=T,cex=1,angle_col = "45")

######### adonis section ###########
#formatting baseline metadata and running adonis; check to make sure rownames are identical in both tables
rownames(baseline_meta) = baseline_meta$Individual
names(baseline_meta)[5] = 'age'
names(baseline_meta)[2] = 'feeding'
identical(rownames(baseline_meta),rownames(filttest_baseline))
baseline_month_adonis = adonis2(filttest_baseline ~ baseline_meta$GENDER + baseline_meta$`Delivery mode` + baseline_meta$age + baseline_meta$feeding, data = baseline_meta, method = "bray",na.rm=TRUE,by="margin")

#formatting four month group metadata and running adonis; check to make sure rownames are identical in both tables
rownames(four_meta) = four_meta$Individual
names(four_meta)[5] = 'age'
names(four_meta)[2] = 'feeding'
identical(rownames(four_meta),rownames(filttest_four))
four_month_adonis = adonis2(filttest_four ~ four_meta$`Antibiotic treatment to infant 0-4M times` + four_meta$GENDER + four_meta$`Delivery mode` + four_meta$age + four_meta$feeding, data = four_meta, method = "bray",na.rm=TRUE,by="margin")

#formatting twelve month group metadata and running adonis; check to make sure rownames are identical in both tables
rownames(twelve_meta) = twelve_meta$Individual
names(twelve_meta)[5] = 'age'
names(twelve_meta)[2] = 'feeding'
identical(rownames(twelve_meta),rownames(filttest_twelve))
twelve_month_adonis = adonis2(filttest_twelve ~ twelve_meta$`Antibiotic treatment to infant 4-12M times` + twelve_meta$GENDER + twelve_meta$`Delivery mode` + twelve_meta$age + twelve_meta$feeding, data = twelve_meta, method = "bray",na.rm=TRUE,by="margin")
twelve_month_adonis_FeedHist = adonis2(filttest_twelve ~ twelve_meta$`Antibiotic treatment to infant 4-12M times` + twelve_meta$GENDER + twelve_meta$`Delivery mode` + twelve_meta$age + twelve_meta$feeding_history, data = twelve_meta, method = "bray",na.rm=TRUE,by="margin")

combined_adonis_df = as.data.frame(rbind(baseline_month_adonis$aov.tab,four_month_adonis$aov.tab,twelve_month_adonis$aov.tab))

#combining and saving adonis results
adonis_results_combined = as.data.frame(rbind(baseline_month_adonis[4,c(3,5)],four_month_adonis[5,c(3,5)],twelve_month_adonis[5,c(3,5)]))
rownames(adonis_results_combined) = gsub("\\$.*","",rownames(adonis_results_combined))
write.csv(adonis_results_combined, paste0(output_dir,"/backhed_adonis_results.csv"))

######### maaslin2 section ##########
#formatting metadata tables as input into maaslin2
baseline_meta$Individual = NULL
four_meta$Individual= NULL
twelve_meta$Individual = NULL
names(baseline_meta)[3] = 'Delivery_mode'
names(four_meta) <- gsub(x = names(four_meta), pattern = " ", replacement = "_")
names(four_meta) <- gsub(x = names(four_meta), pattern = "-", replacement = "_")
names(twelve_meta) <- gsub(x = names(twelve_meta), pattern = " ", replacement = "_")
names(twelve_meta) <- gsub(x = names(twelve_meta), pattern = "-", replacement = "_")
twelve_meta$feeding_history = NULL

Maaslin2(input_data = filttest_baseline, input_metadata = baseline_meta, output = paste0(output_dir,"/backhed_maaslin2_baseline_8.31"),reference=c("feeding,exclusively breastfeeding"),max_significance = 1.0)
Maaslin2(input_data = filttest_four, input_metadata = four_meta, output = paste0(output_dir,"/backhed_maaslin2_four_8.31"),reference=c("feeding,exclusively breastfeeding"),max_significance = 1.0)
Maaslin2(input_data = filttest_twelve, input_metadata = twelve_meta, output = paste0(output_dir,"/backhed_maaslin2_twelve_8.31"),reference=c("feeding,any breastfeeding"),max_significance = 1.0)

#saving age-specific output files
data_list <- list(
  BackhedBaseline_bugs = filttest_baseline,
  BackhedBaseline_meta = baseline_meta,
  BackhedFour_bugs = filttest_four,
  BackhedFour_meta = four_meta,
  BackhedTwelve_bugs = filttest_twelve,
  BackhedTwelve_meta = twelve_meta
)

lapply(names(data_list), function(name) {
  write.csv(data_list[[name]], file = file.path(output_dir, paste0(name, ".csv")), row.names = TRUE)
})

### bug abundance line plots ###
#combine dataframes from baseline, 4, and 12 timepoints (separated because filtered
# by time point)
common <- intersect(names(filttest_baseline), names(filttest_four))
common = intersect(common,names(filttest_twelve))
all_bugs = rbind(filttest_baseline[,common], filttest_four[,common],filttest_twelve[,common])

common_meta = intersect(names(baseline_meta),names(four_meta))
common_meta = intersect(common_meta,names(twelve_meta))
all_meta = rbind(baseline_meta[,common_meta],four_meta[,common_meta],twelve_meta[,common_meta])
all_meta$time = gsub(".*_", "\\1", rownames(all_meta))
all_bugs$time = gsub(".*_", "\\1", rownames(all_bugs))
all_bugs$person = gsub("_.*", "\\1", rownames(all_bugs))

all_bugs$time = as.character(all_bugs$time)
all_bugs$time = factor(all_bugs$time, levels=c("B", "4M", "12M"))
all_meta$feeding[all_meta$feeding == "exclusively breastfeeding"] = "any breastfeeding"
all_meta$feeding[all_meta$feeding == "mixed feeding"] = "any breastfeeding"
all_meta$feeding[all_meta$feeding == "exclusively formula feeding"] = "no breastfeeding"
all_bugs$feeding = all_meta$feeding

tax = all_bugs %>% 
  group_by(person) %>% 
  mutate(flag2 = if_else(lag(time) == 'B' & time == '4M' &  lag(feeding) != feeding,1, 0, 0))

tax = tax %>% 
  group_by(person) %>% 
  mutate(flag3 = if_else(lag(time) == '4M' & time == '12M' &  lag(feeding) != feeding,1, 0, 0))

tax = tax %>%
  mutate(Feeding_Change = case_when(
    flag2 == 1 ~ "B to 4 change",
    flag3 == 1 ~ "4 to 12 change"
  ))
tax$Feeding_Change = tax$Feeding_Change %>% replace_na("no change")
tax$Feeding_Change = factor(tax$Feeding_Change)
tax = tax %>% ungroup()
tax = as.data.frame(tax)

tax$time = as.character(tax$time)
tax$time[tax$time == "B"] = 0
tax$time[tax$time == "4M"] = 4
tax$time[tax$time == "12M"] = 12
tax$time = as.numeric(tax$time)

split_test = split(tax,tax$time)
df_baseline = as.data.frame(split_test[1])
df_four = as.data.frame(split_test[2])
df_twelve = as.data.frame(split_test[3])
names(df_baseline)[19] = "person"
names(df_four)[19] = "person"
names(df_twelve)[19] = "person"

first_segment = merge(df_baseline, df_four, by="person", all = T)
second_segment = merge(df_four, df_twelve, by="person", all = T)

myplots <- vector('list', ncol(tax[,1:17]))
for (i in 1:ncol(tax[,1:17])) {
  myplots[[i]] <- local({
    i <- i
    bugname = names(tax[i])
    firstcol = gsub(" ", "", paste("X0.",bugname))
    secondcol = gsub(" ", "", paste("X4.",bugname))
    thirdcol = gsub(" ", "", paste("X12.",bugname))
    plt = ggplot(tax, aes(time,tax[,i])) +
      geom_point(aes(color=feeding), size=5, alpha=1) +
      scale_x_continuous(breaks = c(0, 4, 12)) +    
      geom_segment(data = first_segment,
                   aes(x = first_segment$X0.time,
                       y = first_segment[,firstcol],
                       xend = first_segment$X4.time,
                       yend = first_segment[,secondcol],
                       col = first_segment$X4.Feeding_Change)) +
      geom_segment(data = second_segment,
                   aes(x = second_segment$X4.time,
                       y = second_segment[,secondcol],
                       xend = second_segment$X12.time,
                       yend = second_segment[,thirdcol],
                       col = second_segment$X12.Feeding_Change)) +
      scale_color_manual(values=c("black","pink","black","purple","#E1E1E1"))+
      ylab(gsub("_", " ", bugname)) +
      theme_hc() +
      theme(legend.position = "none") +
      theme(axis.title.x=element_blank())
    print(plt)
})
}

cairo_pdf("/backhed_abund_over_time.pdf", width=12, height = 16)
ggarrange(plotlist = myplots, ncol=3,nrow=6)
dev.off()

cairo_pdf("figure3/bacteroides_dorei_AbundOverTime.pdf", width=10, height = 6)
myplots[1]
dev.off()

cairo_pdf("figure3/collins_aero_AbundOverTime.pdf", width=10, height = 6)
myplots[12]
dev.off()

cairo_pdf("figure3/bifido_longum_AbundOverTime.pdf", width=10, height = 6)
myplots[15]
dev.off()

