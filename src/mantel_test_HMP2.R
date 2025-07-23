#This script performs the Mantel test (microbiome distance vs. diet distance).HMP2 samples are binned and averaged by week bins;
library(tidyverse)
library(dplyr)
library(stringr)
library(reshape2)

input_dir = file.path("/home","tobynb","diet_test","input","HMP2","Intermediate",fsep = "/")
output_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")

#load microbiome and data files
species <- data.table::fread(paste0(output_dir,"/hmp_species_df.csv"), header=T, check.names=F) %>%
  dplyr::rename(Sample_ID=1) %>%
	arrange(Sample_ID)

#####################################################################################################################################
# Mantel test loop - can loop through several metadata matrices (e.g., if you wanted to compute for different phenotypes separately #
#####################################################################################################################################

input_files <- c(paste0(output_dir,"/df_CD.csv"), paste0(output_dir,"/df_nonIBD.csv"), paste0(output_dir,"/df_UC.csv"))

mantel_res <- list()
options(digits = 5) 

for (input_file in input_files){
	
  #separates non-diet covariates from dietary columns
	metadata <- data.table::fread(input_file, header=T, check.names=F) %>%
	  dplyr::rename(Sample_ID=1) %>%
		select(Sample_ID, ncol(.)-1, ncol(.)-2, ncol(.)-4)

	#groups dietary variables by individual and week number to average the dietary profiles in week bin
	diet_mean <- data.table::fread(input_file, header=T, check.names=F) %>%
	  dplyr::rename(Sample_ID=1) %>%
	  arrange(Sample_ID) %>%
	  reshape2::melt(id=c("Participant_ID", "Sample_ID", "diagnosis", "bin")) %>%
	  dplyr::group_by(bin, Participant_ID, variable) %>%
	  dplyr::summarise(value = mean(c(as.numeric(value)))) %>%
	  ungroup()

	#groups the microbiome profiles by individual and week bins and averages the samples in each week bin
	species_mean <- species %>%
		filter(Sample_ID %in% metadata$Sample_ID) %>%
		reshape2::melt(id=c("Participant_ID", "Sample_ID", "bin")) %>%
	  dplyr::group_by(bin, Participant_ID, variable) %>%
	  dplyr::summarise(value = mean(c(value))) %>%
		ungroup()

	for (i in levels(as.factor(metadata$bin))) {
	
	  #gets list of microbes; excludes bugs that aren't in any samples for that given week (abundance = 0)
		bugs <- species_mean %>%
			filter(bin == i) %>%
			filter(value > 0) %>%
			select(variable) %>%
			distinct()
	
		#gets list of participants that provided samples for each week
		participants <- species_mean %>%
			filter(bin == i) %>%
		  dplyr::group_by(Participant_ID) %>%
		  dplyr::summarise(Total = sum(value)) %>%
			filter(Total > 0)
	
		#subsets microbiome df to include a given week's (averaged) samples and participants
		bugs_df <- species_mean %>%
			filter(bin == i) %>%
			filter(variable %in% bugs$variable) %>%
			filter(Participant_ID %in% participants$Participant_ID) %>%
			reshape2::dcast(Participant_ID ~ variable, value.var="value") %>%
			arrange(Participant_ID) %>%
			column_to_rownames("Participant_ID")
	
		#subsets diet df to include a given week's (averaged) diet profiles and participants
		diet_df <- diet_mean %>%
			filter(bin == i) %>%
			reshape2::dcast(Participant_ID ~ variable, value.var="value") %>%
			arrange(Participant_ID) %>%
			filter(Participant_ID %in% row.names(bugs_df)) %>%
			column_to_rownames("Participant_ID")
	  
		#removes non-diet (age, abx info, sample number) metadata
		diet_df = diet_df[1:(length(diet_df)-3)]
		
		#ensures samples are in the same order
		print(identical(rownames(diet_df),rownames(bugs_df)))
		
		#calculates distance metrics for microbiome and diet dataframes
		bugs.norm.dist <- vegan::vegdist(bugs_df, method="bray")
	  diet.norm.dist <- vegan::vegdist(diet_df, method="manhattan")
		
	  #performs mantel test
		MT <- vegan::mantel(bugs.norm.dist, diet.norm.dist, method="spearman")
	
		#grabs effect size and significance values, places in a dataframe which we'll save and use for visualization
		r <- MT$statistic
		p <- MT$signif
		mt_df <- data.frame(Bin=i, r=r, p=p, n=nrow(bugs_df), group=input_file) %>%
			mutate(group = gsub(".*df_|\\.csv", "", group))
	
		mantel_res <- rbind(mantel_res, mt_df)

	}

}

#formatting data frame which will be used in downstream plotting
mantel_res$Bin <- factor(mantel_res$Bin, levels=c(
  "(0, 7)",
  "(8, 15)",
  "(16, 24)",
  "(25, 33)",
  "(34, 42)",
  "(43, 57)"
)) 

#formatting data frame which will be used in downstream plotting
mantel_res$group <- factor(mantel_res$group, levels=c("nonIBD", "CD", "UC"))
#mantel_res$group = gsub("_new","",mantel_res$group)
ann <- mantel_res %>%
  filter(p < 0.05) %>%
  mutate(signif = case_when( p <= 0.05 & p >= 0.01 ~ "*", p < 0.01 ~ "**" )) %>%
  merge(., mantel_res, by=colnames(mantel_res), all=TRUE)
mantel_res$Bin<- gsub(', ', '_', mantel_res$Bin )
write.table(mantel_res, paste0(output_dir,"/HMP2_mantel_results.csv"), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)