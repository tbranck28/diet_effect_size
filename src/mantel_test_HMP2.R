setwd("/Users/tobynbranck/Documents/Effect_size/analysis/")
library(tidyverse)
library(dplyr)
library(stringr)
options(digits = 5) 
species <- data.table::fread("hmp_species_df.csv", header=T, check.names=F) %>%
  dplyr::rename(Sample_ID=1) %>%
	arrange(Sample_ID)

####################
# loop Mantel test #
####################

input_files <- c("df_CD.csv", "df_nonIBD.csv", "df_UC.csv")

mantel_res <- list()
#options(digits = 30)
for (input_file in input_files){
	
	metadata <- data.table::fread(input_file, header=T, check.names=F) %>%
	  dplyr::rename(Sample_ID=1) %>%
		select(Sample_ID, ncol(.)-1, ncol(.)-2, ncol(.)-4)

	diet_mean <- data.table::fread(input_file, header=T, check.names=F) %>%
	  dplyr::rename(Sample_ID=1) %>%
	  arrange(Sample_ID) %>%
	  reshape2::melt(id=c("Participant_ID", "Sample_ID", "diagnosis", "bin")) %>%
	  dplyr::group_by(bin, Participant_ID, variable) %>%
	  dplyr::summarise(value = mean(c(as.numeric(value)))) %>%
	  ungroup()

	species_mean <- species %>%
		filter(Sample_ID %in% metadata$Sample_ID) %>%
		reshape2::melt(id=c("Participant_ID", "Sample_ID", "bin")) %>%
	  dplyr::group_by(bin, Participant_ID, variable) %>%
	  dplyr::summarise(value = mean(c(value))) %>%
		ungroup()

	for (i in levels(as.factor(metadata$bin))) {
	
		bugs <- species_mean %>%
			filter(bin == i) %>%
			filter(value > 0) %>%
			select(variable) %>%
			distinct()
	
		participants <- species_mean %>%
			filter(bin == i) %>%
		  dplyr::group_by(Participant_ID) %>%
		  dplyr::summarise(Total = sum(value)) %>%
			filter(Total > 0)
	
		bugs_df <- species_mean %>%
			filter(bin == i) %>%
			filter(variable %in% bugs$variable) %>%
			filter(Participant_ID %in% participants$Participant_ID) %>%
			reshape2::dcast(Participant_ID ~ variable, value.var="value") %>%
			arrange(Participant_ID) %>%
			column_to_rownames("Participant_ID")
	
		diet_df <- diet_mean %>%
			filter(bin == i) %>%
			reshape2::dcast(Participant_ID ~ variable, value.var="value") %>%
			arrange(Participant_ID) %>%
			filter(Participant_ID %in% row.names(bugs_df)) %>%
			column_to_rownames("Participant_ID")
	  diet_df = diet_df[1:(length(diet_df)-3)]
		
		bugs.norm.dist <- vegan::vegdist(bugs_df, method="bray")
	  diet.norm.dist <- vegan::vegdist(diet_df, method="manhattan")
		
		MT <- vegan::mantel(bugs.norm.dist, diet.norm.dist, method="spearman")
	
		r <- MT$statistic
		p <- MT$signif
	  p
		mt_df <- data.frame(Bin=i, r=r, p=p, n=nrow(bugs_df), group=input_file) %>%
			mutate(group = gsub("df_|\\.csv", "", group))
	
		mantel_res <- rbind(mantel_res, mt_df)

	}

}

################
# Plot results #
################

# set breaks

my_breaks = c(
	1/max(mantel_res$p),
	(1 /max(mantel_res$p) + 1/min(mantel_res$p))/2,
	1/min(mantel_res$p)
	)

# define labels

my_labels = c(
	round(max(mantel_res$p), 2),
	round((max(mantel_res$p) + min(mantel_res$p))/2, 2),
	round(min(mantel_res$p), 2)
	)

# set colour palette

pal <- c(CD = "#a42a2b", UC = "#ffa400", nonIBD = "#6294eb")

# set the order of the x-axis

mantel_res$Bin <- factor(mantel_res$Bin, levels=c(
  "(0, 7)",
  "(8, 15)",
  "(16, 24)",
  "(25, 33)",
  "(34, 42)",
  "(43, 57)"
)) 

mantel_res$group <- factor(mantel_res$group, levels=c("nonIBD", "CD", "UC"))

ann <- mantel_res %>%
  filter(p < 0.05) %>%
  mutate(signif = case_when( p <= 0.05 & p >= 0.01 ~ "*", p < 0.01 ~ "**" )) %>%
  merge(., mantel_res, by=colnames(mantel_res), all=TRUE)

mantel_res$Bin<- gsub(', ', '_', mantel_res$Bin )
write.table(mantel_res, "/Users/tobynbranck/Documents/Effect_size/analysis/HMP2_mantel_results.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# plot
theme_update(plot.title = element_text(hjust = 0.5))
mantel_p <- ggplot(mantel_res, aes(x=Bin, y=r, colour=group, fill=group)) +
	geom_hline(yintercept=0, linetype="dashed", colour="grey") +
	geom_line(aes(group=group), size=0.75) +
	geom_point(aes(size=n), pch=21, fill="white") +
	geom_point(aes(size=n, alpha=1/p, fill=group), pch=21, colour="black") +
  ggtitle("Mantel r: HMP2 Adults") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() +
	theme(
		aspect.ratio=0.75, 
		axis.text.x=element_text(angle=45, hjust=1)) +
	labs(x="Time bins (weeks)", y="r", alpha="p-value", colour="Group", fill="Group") +
	scale_alpha_continuous(breaks=my_breaks, labels=my_labels) +
	scale_colour_manual(values=pal, labels=c("CD"="CD", "nonIBD"="Non-IBD", "UC"="UC")) +
	scale_fill_manual(values=pal, labels=c("CD"="CD", "nonIBD"="Non-IBD", "UC"="UC")) +
  geom_text(aes(label=p),hjust=-.3, vjust=0, size=3)

ggsave(plot=mantel_p, file="mantel_test_HMP2Adults.png", dpi=300, height=5, width=5/0.75)

mantel_bar <- ggplot(mantel_res, aes(x=Bin, y=r, fill=group)) +
  geom_bar(stat="identity", position="dodge", colour="black") +
  theme_bw() +
  theme(
    axis.text.x=element_text(angle=45, hjust=1),
    legend.position="top") +
  labs(x="Age (days)", y="r", alpha="p-value", colour="Group", fill="Group") +
  scale_alpha_continuous(breaks=my_breaks, labels=my_labels) +
  scale_colour_manual(values=pal, labels=c("CD"="CD", "nonIBD"="Non-IBD", "UC"="UC")) +
  scale_fill_manual(values=pal, labels=c("CD"="CD", "nonIBD"="Non-IBD", "UC"="UC")) +
  geom_text(data=ann, aes(x=Bin, y=r, label=signif), position=position_dodge(width=0.9), size=5)

ggsave(plot=mantel_bar, file="mantel_test_bar_HMP2Adults.png", dpi=300, height=5, width=10)

# combine the plots

mantel_combo <- cowplot::plot_grid(mantel_p, mantel_bar, nrow=1, align="h", axis="tblr", scale=0.95)

ggsave(plot=mantel_combo, file="mantel_test_cowplot_HMP2Adults.png", dpi=300, height=5, width=15)

