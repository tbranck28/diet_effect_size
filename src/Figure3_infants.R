library("RNOmni")
library("ggplot2")
library("grid")
library("vegan")
library(tidyverse)
library(lme4)
library(factoextra)
library(ade4)
library(FactoMineR)
library("dbplyr")
library("plyr")
library("phyloseq")
library(gtools)
library("viridis")
library("colorspace")
library("cowplot")
library("reshape2")
library("Maaslin2")
library(randomForest)
require(caTools)
library(ggh4x)
library(ggridges)
library(ggpubr)

setwd("/Users/tobynbranck/Documents/Effect_size/analysis")
#making combined adonis results table from the two infant studies
backhed_adonis = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/backhed_adonis_results.csv", header=T, check.names = F))
murphy_adonis = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/infants_murphy_adonis_results.csv", header=T, check.names = F))
combined_adonis = rbind(backhed_adonis,murphy_adonis)
combined_adonis$FDR = p.adjust(combined_adonis$`Pr(>F)`,method="fdr")
combined_adonis$Study = c("Backhed infants","Backhed infants","Backhed infants","Murphy infants","Murphy infants","Murphy infants")
combined_adonis$'Time point' = c("baseline","four months","twelve months","baseline","three months","twelve months")
write.csv(combined_adonis, "/Users/tobynbranck/Documents/Effect_size/analysis/fig4_adonis_table.csv")

combined_adonis$label = paste0(combined_adonis$Study," ",combined_adonis$`Time point`)

combined_adonis = combined_adonis %>%
  mutate(label = str_remove_all(label, "infants"))
pdf("infants_adonis_plot.pdf", width=10, height = 8)
ggplot(data=combined_adonis, aes(x=label, y=R2,label = ifelse(FDR < 0.25, "*", ""))) +
  geom_bar(stat="identity") +
  geom_text(vjust = 0) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()


#import all maaslin2 results for both infant datasets
backhed_baseline = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/backhed_maaslin2_baseline_8.31/all_results.tsv", header=T, check.names = F))
backhed_baseline$Study = 'B0'
backhed_four = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/backhed_maaslin2_four_8.31/all_results.tsv", header=T, check.names = F))
backhed_four$Study = 'B4'
backhed_twelve = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/backhed_maaslin2_twelve_8.31/all_results.tsv", header=T, check.names = F))
backhed_twelve$Study = 'B12'
murphy_baseline = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/murphy_maaslin2_baseline_8.31/all_results.tsv", header=T, check.names = F))
murphy_baseline$Study = 'M0'
murphy_three = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/murphy_maaslin2_3_month_8.31/all_results.tsv", header=T, check.names = F))
murphy_three$Study = 'M3'
murphy_twelve = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/murphy_maaslin2_twelve_months_8.31/all_results.tsv", header=T, check.names = F))
murphy_twelve$Study = 'M12'

maaslin_results_combined = rbind(backhed_baseline,backhed_four,backhed_twelve,murphy_baseline,murphy_three,murphy_twelve)
maaslin_results_combined$feature<-gsub("_", ".", maaslin_results_combined$feature)

keep = c("feeding","baseline_feeding","three_month_feeding","twelve_month_feeding")
maaslin_results_combined = maaslin_results_combined[maaslin_results_combined$metadata %in% keep, ]
maaslin_results_combined$value = gsub("no breast","no_breastfeeding",maaslin_results_combined$value)
maaslin_results_combined$value = str_replace(maaslin_results_combined$value, "mixed$", "mixed_BF_formula")
maaslin_results_combined$value = str_replace(maaslin_results_combined$value, "mixed_feeding$", "mixed_BF_formula")
maaslin_results_combined$value = gsub("bf","BF",maaslin_results_combined$value)
write.csv(maaslin_results_combined, "/Users/tobynbranck/Documents/Effect_size/analysis/infants_maaslin_combined.csv")

maaslin_results_combined$exp = NA
maaslin_results_combined$exp[1:210] = "Backhed et al."
maaslin_results_combined$exp[211:322] = "Murphy et al."
ggplot(maaslin_results_combined, aes(x=Study, y=coef, fill=value)) + 
  geom_boxplot() +
  facet_wrap(~exp, scale="free") +
  scale_y_continuous(limits=c(-6,6)) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

maaslin_results_combined$feature[maaslin_results_combined$feature=="X.Collinsella..massiliensis"] = "Collinsella.massiliensis"
maaslin_results_combined$Study = factor(maaslin_results_combined$Study, levels=c('B0','B4','B12','M0','M3','M12'))
maaslin_results_combined$feeding_time = paste0(maaslin_results_combined$Study,maaslin_results_combined$value)
anova = aov(abs(coef) ~ feeding_time,data=maaslin_results_combined)
summary(anova)

#values calculated/brought in from Figure2 script. Median, mean, and 95% CI values
#for HMP2 nonIBD RC1 beta coef's (the absolute value of the coefficients)
control_median = 0.19553
control_mean = 0.25439 
control_upper = 0.30261
control_lower = 0.20617

cairo_pdf("figure3/infants_betas.pdf", width=8, height = 6)
ggplot(maaslin_results_combined, aes(x=Study, y=abs(coef), color=value)) + 
  geom_boxplot(lwd=2,outlier.shape = NA) +
  geom_rect(aes(xmin =0.5,xmax = 6.5,ymin = control_lower,ymax = control_upper),
            alpha = 0.1, fill = "grey86") +
  geom_segment(aes(x=0.5,xend=6.5,y=control_median,yend=control_median),color="red") +
  geom_segment(aes(x=0.5,xend=6.5,y=control_lower,yend=control_lower),color="darkgrey",linetype=2) +
  geom_segment(aes(x=0.5,xend=6.5,y=control_upper,yend=control_upper),color="darkgrey",linetype=2) +
  #facet_wrap(~exp,scale="free") +
  scale_y_continuous(limits=c(0,6)) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  theme(axis.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white") ) +
  scale_color_manual(values=c("#9966CC","#FFCCFF","#660066")) +
  theme(legend.position = "none") +
  xlab("") +
  geom_jitter(position=position_jitterdodge())
dev.off()

maaslin_results_combined$feature = gsub("\\."," ",maaslin_results_combined$feature)

pdf("figure3/infants_maaslinheatmap.pdf", width=16, height = 16)
ggplot(maaslin_results_combined, aes(value, feature)) + 
  geom_tile(aes(fill = coef)) + 
  facet_grid(~Study,scale="free",space = "free") + 
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  geom_text(aes(label = ifelse(qval <= .25,"*","")),size = 5) +
  theme(axis.title.y=element_blank(),axis.title.x=element_blank()) +
  theme(axis.text=element_text(size=8)) +
  theme(strip.text.x = element_text(size = 8))
dev.off()

cairo_pdf("figure3/infants_maaslinheatmap_supplement.pdf", width=16, height = 16)
ggplot(maaslin_results_combined, aes(value, reorder(feature,desc(feature)))) + 
  geom_tile(aes(fill = coef)) + 
  facet_grid(~Study,scale="free", space = "free") + 
  theme_bw() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  geom_text(aes(label = ifelse(qval <= .25,"*","")),size = 5) +
  theme(axis.title.y=element_blank(),axis.title.x=element_blank()) +
  theme(axis.text=element_text(size=8)) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(text=element_text(family = "Arial")) +
  theme(axis.text.y = element_text(face = "italic"))
dev.off()

maaslin_results_combined_subset = maaslin_results_combined[maaslin_results_combined$qval < .25, ]
bugs_to_keep = unique(maaslin_results_combined_subset$feature)
maaslin_results_combined_subset2 = subset(maaslin_results_combined, maaslin_results_combined$feature %in% bugs_to_keep)

maaslin_results_combined_subset2$Study = str_sub(maaslin_results_combined_subset2$Study,2)
maaslin_results_combined_subset2$Study <- factor(maaslin_results_combined_subset2$Study, levels = c("0","3","4","12"))
maaslin_results_combined_subset2$value <- factor(maaslin_results_combined_subset2$value, levels = c("mixed_BF_formula","exclusively formula","no_breastfeeding"))

cairo_pdf("figure3/infants_maaslinheatmap_subset.pdf", width=10, height = 8)
ggplot(maaslin_results_combined_subset2, aes(value, reorder(feature,desc(feature)))) + 
  geom_tile(aes(fill = coef)) + 
  facet_nested(.~ exp + Study, scales = "free_x",space = "free_x") +
  #facet_grid(~Study,scale="free", space = "free") + 
  theme_bw() +
  scale_fill_gradient2(low = "blue", mid = "gray92", high = "red") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  geom_text(aes(label = ifelse(qval <= .25,"*","")),size = 7,vjust=.8) +
  theme(axis.title.y=element_blank(),axis.title.x=element_blank()) +
  theme(axis.text=element_text(size=10)) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(text=element_text(family = "Arial")) +
  theme(axis.text.y = element_text(face = "italic"))
dev.off()

#### FORESTPLOT of effect size and significance (comparing infants and adults)

maaslin_results_combined_subset = maaslin_results_combined[maaslin_results_combined$qval < .25,]
bug_list = unique(maaslin_results_combined_subset$feature)
maaslin_results_combined_subset = maaslin_results_combined[maaslin_results_combined$feature %in% bug_list,]
#maaslin_results_combined_subset = maaslin_results_combined

hmp2_maaslin2_sig_results_nonibd <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_diseaseStrat_a_nonIBD/all_results.tsv", header=T, check.names=T))
hmp2_maaslin2_sig_results_nonibd$Study = "HMP2_controls"
mlvs_maaslin2_sig_results <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_maaslin2_PCs/all_results.tsv", header=T, check.names=T))
mlvs_maaslin2_sig_results$Study = "MLVS"
adult_maas_results = rbind(hmp2_maaslin2_sig_results_nonibd,mlvs_maaslin2_sig_results)
keep_vars = c("RC1","RC2")
adult_maas_results = adult_maas_results[adult_maas_results$metadata %in% keep_vars,]

adult_maas_results$feature = gsub("X.","",adult_maas_results$feature)
adult_maas_results$feature = gsub("Collinsella..massiliensis","Collinsella.massiliensis",adult_maas_results$feature)

infant_sig_bugs = unique(maaslin_results_combined_subset$feature)
infant_sig_bugs = gsub(" ",".",infant_sig_bugs)
adult_maas_results = adult_maas_results[adult_maas_results$feature %in% infant_sig_bugs,]
adult_maas_results$genus = gsub("\\..*","",adult_maas_results$feature)

hadza_seas_maaslin2_sig_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hadza_seasonal_maaslin2_output_binary/all_results.tsv", header=T, check.names=T))
hadza_seas_maaslin2_sig_results$Study = "Hadza"
hadza_seas_maaslin2_sig_results = hadza_seas_maaslin2_sig_results[hadza_seas_maaslin2_sig_results$metadata == "Binary_season",]
hadza_seas_maaslin2_sig_results = hadza_seas_maaslin2_sig_results[!grepl("f__", hadza_seas_maaslin2_sig_results$feature),]
hadza_seas_maaslin2_sig_results$genus = gsub(".*\\__","",hadza_seas_maaslin2_sig_results$feature)
hadza_seas_maaslin2_sig_results$genus = gsub("_.*","",hadza_seas_maaslin2_sig_results$genus)

maaslin_results_combined_subset$genus = gsub(" .*","",maaslin_results_combined_subset$feature)
adult_maas_results$feature = gsub("\\."," ",adult_maas_results$feature)
test_merged = merge(adult_maas_results, maaslin_results_combined_subset, by="feature",all = F)
names(test_merged)[23] = "genus"

test_merged2 = merge(test_merged, hadza_seas_maaslin2_sig_results, by="genus",all = T)

#keeping only RC1 for simplicity##
test_merged2 = test_merged2[!test_merged2$metadata.x == "RC2",]
test_merged2 = test_merged2[!is.na(test_merged2$feature.x),]
dat2 = melt(test_merged2, id.vars = "feature.x")

inf_df = test_merged2[,c("feature.x","coef.y","qval.y","Study.y")]
adultW_df = test_merged2[,c("feature.x","coef.x","qval.x","Study.x")]
Hadza_df = test_merged2[,c("feature.x","coef","qval","Study")]

names(inf_df) = c("Bug", "Effect size", "qvalue", "Study")
names(adultW_df) = c("Bug", "Effect size", "qvalue", "Study")
names(Hadza_df) = c("Bug", "Effect size", "qvalue", "Study")

merged_df = unique(rbind(inf_df,adultW_df,Hadza_df))
merged_df = na.omit(merged_df)
merged_df$group = NA
merged_df = merged_df %>% 
  mutate(group = case_when(Study == "B0" ~ "Infants",
                        Study == "B4"  ~ "Infants",
                        Study == "B12" ~ "Infants",
                        Study == "M0" ~ "Infants",
                        Study == "M3" ~ "Infants",
                        Study =="M12" ~ "Infants",
                        Study == "HMP2_controls" ~ "Adults - Westernized diet",
                        Study == "MLVS" ~ "Adults - Westernized diet",
                        Study == "Hadza" ~ "Adults - Traditional diet"))

ggplot(merged_df, aes(x=Bug, y= merged_df$`Effect size`))+
  geom_point(aes(shape = Study,color=group),size=3) +
  geom_line() +
  scale_shape_manual(values = 0:9) +
  geom_hline(yintercept = 0, linetype=2, color="grey") +
  coord_flip() +
  xlab("") +
  ylab('Effect size') + 
  theme_classic()

merged_df$Bug = gsub("\\."," ",merged_df$Bug)

cairo_pdf("figure3/effect_size_FOREST_full.pdf",width=8,height=10)
ggplot(merged_df, aes(x=reorder(Bug,desc(Bug)), y= abs(merged_df$`Effect size`)))+
  geom_point(aes(shape = Study,color=group,size=-log10(qvalue)),stroke = 1) +
  scale_size(range=c(2,8)) +
  geom_line(color="grey",alpha=0.5) +
  scale_shape_manual(values = 0:9) +
  coord_flip() +
  xlab("") +
  ylab('Effect size') +
  scale_color_manual(values=c("darkgreen","blue4","darkred")) +
  theme_classic() +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(text = element_text(family = "Arial"))
dev.off()

effectSize_forest = ggplot(merged_df, aes(x=reorder(Bug,desc(Bug)), y= abs(merged_df$`Effect size`)))+
  geom_point(aes(shape = Study,color=group,size=-log10(qvalue)),stroke = 1) +
  scale_size(range=c(2,8)) +
  geom_line(color="grey",alpha=0.5) +
  scale_shape_manual(values = 0:9) +
  coord_flip() +
  xlab("") +
  ylab('Effect size') +
  scale_color_manual(values=c("darkgreen","blue4","darkred")) +
  theme_classic() +
  theme(axis.text.y = element_text(face = "italic"))

ES_legend = get_legend(effectSize_forest)
cairo_pdf("figure3/effect_size_FOREST_Legend.pdf",width=3,height=6)
as_ggplot(ES_legend)
dev.off()

cairo_pdf("figure3/effect_size_FOREST_subset_nolegend.pdf",width=6,height=4)
ggplot(merged_df, aes(x=reorder(Bug,desc(Bug)), y= abs(merged_df$`Effect size`)))+
  geom_point(aes(shape = Study,color=group,size=-log10(qvalue)),stroke = 1) +
  scale_size(range=c(2,8)) +
  geom_line(color="grey",alpha=0.5) +
  scale_shape_manual(values = 0:9) +
  coord_flip() +
  xlab("") +
  ylab('Effect size') +
  scale_color_manual(values=c("darkgreen","blue4","darkred")) +
  theme_classic() +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(legend.position="none")
dev.off()

cairo_pdf("qval_FOREST.pdf",width=7,height=5)
ggplot(merged_df, aes(x=reorder(Bug,desc(Bug)), y= -log10(merged_df$qvalue)))+
  geom_point(aes(shape = Study,color=group),size=3,stroke = 1) +
  geom_line(color="grey",alpha=0.5) +
  scale_shape_manual(values = 0:9) +
  coord_flip() +
  xlab("") +
  ylab('qvalue') +
  scale_color_manual(values=c("darkgreen","blue","red")) +
  theme_classic()
dev.off()

#### EXTRA CODE ####
####plotting bug-specific relative abundance plots####
backhed_baseline_boxplots = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/BackhedBaseline_bugs.csv", header=T, check.names = F))
backhed_four_boxplots = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/BackhedFour_bugs.csv", header=T, check.names = F))
backhed_twelve_boxplots = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/BackhedTwelve_bugs.csv", header=T, check.names = F))
murphy_baseline_boxplots = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/MurphyBaseline_bugs.csv", header=T, check.names = F))
murphy_three_boxplots = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/MurphyThreeMonth_bugs.csv", header=T, check.names = F))
murphy_twelve_boxplots = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/MurphyTwelveMonth_bugs.csv", header=T, check.names = F))

backhed_baseline_meta = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/BackhedBaseline_meta.csv", header=T, check.names = F))
backhed_four_meta = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/BackhedFour_meta.csv", header=T, check.names = F))
backhed_twelve_meta = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/BackhedTwelve_meta.csv", header=T, check.names = F))
murphy_baseline_meta = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/MurphyBaseline_meta.csv", header=T, check.names = F))
murphy_three_meta = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/MurphyThreeMonth_meta.csv", header=T, check.names = F))
murphy_twelve_meta = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/MurphyTwelveMonth_meta.csv", header=T, check.names = F))

backhed_baseline_boxplots$study = 'Backhed et al.'
murphy_baseline_boxplots$study = 'Murphy et al.'

names(backhed_baseline_boxplots) = gsub("_"," ",names(backhed_baseline_boxplots))
baseline_df = rbind.fill(backhed_baseline_boxplots,murphy_baseline_boxplots)
baseline_df$study[1:80] = "Backhed et al."
baseline_df$study[81:166] = "Murphy et al"
baseline_df$feeding = NA
baseline_df$feeding[1:80] = backhed_baseline_meta$feeding
baseline_df$feeding[81:166] = murphy_baseline_meta$baseline_feeding
baseline_df$feeding[baseline_df$feeding=='exclusively breastfeeding'] = 'exclusively breastfed'
baseline_df$time = "< 1 week"
pdf("figure3/infants_staph_epi_boxplots.pdf", width=6, height = 4)
ggplot(baseline_df, aes(x = feeding, y = `Staphylococcus epidermidis`)) +
  geom_boxplot(size=1) +
  geom_jitter() +
  facet_grid(.~study,scales = "free_x",space = "free_x") +
  theme_classic() + 
  #theme(axis.text.x = element_text(angle=45,hjust = 1)) +
  xlab(NULL)
dev.off()

baseline_df$feeding[baseline_df$feeding=="mixed_feeding"] = "mixed feeding"

staph_epi_df = as.data.frame(baseline_df[,"Staphylococcus epidermidis"])
rownames(staph_epi_df) = baseline_df$V1
staph_epi_df$feeding = baseline_df$feeding
names(staph_epi_df)[1] = "Staphylococcus epidermidis"

staph_epi_df$study = NA
staph_epi_df$study[1:80] = "Backhed et al."
staph_epi_df$study[81:166] = "Murphy et al."

staph_epi_df$concat = paste(staph_epi_df$study,staph_epi_df$feeding)
zeros = aggregate(staph_epi_df$`Staphylococcus epidermidis` ~ staph_epi_df$concat, staph_epi_df,function(x) sum(x==0,na.rm=TRUE))
group_totals = aggregate(staph_epi_df$concat,by=list(staph_epi_df$concat),FUN=length)
zeros_df = cbind(zeros,group_totals$x)
names(zeros_df) = c("group","num_zeros","group_tot")
zeros_df = transform(zeros_df, proportion = num_zeros / group_tot)
zeros_df$proportion = zeros_df$proportion*100
zeros_df$study = c("Backhed et al.","Backhed et al.","Murphy et al.","Murphy et al.")

staph_epi_df[staph_epi_df==0] = NA
staph_epi_df = staph_epi_df[complete.cases(staph_epi_df),]

staph_epi_df_backhed = subset(staph_epi_df, study == "Backhed et al.")
staph_epi_df_murphy = subset(staph_epi_df, study == "Murphy et al.")

tmp_back = subset(staph_epi_df_backhed, feeding == "exclusively breastfed")
median_value_back = median(tmp_back$`Staphylococcus epidermidis`)
staph_epi_df$`Staphylococcus epidermidis`[1:63] = log10(staph_epi_df$`Staphylococcus epidermidis`[1:63]/median_value_back)

tmp_murph = subset(staph_epi_df_murphy, feeding == "exclusively breastfed")
median_value_murph = median(tmp_murph$`Staphylococcus epidermidis`)
staph_epi_df$`Staphylococcus epidermidis`[64:140] = log10(staph_epi_df$`Staphylococcus epidermidis`[64:140]/median_value_murph)

density_s_epi_joint = ggplot(staph_epi_df, aes(x = as.numeric(staph_epi_df$`Staphylococcus epidermidis`), y= factor(concat,levels=c("Murphy et al. exclusively breastfed","Murphy et al. mixed feeding","Backhed et al. exclusively breastfed","Backhed et al. mixed feeding")), fill=factor(concat,levels=c("Murphy et al. exclusively breastfed","Murphy et al. mixed feeding","Backhed et al. exclusively breastfed","Backhed et al. mixed feeding")))) + 
  geom_density_ridges(alpha = 0.8,scale=3) + 
  theme_classic(base_size = 14) + 
  labs(fill = "Study", title = "Staphylococcus epidermidis") + 
  ylab("Density") + 
  xlab("Relative abundance") + 
  theme(plot.title = element_text(face = "italic")) + 
  theme(legend.position = "NONE") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  + 
  scale_fill_manual(values=c("gray90","#999999", "#FFCCCC","#996666"))

zeroplots_s_epi_joint = ggplot(zeros_df, aes(x = factor(group,levels=c("Backhed et al. mixed feeding", "Backhed et al. exclusively breastfed", "Murphy et al. mixed feeding","Murphy et al. exclusively breastfed")), y = proportion, fill = factor(group,levels=c("Backhed et al. mixed feeding", "Backhed et al. exclusively breastfed", "Murphy et al. mixed feeding","Murphy et al. exclusively breastfed")))) + 
  geom_bar(position = position_dodge(0.1), stat = "identity") + theme_classic(base_size = 14)+ 
  ylab("Zeros (%)") + 
  xlab("") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())  + 
  theme(legend.position = c(0.65,0.88),legend.text=element_text(size=8),legend.title = element_blank()) +
  ylim(0,100) +
  scale_fill_manual(values=c("#996666","#FFCCCC","#999999","gray90"))
  
#grid.arrange(zeroplots_s_epi_joint, density_s_epi_joint, ncol = 2, nrow = 1, widths = c(2, 3))
gA = ggplotGrob(zeroplots_s_epi_joint)
gB = ggplotGrob(density_s_epi_joint)
grid::grid.newpage()

pdf("figure3/infants_staph_epi_density.pdf", width=8, height = 4)
grid::grid.draw(cbind(gA,gB))
dev.off()

#### Four month bug specific abundances #########
backhed_four_boxplots$study = 'Backhed et al.'
murphy_three_boxplots$study = 'Murphy et al.'

names(backhed_four_boxplots) = gsub("_"," ",names(backhed_four_boxplots))
mid_df = rbind.fill(backhed_four_boxplots,murphy_three_boxplots)
mid_df$study[1:91] = "Backhed et al."
mid_df$study[92:165] = "Murphy et al"
mid_df$feeding = NA
mid_df$feeding[1:91] = backhed_four_meta$feeding
mid_df$feeding[92:165] = murphy_three_meta$three_month_feeding
mid_df$feeding[mid_df$feeding=='exclusively breastfeeding'] = 'exclusively breastfed'
mid_df$time = "mid (3 and 4 months)"

pdf("figure3/infants_bifido_longum_boxplots.pdf", width=6, height = 4)
ggplot(mid_df, aes(x = factor(feeding,levels=c("exclusively breastfed","mixed feeding","exclusively formula feeding","no_breastfeeding")), y = `Bifidobacterium longum`)) +
  geom_boxplot(size=1) +
  geom_jitter() +
  facet_grid(.~study,scales = "free_x",space = "free_x") +
  theme_classic() + 
  #theme(axis.text.x = element_text(angle=45,hjust = 1)) +
  xlab(NULL)
dev.off()

pdf("figure3/infants_gpam_boxplots.pdf", width=6, height = 4)
ggplot(mid_df, aes(x = factor(feeding,levels=c("exclusively breastfed","mixed feeding","exclusively formula feeding","no_breastfeeding")), y = `Gordonibacter pamelaeae`)) +
  geom_boxplot(size=1) +
  geom_jitter() +
  facet_grid(.~study,scales = "free_x",space = "free_x") +
  theme_classic() + 
  #theme(axis.text.x = element_text(angle=45,hjust = 1)) +
  xlab(NULL)
dev.off()

mid_df$feeding[mid_df$feeding=="mixed_feeding"] = "mixed feeding"

bifido_longum_df = as.data.frame(mid_df[,"Bifidobacterium longum"])
rownames(bifido_longum_df) = mid_df$V1
bifido_longum_df$feeding = mid_df$feeding
names(bifido_longum_df)[1] = "Bifidobacterium longum"

bifido_longum_df$study = NA
bifido_longum_df$study[1:91] = "Backhed et al."
bifido_longum_df$study[92:165] = "Murphy et al."

bifido_longum_df$concat = paste(bifido_longum_df$study,bifido_longum_df$feeding)
zeros = aggregate(bifido_longum_df$`Bifidobacterium longum` ~ bifido_longum_df$concat, bifido_longum_df,function(x) sum(x==0,na.rm=TRUE))
group_totals = aggregate(bifido_longum_df$concat,by=list(bifido_longum_df$concat),FUN=length)
zeros_df = cbind(zeros,group_totals$x)
names(zeros_df) = c("group","num_zeros","group_tot")
zeros_df = transform(zeros_df, proportion = num_zeros / group_tot)
zeros_df$proportion = zeros_df$proportion*100
zeros_df$study = c("Backhed et al.","Backhed et al.","Backhed et al.","Murphy et al.","Murphy et al.","Murphy et al.")

bifido_longum_df[bifido_longum_df==0] = NA
bifido_longum_df = bifido_longum_df[complete.cases(bifido_longum_df),]

bifido_df_backhed = subset(bifido_longum_df, study == "Backhed et al.")
bifido_df_murphy = subset(bifido_longum_df, study == "Murphy et al.")

tmp_back_bifido = subset(bifido_df_backhed, feeding == "exclusively breastfed")
median_value_back_bifido = median(tmp_back_bifido$`Bifidobacterium longum`)
bifido_longum_df$`Bifidobacterium longum`[1:82] = log10(bifido_longum_df$`Bifidobacterium longum`[1:82]/median_value_back_bifido)

tmp_murph_bifido = subset(bifido_df_murphy, feeding == "exclusively breastfed")
median_value_murph_bifido = median(tmp_murph_bifido$`Bifidobacterium longum`)
bifido_longum_df$`Bifidobacterium longum`[83:131] = log10(bifido_longum_df$`Bifidobacterium longum`[83:131]/median_value_murph_bifido)

density_bifido_longum_joint = ggplot(bifido_longum_df, aes(x = as.numeric(bifido_longum_df$`Bifidobacterium longum`), y= factor(concat,levels=c("Murphy et al. exclusively breastfed","Murphy et al. mixed feeding","Murphy et al. no_breastfeeding","Backhed et al. exclusively breastfed","Backhed et al. mixed feeding","Backhed et al. exclusively formula feeding")), fill=factor(concat,levels=c("Murphy et al. exclusively breastfed","Murphy et al. mixed feeding","Murphy et al. no_breastfeeding","Backhed et al. exclusively breastfed","Backhed et al. mixed feeding","Backhed et al. exclusively formula feeding")))) + 
  geom_density_ridges(alpha = 0.8,scale=3) + 
  theme_classic(base_size = 14) + 
  labs(fill = "Study", title = "Bifidobacterium longum") + 
  ylab("Density") + 
  xlab("Relative abundance") + 
  theme(plot.title = element_text(face = "italic")) + 
  theme(legend.position = "NONE") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  + 
  scale_fill_manual(values=c("gray90","#999999","#333333","#FFCCCC","#996666","#660000"))

zeroplots_bifido_longum_joint = ggplot(zeros_df, aes(x = factor(group,levels=c("Backhed et al. exclusively formula feeding","Backhed et al. mixed feeding", "Backhed et al. exclusively breastfed", "Murphy et al. no_breastfeeding","Murphy et al. mixed feeding","Murphy et al. exclusively breastfed")), y = proportion, fill = factor(group,levels=c("Backhed et al. exclusively formula feeding","Backhed et al. mixed feeding", "Backhed et al. exclusively breastfed", "Murphy et al. no_breastfeeding","Murphy et al. mixed feeding","Murphy et al. exclusively breastfed")))) + 
  geom_bar(position = position_dodge(0.1), stat = "identity") + theme_classic(base_size = 14)+ 
  ylab("Zeros (%)") + 
  xlab("") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())  + 
  theme(legend.position = c(0.65,0.9),legend.text=element_text(size=8),legend.title = element_blank()) +
  ylim(0,100) +
  #scale_fill_manual(values=c("gray90","#999999","#333333","#FFCCCC","#996666","#660000"))
  scale_fill_manual(values=c("#660000","#996666","#FFCCCC","#333333","#999999","gray90"))

gA_bifido_longum = ggplotGrob(zeroplots_bifido_longum_joint)
gB_bifido_longum = ggplotGrob(density_bifido_longum_joint)
grid::grid.newpage()

pdf("figure3/infants_bifido_longum_density.pdf", width=8, height = 4)
grid::grid.draw(cbind(gA_bifido_longum,gB_bifido_longum))
dev.off()

###gordonibacter pamelaeae (mid time points)
gpam_df = as.data.frame(mid_df[,"Gordonibacter pamelaeae"])
rownames(gpam_df) = mid_df$V1
gpam_df$feeding = mid_df$feeding
names(gpam_df)[1] = "Gordonibacter pamelaeae"

gpam_df$study = NA
gpam_df$study[1:91] = "Backhed et al."
gpam_df$study[92:165] = "Murphy et al."

gpam_df$concat = paste(gpam_df$study,gpam_df$feeding)
zeros = aggregate(gpam_df$`Gordonibacter pamelaeae` ~ gpam_df$concat, gpam_df,function(x) sum(x==0,na.rm=TRUE))
group_totals = aggregate(gpam_df$concat,by=list(gpam_df$concat),FUN=length)
zeros_df = cbind(zeros,group_totals$x)
names(zeros_df) = c("group","num_zeros","group_tot")
zeros_df = transform(zeros_df, proportion = num_zeros / group_tot)
zeros_df$proportion = zeros_df$proportion*100
zeros_df$study = c("Backhed et al.","Backhed et al.","Backhed et al.","Murphy et al.","Murphy et al.","Murphy et al.")

gpam_df[gpam_df==0] = NA
gpam_df = gpam_df[complete.cases(gpam_df),]

gpam_df_backhed = subset(gpam_df, study == "Backhed et al.")
gpam_df_murphy = subset(gpam_df, study == "Murphy et al.")

tmp_back_gpam = subset(gpam_df_backhed, feeding == "exclusively breastfed")
median_value_back_gpam = median(tmp_back_gpam$`Gordonibacter pamelaeae`)
gpam_df$`Gordonibacter pamelaeae`[1:56] = log10(gpam_df$`Gordonibacter pamelaeae`[1:56]/median_value_back_gpam)

tmp_murph_gpam = subset(gpam_df_murphy, feeding == "exclusively breastfed")
median_value_murph_gpam = median(tmp_murph_gpam$`Gordonibacter pamelaeae`)
gpam_df$`Gordonibacter pamelaeae`[57:86] = log10(gpam_df$`Gordonibacter pamelaeae`[57:86]/median_value_murph_gpam)

density_gpam_joint = ggplot(gpam_df, aes(x = as.numeric(gpam_df$`Gordonibacter pamelaeae`), y= factor(concat,levels=c("Murphy et al. exclusively breastfed","Murphy et al. mixed feeding","Murphy et al. no_breastfeeding","Backhed et al. exclusively breastfed","Backhed et al. mixed feeding","Backhed et al. exclusively formula feeding")), fill=factor(concat,levels=c("Murphy et al. exclusively breastfed","Murphy et al. mixed feeding","Murphy et al. no_breastfeeding","Backhed et al. exclusively breastfed","Backhed et al. mixed feeding","Backhed et al. exclusively formula feeding")))) + 
  geom_density_ridges(alpha = 0.8,scale=3) + 
  theme_classic(base_size = 14) + 
  labs(fill = "Study", title = "Gordonibacter pamelaeae") + 
  ylab("Density") + 
  xlab("Relative abundance") + 
  theme(plot.title = element_text(face = "italic")) + 
  theme(legend.position = "NONE") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  + 
  scale_fill_manual(values=c("gray90","#999999","#333333","#FFCCCC","#996666","#660000"))

zeroplots_gpam_joint = ggplot(zeros_df, aes(x = factor(group,levels=c("Backhed et al. exclusively formula feeding","Backhed et al. mixed feeding", "Backhed et al. exclusively breastfed", "Murphy et al. no_breastfeeding","Murphy et al. mixed feeding","Murphy et al. exclusively breastfed")), y = proportion, fill = factor(group,levels=c("Backhed et al. exclusively formula feeding","Backhed et al. mixed feeding", "Backhed et al. exclusively breastfed", "Murphy et al. no_breastfeeding","Murphy et al. mixed feeding","Murphy et al. exclusively breastfed")))) + 
  geom_bar(position = position_dodge(0.1), stat = "identity") + theme_classic(base_size = 14)+ 
  ylab("Zeros (%)") + 
  xlab("") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())  + 
  theme(legend.position = c(0.65,0.9),legend.text=element_text(size=8),legend.title = element_blank()) +
  ylim(0,100) +
  #scale_fill_manual(values=c("gray90","#999999","#333333","#FFCCCC","#996666","#660000"))
  scale_fill_manual(values=c("#660000","#996666","#FFCCCC","#333333","#999999","gray90"))

gA_gpam = ggplotGrob(zeroplots_gpam_joint)
gB_gpam = ggplotGrob(density_gpam_joint)
grid::grid.newpage()

pdf("figure3/infants_gpam_density.pdf", width=8, height = 4)
grid::grid.draw(cbind(gA_gpam,gB_gpam))
dev.off()

### twelve months ###
backhed_twelve_boxplots$study = 'Backhed et al.'
murphy_twelve_boxplots$study = 'Murphy et al.'

names(backhed_twelve_boxplots) = gsub("_"," ",names(backhed_twelve_boxplots))
twelve_df = rbind.fill(backhed_twelve_boxplots,murphy_twelve_boxplots)
twelve_df$study[1:92] = "Backhed et al."
twelve_df$study[93:183] = "Murphy et al"
twelve_df$feeding = NA
twelve_df$feeding[1:92] = backhed_twelve_meta$feeding
twelve_df$feeding[93:183] = murphy_twelve_meta$twelve_month_feeding
twelve_df$feeding[twelve_df$feeding=='exclusively breastfeeding'] = 'exclusively breastfed'
twelve_df$time = "twelve months"
twelve_df$feeding[twelve_df$feeding=="no_breastfeeding"] = "no breastfeeding"

pdf("figure3/infants_blautiacoccoides_boxplots.pdf", width=6, height = 4)
ggplot(twelve_df, aes(x = factor(feeding,levels=c("any breastfeeding","mixed_bf_solid_formula","no breastfeeding")), y = `Blautia coccoides`)) +
  geom_boxplot(size=1) +
  geom_jitter() +
  facet_grid(.~study,scales = "free_x",space = "free_x") +
  theme_classic() + 
  xlab(NULL)
dev.off()

twelve_df$feeding[twelve_df$feeding=="mixed_feeding"] = "mixed feeding"

blautiacoccoides_df = as.data.frame(twelve_df[,"Blautia coccoides"])
rownames(blautiacoccoides_df) = twelve_df$V1
blautiacoccoides_df$feeding = twelve_df$feeding
names(blautiacoccoides_df)[1] = "Blautia coccoides"

blautiacoccoides_df$study = NA
blautiacoccoides_df$study[1:92] = "Backhed et al."
blautiacoccoides_df$study[93:183] = "Murphy et al."

blautiacoccoides_df$concat = paste(blautiacoccoides_df$study,blautiacoccoides_df$feeding)
zeros = aggregate(blautiacoccoides_df$`Blautia coccoides` ~ blautiacoccoides_df$concat, blautiacoccoides_df,function(x) sum(x==0,na.rm=TRUE))
group_totals = aggregate(blautiacoccoides_df$concat,by=list(blautiacoccoides_df$concat),FUN=length)
zeros_df = cbind(zeros,group_totals$x)
names(zeros_df) = c("group","num_zeros","group_tot")
zeros_df = transform(zeros_df, proportion = num_zeros / group_tot)
zeros_df$proportion = zeros_df$proportion*100
zeros_df$study = c("Backhed et al.","Backhed et al.","Murphy et al.","Murphy et al.")

blautiacoccoides_df[blautiacoccoides_df==0] = NA
blautiacoccoides_df = blautiacoccoides_df[complete.cases(blautiacoccoides_df),]

blautia_backhed = subset(blautiacoccoides_df, study == "Backhed et al.")
blautia_murphy = subset(blautiacoccoides_df, study == "Murphy et al.")

#tmp_back_blautia = subset(blautia_backhed, feeding == "any breastfeeding")
#median_value_back_blautia = 0.00000000001

tmp_murph_blautia = subset(blautia_murphy, feeding == "mixed_bf_solid_formula")
median_value_murph_blautia = median(tmp_murph_blautia$`Blautia coccoides`)
blautiacoccoides_df$`Blautia coccoides`[33:69] = log10(blautiacoccoides_df$`Blautia coccoides`[33:69]/median_value_murph_blautia)

blautiacoccoides_df$`Blautia coccoides`[1:32] = log10(blautiacoccoides_df$`Blautia coccoides`[1:32]/median_value_murph_blautia)

density_blautiacoccoides_joint = ggplot(blautiacoccoides_df, aes(x = as.numeric(blautiacoccoides_df$`Blautia coccoides`), y= factor(concat,levels=c("Murphy et al. mixed_bf_solid_formula","Murphy et al. no_breastfeeding","Backhed et al. any breastfeeding","Backhed et al. no breastfeeding")), fill=factor(concat,levels=c("Murphy et al. mixed_bf_solid_formula","Murphy et al. no_breastfeeding","Backhed et al. any breastfeeding","Backhed et al. no breastfeeding")))) + 
  geom_density_ridges(alpha = 0.8,scale=3) + 
  theme_classic(base_size = 14) + 
  labs(fill = "Study", title = "Blautia coccoides") + 
  ylab("Density") + 
  xlab("Relative abundance") + 
  theme(plot.title = element_text(face = "italic")) + 
  theme(legend.position = "NONE") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  + 
  scale_fill_manual(values=c("#999999","#333333","#336666"))

zeroplots_blautiacoccoides_joint = ggplot(zeros_df, aes(x = factor(group,levels=c("Backhed et al. no breastfeeding","Backhed et al. any breastfeeding", "Murphy et al. no_breastfeeding","Murphy et al. mixed_bf_solid_formula")), y = proportion, fill = factor(group,levels=c("Backhed et al. no breastfeeding","Backhed et al. any breastfeeding", "Murphy et al. no_breastfeeding","Murphy et al. mixed_bf_solid_formula")))) + 
  geom_bar(position = position_dodge(0.1), stat = "identity") + theme_classic(base_size = 14)+ 
  ylab("Zeros (%)") + 
  xlab("") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())  + 
  theme(legend.position = c(0.85,0.9),legend.text=element_text(size=6),legend.title = element_blank()) +
  ylim(0,100) +
  scale_fill_manual(values=c("#336666","#99CCCC","#333333","#999999"))

gA_blautiacoccoides = ggplotGrob(zeroplots_blautiacoccoides_joint)
gB_blautiacoccoides = ggplotGrob(density_blautiacoccoides_joint)
grid::grid.newpage()

pdf("figure3/infants_blautiacoccoides_density.pdf", width=8, height = 4)
grid::grid.draw(cbind(gA_blautiacoccoides,gB_blautiacoccoides))
dev.off()

### spearman for effect sizes measured in both datasets ####
combined_correlation = maaslin_results_combined
combined_lists = split(combined_correlation,combined_correlation$Study)
b0 = as.data.frame(combined_lists[1])
b4 = as.data.frame(combined_lists[2])
b12 = as.data.frame(combined_lists[3])
m0 = as.data.frame(combined_lists[4])
m3 = as.data.frame(combined_lists[5])
m12 = as.data.frame(combined_lists[6])

names(b0) = gsub("B0.","",names(b0))
names(b4) = gsub("B4.","",names(b4))
names(b12) = gsub("B12.","",names(b12))
names(m0) = gsub("M0.","",names(m0))
names(m3) = gsub("M3.","",names(m3))
names(m12) = gsub("M12.","",names(m12))

#for the purposes of this, renaming exclusively formula to no breastfeeding to 
#be able to compare
b4$value[b4$value=='exclusively formula'] = 'no_breastfeeding'

baseline_corr = merge(b0, m0, by=c("feature","value"))
mid_corr = merge(b4, m3, by=c("feature","value"))
late_corr = merge(b12,m12, by=c("feature","value"))

baseline_spear = cor(baseline_corr$coef.x,baseline_corr$coef.y,method="spearman")
mid_spear = cor(mid_corr$coef.x,mid_corr$coef.y,method = "spearman")
late_spear = cor(late_corr$coef.x,late_corr$coef.y,method = "spearman")

basecorr = ggplot(baseline_corr,aes(x=coef.x,y=coef.y)) +
  ggtitle("< 1 week") +
  geom_point() +
  stat_smooth(method = 'lm') +
  geom_text(x=1.5, y=-3, label="spearman r = 0.22") +
  xlab(NULL) + ylab(NULL)
midcorr = ggplot(mid_corr,aes(x=coef.x,y=coef.y)) +
  ggtitle("3-4 months") +
  geom_point() +
  stat_smooth(method = 'lm') +
  geom_text(x=2.2, y=-5, label="spearman r = -0.03") +
  xlab(NULL) + ylab(NULL)
latecorr = ggplot(late_corr,aes(x=coef.x,y=coef.y)) +
  ggtitle("Twelve months") +
  geom_point() +
  stat_smooth(method = 'lm') +
  geom_text(x=2, y=-2.65, label="spearman r = 0.32") +
  xlab(NULL) + ylab(NULL)


fig = ggarrange(basecorr, midcorr, latecorr, ncol = 3, nrow = 1,labels=NULL)
pdf("figure3/corr_plots.pdf",width=10,height=4)
annotate_figure(fig, left = textGrob("Murphy et al. beta-coefficients", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Backhed et al. beta-coefficients", gp = gpar(cex = 1.3)))
dev.off()


