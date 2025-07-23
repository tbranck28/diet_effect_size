### Figure 3 and accompanying supplemental plots ###
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
library(ggh4x)
library(ggridges)
library(ggpubr)

#set up directory path
output_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")

######################################
### TABLE 2: Infant adonis results ### 
######################################
#making combined adonis results table from the two infant studies
backhed_adonis = as.data.frame(data.table::fread(paste0(output_dir,"/backhed_adonis_results.csv"), header=T, check.names = F))
backhed_adonis$Study = "Backhed infants"
murphy_adonis = as.data.frame(data.table::fread(paste0(output_dir,"/infants_murphy_adonis_results.csv"), header=T, check.names = F))
murphy_adonis$Study = "Murphy infants"
murphy_adonis$F = NULL
combined_adonis = rbind(backhed_adonis,murphy_adonis)
combined_adonis$FDR = p.adjust(combined_adonis$`Pr(>F)`,method="fdr")
combined_adonis$'Time point' = gsub("_"," ",combined_adonis$V1)
combined_adonis$'Time point' = gsub("meta","months",combined_adonis$'Time point')
write.csv(combined_adonis, paste0(output_dir,"/fig4_adonis_table.csv"))

#########################################
### FIGURE 3A: Infant maaslin results ### 
#########################################
#import all maaslin2 results for both infant datasets
backhed_baseline = as.data.frame(data.table::fread(paste0(output_dir,"/backhed_maaslin2_baseline_8.31/all_results.tsv"), header=T, check.names = F))
backhed_baseline$Study = 'B0'
backhed_four = as.data.frame(data.table::fread(paste0(output_dir,"/backhed_maaslin2_four_8.31/all_results.tsv"), header=T, check.names = F))
backhed_four$Study = 'B4'
backhed_twelve = as.data.frame(data.table::fread(paste0(output_dir,"/backhed_maaslin2_twelve_8.31/all_results.tsv"), header=T, check.names = F))
backhed_twelve$Study = 'B12'
murphy_baseline = as.data.frame(data.table::fread(paste0(output_dir,"/murphy_maaslin2_baseline_8.31/all_results.tsv"), header=T, check.names = F))
murphy_baseline$Study = 'M0'
murphy_three = as.data.frame(data.table::fread(paste0(output_dir,"/murphy_maaslin2_3_month_8.31/all_results.tsv"), header=T, check.names = F))
murphy_three$Study = 'M3'
murphy_twelve = as.data.frame(data.table::fread(paste0(output_dir,"/murphy_maaslin2_twelve_months_8.31/all_results.tsv"), header=T, check.names = F))
murphy_twelve$Study = 'M12'

maaslin_results_combined = rbind(backhed_baseline,backhed_four,backhed_twelve,murphy_baseline,murphy_three,murphy_twelve)
maaslin_results_combined$feature<-gsub("_", ".", maaslin_results_combined$feature)

keep = c("feeding","baseline_feeding","three_month_feeding","twelve_month_feeding")
maaslin_results_combined = maaslin_results_combined[maaslin_results_combined$metadata %in% keep, ]
maaslin_results_combined$value = gsub("no breast","no_breastfeeding",maaslin_results_combined$value)
maaslin_results_combined$value = str_replace(maaslin_results_combined$value, "mixed$", "mixed_BF_formula")
maaslin_results_combined$value = str_replace(maaslin_results_combined$value, "mixed_feeding$", "mixed_BF_formula")
maaslin_results_combined$value = gsub("bf","BF",maaslin_results_combined$value)
write.csv(maaslin_results_combined, paste0(output_dir,"/infants_maaslin_combined.csv"))

maaslin_results_combined = maaslin_results_combined %>% mutate(exp = case_when(Study == "B0" ~ "Backhed et al.",
                                                                               Study == "B4" ~ "Backhed et al.",
                                                                               Study == "B12" ~ "Backhed et al.",
                                                                               Study == "M0" ~ "Murphy et al.",
                                                                               Study == "M3" ~ "Murphy et al.",
                                                                               Study == "M12" ~ "Murphy et al.")) 

ggplot(maaslin_results_combined, aes(x=Study, y=coef, fill=value)) + 
  geom_boxplot() +
  facet_wrap(~exp, scale="free") +
  scale_y_continuous(limits=c(-6,6)) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

#string formatting for Collinsella
maaslin_results_combined$feature[maaslin_results_combined$feature=="X.Collinsella..massiliensis"] = "Collinsella.massiliensis"
maaslin_results_combined$Study = factor(maaslin_results_combined$Study, levels=c('B0','B4','B12','M0','M3','M12'))
maaslin_results_combined$feeding_time = paste0(maaslin_results_combined$Study,maaslin_results_combined$value)
anova = aov(abs(coef) ~ feeding_time,data=maaslin_results_combined)
summary(anova)

#values calculated/brought in from Figure2 script. Median, mean, and 95% CI values
#for HMP2 nonIBD RC1 beta coef's (the absolute value of the coefficients)
control_median = 0.18249
control_mean = 0.23162 
control_upper = 0.27333
control_lower = 0.18991

cairo_pdf(paste0(output_dir,"/infants_betas.pdf"), width=8, height = 6)
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

#####################################################################
### FIGURE 3B and Supplemental Figure 10: Infant maaslin heatmaps ### 
#####################################################################

maaslin_results_combined$feature = gsub("\\."," ",maaslin_results_combined$feature)

cairo_pdf(paste0(output_dir,"/infants_maaslinheatmap_supplement.pdf"), width=16, height = 16)
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

cairo_pdf(paste0(output_dir,"/infants_maaslinheatmap_subset.pdf"), width=10, height = 8)
ggplot(maaslin_results_combined_subset2, aes(value, reorder(feature,desc(feature)))) + 
  geom_tile(aes(fill = coef)) + 
  facet_nested(.~ exp + Study, scales = "free_x",space = "free_x") +
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

########################################################################################################################
### FIGURE 3C and Supplemental Figure 12: FOREST PLOT of effect size and significance (comparing infants and adults) ### 
########################################################################################################################
#subsetting maaslin results to keep only the significant tests
maaslin_results_combined_subset = maaslin_results_combined[maaslin_results_combined$qval < .25,] #remove this line if you don't care about keeping only significant associations
#maaslin_results_combined_subset = maaslin_results_combined

bug_list = unique(maaslin_results_combined_subset$feature)
maaslin_results_combined_subset = maaslin_results_combined[maaslin_results_combined$feature %in% bug_list,]
#maaslin_results_combined_subset = maaslin_results_combined

#import healthy adult human maaslin results and combining together
hmp2_maaslin2_sig_results_nonibd <- as.data.frame(data.table::fread(paste0(output_dir,"/hmp2_maaslin2_PCs_diseaseStrat_a_nonIBD/all_results.tsv"), header=T, check.names=T))
hmp2_maaslin2_sig_results_nonibd$Study = "HMP2_controls"
mlvs_maaslin2_sig_results <- as.data.frame(data.table::fread(paste0(output_dir,"/mlvs_maaslin2_PCs/all_results.tsv"), header=T, check.names=T))
mlvs_maaslin2_sig_results$Study = "MLVS"
adult_maas_results = rbind(hmp2_maaslin2_sig_results_nonibd,mlvs_maaslin2_sig_results)

#subsetting to keep only diet vs microbiome feature tests
keep_vars = c("RC1","RC2")
adult_maas_results = adult_maas_results[adult_maas_results$metadata %in% keep_vars,]

adult_maas_results$feature = gsub("X.","",adult_maas_results$feature)
adult_maas_results$feature = gsub("Collinsella..massiliensis","Collinsella.massiliensis",adult_maas_results$feature)

#get list of significant bugs from the infant maaslin tests and subset the healthy adult human maaslin results
infant_sig_bugs = unique(maaslin_results_combined_subset$feature)
infant_sig_bugs = gsub(" ",".",infant_sig_bugs)
adult_maas_results = adult_maas_results[adult_maas_results$feature %in% infant_sig_bugs,]
adult_maas_results$genus = gsub("\\..*","",adult_maas_results$feature)

#import hadza maaslin results
hadza_seas_maaslin2_sig_results = as.data.frame(data.table::fread(paste0(output_dir,"/hadza_seasonal_maaslin2_output_binary/all_results.tsv"), header=T, check.names=T))
hadza_seas_maaslin2_sig_results$Study = "Hadza"

#keeping only diet tests
hadza_seas_maaslin2_sig_results = hadza_seas_maaslin2_sig_results[hadza_seas_maaslin2_sig_results$metadata == "Binary_season",]

#dataframe formatting
hadza_seas_maaslin2_sig_results = hadza_seas_maaslin2_sig_results[!grepl("f__", hadza_seas_maaslin2_sig_results$feature),]
hadza_seas_maaslin2_sig_results$genus = gsub(".*\\__","",hadza_seas_maaslin2_sig_results$feature)
hadza_seas_maaslin2_sig_results$genus = gsub("_.*","",hadza_seas_maaslin2_sig_results$genus)
maaslin_results_combined_subset$genus = gsub(" .*","",maaslin_results_combined_subset$feature)
adult_maas_results$feature = gsub("\\."," ",adult_maas_results$feature)

#merge the significant bug-subsetted healthy adult maaslin results with the subsetted infant maaslin results
test_merged = merge(adult_maas_results, maaslin_results_combined_subset, by="feature",all = F)
names(test_merged)[23] = "genus"

#adding the hadza results
test_merged2 = merge(test_merged, hadza_seas_maaslin2_sig_results, by="genus",all = T)

#keeping only RC1 for simplicity##
test_merged2 = test_merged2[!test_merged2$metadata.x == "RC2",]
test_merged2 = test_merged2[!is.na(test_merged2$feature.x),]

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

cairo_pdf(paste0(output_dir,"/effect_size_FOREST_subset.pdf"),width=8,height=6)
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
cairo_pdf(paste0(output_dir,"/effect_size_FOREST_Legend.pdf"),width=3,height=6)
as_ggplot(ES_legend)
dev.off()

cairo_pdf(paste0(output_dir,"/effect_size_FOREST_subset_nolegend.pdf"),width=6,height=4)
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

cairo_pdf(paste0(output_dir,"/qval_FOREST.pdf"),width=7,height=5)
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

#########################################################################
### Supplemental Figure 11: Spearman Murphy vs. Backhed maaslin betas ### 
#########################################################################
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
pdf(paste0(output_dir,"/corr_plots.pdf"),width=10,height=4)
annotate_figure(fig, left = textGrob("Murphy et al. beta-coefficients", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Backhed et al. beta-coefficients", gp = gpar(cex = 1.3)))
dev.off()
