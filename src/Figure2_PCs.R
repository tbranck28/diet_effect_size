setwd("/Users/tobynbranck/Documents/Effect_size/analysis/")
library(ggpubr)
library(tidyverse)
library(gridExtra)
library(patchwork)
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
require(cowplot)

###### Fig 2A: MANTEL PLOT################################################################################
mantel_result_hmp2 <- data.table::fread("HMP2_mantel_results.csv", header=T, check.names=F)
mantel_result_mlvs <- data.table::fread("MLVS_mantel_results.csv", header=T, check.names=F)

mantel_result_mlvs$WeekNumber = NULL
combined = rbind(mantel_result_hmp2, mantel_result_mlvs)
combined[,c("group")][is.na(combined[,c("group")])] <- "nonIBD"
combined$Label = c("HMP2", "HMP2", "HMP2", "HMP2", "HMP2", "HMP2", "HMP2", "HMP2", "HMP2", "HMP2", "HMP2", "HMP2", "HMP2", "HMP2", "HMP2", "HMP2", "HMP2", "HMP2", "MLVS_Week", "MLVS_Week")
combined$Label = paste(combined$Label, combined$Bin)

pal <- c(CD = "indianred3", UC = "lemonchiffon3", nonIBD = "azure3")
combined$fdr = p.adjust(combined$p,method="fdr")
ann <- combined %>%
  filter(p < 0.05) %>%
  mutate(signif = case_when( fdr <= 0.1 ~ "*")) %>%
  merge(., combined, by=colnames(combined), all=TRUE)
bar <- ggplot(combined, aes(x=factor(Label, levels=c("HMP2 (0_7)", "HMP2 (8_15)", "HMP2 (16_24)", "HMP2 (25_33)","HMP2 (34_42)","HMP2 (43_57)","MLVS_Week 1", "MLVS_Week 2")), y=r, fill=group)) +
  geom_bar(stat="identity", position="dodge",color="black") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1),legend.position="top") +
  labs(x="Time (Week Bins)", y='Mantel r', alpha="p-value", colour="Group", fill="Diagnosis") +
  scale_colour_manual(values=pal, labels=c("CD"="CD", "nonIBD"="nonIBD", "UC"="UC")) +
  scale_fill_manual(values=pal, labels=c("CD"="CD", "nonIBD"="nonIBD", "UC"="UC")) +
  geom_text(data=ann, aes(x=Label, y=r, label=signif), position=position_dodge(width=0.9), size=10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=14)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.title=element_text(size=14))
cairo_pdf(file = "figure2/Combined_plot_HMP2_MLVS_stratified.pdf", width=7, height = 5)
bar
dev.off()

############# Supplemental adonis plot (individual components) ###########################################################################
hmp2_CD_adonis <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2CD_adonis_results.csv", header=T, check.names=T))
hmp2_UC_adonis <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2UC_adonis_results.csv", sep=",", header=T, check.names=F))
hmp2_nonIBD_adonis <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2nonIBD_adonis_results.csv", sep=",", header=T, check.names=F))
mlvs_adonis <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_adonis_results.csv", sep=",", header=T, check.names=F))
hmp2_CD_adonis$'P-Value'= hmp2_CD_adonis$P.Value
hmp2_CD_adonis$P.Value = NULL
hmp2_CD_adonis$Study = 'HMP2_CD'
hmp2_UC_adonis$Study = 'HMP2_UC'
hmp2_nonIBD_adonis$Study = 'HMP2_nonIBD'
mlvs_adonis$Study = 'MLVS'
df3 <- rbind(hmp2_CD_adonis, hmp2_UC_adonis, hmp2_nonIBD_adonis, mlvs_adonis)
#had<-data.frame(0.093,NA,"*","Season-dependent traditional diet", 0.001, "Hadza")
#names(had)<-c("R2","p_adj","stars","diet_variable","P-Value","Study")
#newdf <- rbind(df3, had)
newdf = df3
newdf$p_adj = NA
newdf$p_adj = p.adjust(newdf$`P-Value`, method = "fdr")
newdf$stars = cut(newdf$p_adj, c(0,0.1), labels = c("*"))
newdf_truncate = newdf
#adonis_bar = ggplot(data=newdf, aes(x=reorder(diet_variable,R2),y=R2,fill=factor(Study, levels = c("Hadza","HMP2_CD","HMP2_UC", "HMP2_nonIBD","MLVS")))) +
adonis_bar = ggplot(data=newdf, aes(x=reorder(diet_variable,R2),y=R2,fill=factor(Study, levels = c("HMP2_CD","HMP2_UC", "HMP2_nonIBD","MLVS")))) +
  geom_bar(position=position_dodge(.9),stat="identity") + 
  coord_flip() +
  labs(fill = "Study") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("") +
  scale_y_continuous(expand = c(0,0.005)) +
  ylab(expression(R^'2')) +
  scale_fill_manual(values=c("indianred3","lemonchiffon3","azure3","azure4")) +
  theme(legend.position = c(.8, .2)) +
  geom_text(data = newdf,x=newdf$diet_variable,label = newdf$stars, y = newdf$R2 + .002, vjust=.9)

cairo_pdf(file = "figure2/adonis_plot_IndividualComponents_SUPP.pdf", width=8, height = 7)
adonis_bar
dev.off()
#################################################################################################################

######### Fig.2B: adonis from PCs #########
mlvs_adonis <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_adonisPCs_results.csv", sep=",", header=T, check.names=F))
hmp2_adonis_pcs <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/HMP2_adonisPCs_byDisease_results.csv", header=T, check.names=T))
hmp2_adonis_pcs$adonis_label = hmp2_adonis_pcs$cat
hmp2_adonis = hmp2_adonis_pcs[,1:5]
mlvs_adonis$adonis_label <- paste("MLVS", mlvs_adonis$adonis_label, sep=" ")
hmp2_adonis$adonis_label <- paste("HMP2", hmp2_adonis$adonis_label, sep=" ")
 
df3 <- rbind(hmp2_adonis, mlvs_adonis)
had<-data.frame("Hadza",0.093,0.001, NA,"*")
names(had)<-names(df3)
newdf <- rbind(df3, had)
newdf$study = c("HMP2","HMP2","HMP2","HMP2","HMP2","HMP2","MLVS","MLVS","Hadza")
 
selection = as.vector(unique(newdf$adonis_label[1:6]))
#newdf[nrow(newdf) + 1,] = c("Backhed_infant", "0.03497","0.046", "NA", "*", "Backhed")

newdf$adonis_rsq = as.numeric(newdf$adonis_rsq)
newdf$label2 = gsub(" RC.*", "",newdf$adonis_label)
newdf$label2[newdf$label2=='HMP2 nonIBD'] = "nonIBD"
newdf$label2[newdf$label2=='MLVS'] = "nonIBD"
newdf$p_adj = p.adjust(newdf$adonis_pval, method = "fdr")
newdf$stars = cut(newdf$p_adj, c(0,0.1), labels = c("*"))
scaleFUN <- function(x) sprintf("%.2f", x)

adonis_bar_PCs = ggplot(data=newdf, aes(x=reorder(adonis_label,-adonis_rsq),y=adonis_rsq,fill=factor(label2, levels = c("Hadza","HMP2 CD","HMP2 UC", "nonIBD")))) +
  geom_bar(position=position_dodge(.9),stat="identity", color="black") + 
  labs(fill = "Study") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Diet") +
  ylab(expression(R^'2')) +
  scale_fill_manual(values=c("darkseagreen","indianred3","lemonchiffon3","azure3")) +
  theme(legend.position = c(.25, .75),legend.text=element_text(size=12)) +
  scale_y_continuous(limits = c(0,.11),labels = scaleFUN) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = newdf,x=newdf$adonis_label,label = newdf$stars, y = as.numeric(newdf$adonis_rsq) + .004, vjust=.8,size=10) +
  theme(legend.key.size = unit(0.4,"cm"),legend.position="top") +
  theme(axis.text.x = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=14)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.title=element_text(size=14)) +
  theme(legend.position = "none")
cairo_pdf(file = "figure2/adonis_plot_noInset.pdf", width=6, height = 6)
adonis_bar_PCs
dev.off()

newdf_west = newdf[1:8,]
newdf_west$adonis_rsq = as.numeric(newdf_west$adonis_rsq)
newdf_west$adonis_pval = as.numeric(newdf_west$adonis_pval)
adonis_western = ggplot(data=newdf_west, aes(x=reorder(adonis_label,-adonis_rsq),y=adonis_rsq,fill=factor(label2, levels = c("HMP2 CD","HMP2 UC", "HMP2 nonIBD","MLVS")))) +
  geom_bar(position=position_dodge(.9),stat="identity",color="black") + 
  labs(fill = "Study") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Diet") +
  ylab(expression(R^'2')) +
  scale_fill_manual(values=c("indianred3","lemonchiffon3","azure3","azure4")) +
  theme(legend.position = "none") +
  ylim(0,.008) +
  expand_limits(y=c(0,.008)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.title.x=element_blank()) +
  theme(axis.text.y = element_text(size=8),axis.title.y=element_text(size=8)) +
  geom_text(data = newdf_west,x=newdf_west$adonis_label,label = newdf_west$stars, y = newdf_west$adonis_rsq + .0002, vjust=.65,size=10) +
  theme(plot.background = element_rect(colour = "black", fill=NA, size=.5)) +
  theme(axis.text.x = element_text(size = 10)) + 
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.title.y = element_text(size = 10)) +
  theme(legend.position="none")

figure = ggarrange(adonis_bar + rremove("xlab"),adonis_western + rremove("ylab") + rremove("xlab"), ncol=2,nrow=1)
adonis_bar_2plots = annotate_figure(figure, bottom = textGrob("Diet"))
tss = adonis_bar_PCs + inset_element(adonis_western, .5,.5,1,1, align_to = 'plot')

cairo_pdf(file = "figure2/adonis_plot_inset.pdf", width=6, height = 6)
tss
dev.off()

###### Code to truncate number of individual variables in adonis plot #############################################################
#subset so that can truncate plot
truncdf = newdf_truncate
 truncdf = truncdf[order(truncdf$Study,-truncdf$R2),]
 Names_list = list("Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)",
                   "Tea or coffee no sugar and no sugar replacement",
                   "Probiotic",
                   "Galactose_fo",
                   "VitD_calciferol_fo",
                   "Iron_from_fs",
                   "Total_PUFA_fs",
                   "Sucrose_from_fs",
                   "Total_MUFA_fs",
                   "Total_Trans_fo",
                   "Total_SFA_fs",
                   "Total_Folate_fo",
                   "Choline_fs",
                   "VitD_calciferol_fs",
                   "Alcohol_fo",
                   "VitK_phylloquinone_fo",
                   "Energy_fs",
                   "Shellfish (shrimp, lobster, scallops, etc.)",
                   "Water",
                   "Processed meat (other red or white meat such as lunch meat, ham, salami, bologna")

 truncdf = truncdf[ ! truncdf$diet_variable %in% Names_list, ]
 truncdf = as.data.frame(truncdf)

 truncdf$diet_variable[truncdf$diet_variable=="Insol_Dietary_Fiber_fs"] <- "Insoluble Dietary Fiber"
 truncdf$diet_variable[truncdf$diet_variable=="Calcium_fs"] <- "Calcium"
 truncdf$diet_variable[truncdf$diet_variable=="VitB12_cobalamin_fs"] <- "Vitamin B12"
 truncdf$diet_variable[truncdf$diet_variable=="Sol_Dietary_Fiber_fo"] <- "Soluble Dietary Fiber"
 truncdf$diet_variable[truncdf$diet_variable=="Glucose_fs"] <- "Glucose"
 truncdf$diet_variable[truncdf$diet_variable=="Total_Carb_from_fs"] <- "Carbohydrates"
 truncdf$diet_variable[truncdf$diet_variable=="Starch_fo"] <- "Starch"
 truncdf$diet_variable[truncdf$diet_variable=="Total_Protein_fs"] <- "Protein"
 truncdf$diet_variable[truncdf$diet_variable=="Red meat (beef, hamburger, pork, lamb)"] <- "Red Meat"
 truncdf$diet_variable[truncdf$diet_variable=="White meat (chicken, turkey, etc.)"] <- "White meat"
 truncdf$diet_variable[truncdf$diet_variable=="Vegetables (salad, tomatoes, onions, greens, carrots, peppers, green beans, etc)"] <- "Vegetables"
 truncdf$diet_variable[truncdf$diet_variable=="Starch (white rice, bread, pizza, potatoes, yams, cereals, pancakes, etc.)"] <- "Starch (bread, rice, etc.)"
 truncdf$diet_variable[truncdf$diet_variable=="Sweets (pies, jam, chocolate, cake, cookies, etc.)"] <- "Sweets"
 truncdf$diet_variable[truncdf$diet_variable=="Beans (tofu, soy, soy burgers, lentils, Mexican beans, lima beans etc)"] <- "Beans"
 truncdf$diet_variable[truncdf$diet_variable=="Fruits (no juice) (Apples, raisins, bananas, oranges, strawberries, blueberries"] <- "Fruits (no juice)"
 truncdf$diet_variable[truncdf$diet_variable=="Dairy (milk, cream, ice cream, cheese, cream cheese)"] <- "Dairy"
 truncdf$diet_variable[truncdf$diet_variable=="Soft drinks, tea or coffee with sugar (corn syrup, maple syrup, cane sugar, etc)"] <- "Soft drinks, tea or coffee with sugar"
 truncdf$diet_variable[truncdf$diet_variable=="Yogurt or other foods containing active bacterial cultures (kefir, sauerkraut)"] <- "Yogurt"
 truncdf$diet_variable[truncdf$diet_variable=="Whole grains (wheat, oats, brown rice, rye, quinoa, wheat bread, wheat pasta)"] <- "Whole grains"
 truncdf$diet_variable[truncdf$diet_variable=="Fruit juice (orange, apple, cranberry, prune etc.)"] <- "Fruit juice"
 truncdf$diet_variable[truncdf$diet_variable=="Fish (fish nuggets, breaded fish, fish cakes, salmon, tuna, etc.)"] <- "Fish"
 truncdf$diet_variable[truncdf$diet_variable=="Diet soft drinks, tea or coffee with sugar (Stevia, Equal, Splenda etc)"] <- "Diet soft drinks"

 adonis_bar_truncated = ggplot(data=truncdf, aes(x=reorder(diet_variable,R2),y=R2,fill=factor(truncdf$Study, levels = c("Hadza","HMP2_CD","HMP2_UC","HMP2_nonIBD","MLVS")))) +
   geom_bar(position=position_dodge(.9),stat="identity") +
   coord_flip() +
   labs(fill = "Study") +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   xlab("Diet feature") +
  scale_y_continuous(expand = c(0,0.005)) +
   ylab(expression(R^'2')) +
   scale_fill_manual(values=c('darkseagreen',"indianred3","lemonchiffon3",'lightcyan3','#56B4E9')) +
   theme(legend.position = c(.8, .2)) +
   geom_text(data = truncdf,x=truncdf$diet_variable,label = truncdf$stars, y = truncdf$R2 + .002, vjust=.9) +
   theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14)) +
   theme(legend.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))

ggsave(plot=adonis_bar_truncated, file="adonis_plot_NaturePres.png", dpi=300, height=8, width=10)

########################### Supplemental figure: maaslin2 heatmap ##################################################################################
library(reshape2)
hmp2_maaslin2_sig_results_cd <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_diseaseStrat_CD/all_results.tsv", header=T, check.names=T))
hmp2_maaslin2_sig_results_cd$metadata <- paste("HMP2_cd", hmp2_maaslin2_sig_results_cd$metadata, sep="_")

hmp2_maaslin2_sig_results_uc <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_diseaseStrat_UC/all_results.tsv", header=T, check.names=T))
hmp2_maaslin2_sig_results_uc$metadata <- paste("HMP2_uc", hmp2_maaslin2_sig_results_uc$metadata, sep="_")

hmp2_maaslin2_sig_results_nonibd <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_diseaseStrat_a_nonIBD/all_results.tsv", header=T, check.names=T))
hmp2_maaslin2_sig_results_nonibd$metadata <- paste("HMP2_nonIBD", hmp2_maaslin2_sig_results_nonibd$metadata, sep="_")

mlvs_maaslin2_sig_results <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_maaslin2_PCs/all_results.tsv", header=T, check.names=T))
mlvs_maaslin2_sig_results$metadata <- paste("MLVS", mlvs_maaslin2_sig_results$metadata, sep="_")

hadza_seas_maaslin2_sig_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hadza_seasonal_maaslin2_output_binary/all_results.tsv", header=T, check.names=T))
hadza_seas_maaslin2_sig_results$metadata[hadza_seas_maaslin2_sig_results$metadata == 'Season'] = 'Hadza_Season'
hadza_seas_maaslin2_sig_results$metadata[hadza_seas_maaslin2_sig_results$metadata == 'Host_Age'] = 'Hadza_host_age'

nhp_sig_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/dada2_test/nonhumanPrimates_maaslin2_output/all_results.tsv", header=T, check.names=T))
nhp_sig_results$metadata[nhp_sig_results$metadata == 'Host_diet'] = 'NHP_Host_diet'

mouse_sig_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mouse_complexdiet_maaslin2_binned/all_results.tsv", header=T, check.names=T))
mouse_sig_results$metadata <- "Mouse"

testt = rbind(hmp2_maaslin2_sig_results_uc,hmp2_maaslin2_sig_results_cd,hmp2_maaslin2_sig_results_nonibd,mlvs_maaslin2_sig_results)
testt$metadata = gsub("diagnosis2", "dysbiosis", testt$metadata)
filttest = apply(testt, 1, function(r) any(r %in% c("HMP2_cd_antibiotic","HMP2_uc_antibiotic","HMP2_nonIBD_antibiotic", "HMP2_cd_consent_age","HMP2_uc_consent_age", "HMP2_nonIBD_consent_age", "HMP2_cd_dysbiosis", "HMP2_uc_dysbiosis","HMP2_nonIBD_dysbiosis", "MLVS_antibiotic", "MLVS_consent_age")))
test_covariates = testt[filttest,]

test_microbiome = apply(testt, 1, function(r) any(r %in% c("MLVS_RC1","MLVS_RC2","HMP2_cd_RC1",
                                                           "HMP2_cd_RC2","HMP2_uc_RC1","HMP2_uc_RC2","HMP2_nonIBD_RC1","HMP2_nonIBD_RC2")))
test_microbiome = testt[test_microbiome,]
testt$metadata[testt$metadata == 'HMP2_diagnosis2'] <- 'HMP2_dysbiosis'
testt = testt[!(testt$metadata == "MLVS_WeekNum" | testt$metadata == "MLVS_SampNum"),]

test_covariates$metadata[test_covariates$metadata == 'HMP2_diagnosis2'] <- 'HMP2_dysbiosis'

test_covariates$type = "Covariates"
test_microbiome$type = "Taxa"
heatmap_df = rbind(test_covariates,test_microbiome)
heatmap_df$feature = gsub("X.","",heatmap_df$feature)
heatmap_df$feature = gsub("\\."," ",heatmap_df$feature)
sorted = sort(unique(heatmap_df$feature))
heatmap_df$metadata = gsub("cd", "CD", heatmap_df$metadata)
heatmap_df$metadata = gsub("uc", "UC", heatmap_df$metadata)
heatmap_df$metadata = gsub("RC","DP",heatmap_df$metadata)

heatmaps = ggplot(heatmap_df, aes(factor(feature,levels = rev(sorted)), factor(metadata,levels = c("HMP2_CD_consent_age","HMP2_UC_consent_age","HMP2_nonIBD_consent_age","HMP2_CD_antibiotic","HMP2_UC_antibiotic","HMP2_nonIBD_antibiotic","HMP2_CD_dysbiosis","HMP2_UC_dysbiosis","HMP2_nonIBD_dysbiosis","HMP2_CD_DP1",
                                                           "HMP2_CD_DP2","HMP2_UC_DP1","HMP2_UC_DP2","HMP2_nonIBD_DP1","HMP2_nonIBD_DP2","MLVS_antibiotic","MLVS_consent_age","MLVS_DP2","MLVS_DP1")), fill= coef)) + 
  geom_tile(color="black") +
  coord_equal() + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10)) +
  theme(axis.text.y = element_text(size=9),axis.title.y = element_text(size = 10)) +
  scale_fill_gradient2(low="royalblue4", mid="snow", high="red4", midpoint=0, limits = c(-6,6)) +
  geom_text(aes(label = ifelse(qval <= .25,"*","")), vjust = 0.77) + 
  theme(axis.line.y.right = element_line(color = "red")) +
  xlab("") +
  theme(axis.title=element_text(size=14)) +
  labs(fill = expression(paste(beta,' coefficient'))) +
  theme(axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 12, b = 0, l = 0))) +
  #theme(plot.margin=unit(c(-3.4,3,1,1), "cm")) +
  theme(legend.key.size = unit(.35, 'cm')) +
  theme(legend.position = c(-.12, -.05)) +
  theme(axis.title.x=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line =element_blank()) +
  theme(panel.background = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank()) +
  facet_grid(~type,scales="free") +
  theme(strip.text = element_text(size = 12)) +
  theme(text=element_text(family = "Arial")) +
  theme(axis.text.y = element_text(face = "italic"))

cairo_pdf("figure2/maaslin2_heatmaps_fig2SUPP.pdf", width=14, height = 15)
heatmaps
dev.off()

#### Fig. 2C ### BOXPLOT of maaslin2 results
for_boxplot = rbind(testt, hadza_seas_maaslin2_sig_results)
for_boxplot = for_boxplot[- grep("antibiotic", for_boxplot$metadata),]
for_boxplot = for_boxplot[- grep("dysbiosis", for_boxplot$metadata),]
for_boxplot = for_boxplot[- grep("age", for_boxplot$metadata),]
for_boxplot = for_boxplot[- grep("WeekBin", for_boxplot$metadata),]

for_boxplot %>% group_by(reorder(metadata,-abs(coef))) %>% summarise(mean = mean(abs(coef)), sd=sd(abs(coef)))

for_boxplot$metadata <- gsub('_uc_', ' UC ', for_boxplot$metadata)
for_boxplot$metadata <- gsub('_cd_', ' CD ', for_boxplot$metadata)
for_boxplot$metadata <- gsub('_nonIBD_', ' nonIBD ', for_boxplot$metadata)
for_boxplot$metadata <- gsub('MLVS_', 'MLVS ', for_boxplot$metadata)
for_boxplot$metadata <- gsub('Hadza_Season', 'Hadza (seasons) ', for_boxplot$metadata)
for_boxplot$metadata[for_boxplot$metadata == 'Binary_season'] = "Hadza Season"
for_boxplot$metadata = gsub('RC1','Dietary pattern 1',for_boxplot$metadata)
for_boxplot$metadata = gsub('RC2','Dietary pattern 2',for_boxplot$metadata)

write.csv(for_boxplot,"beta_boxplot_df.csv")

beta_coeff_box = ggplot(for_boxplot,aes(x=reorder(metadata,-abs(coef)), y=abs(coef),fill=metadata)) +
  geom_boxplot() +
  labs(y=expression(beta~'coefficients'), x="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) + 
  scale_fill_manual(values=c("darkseagreen","indianred3","indianred3","azure3","azure3","lemonchiffon3","lemonchiffon3","azure4","azure4")) +
  theme(legend.position=c(.85,.6)) +
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.title.y = element_text(size = 10)) +
  theme(legend.title=element_blank()) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.key=element_blank()) + 
  theme(legend.position = 'none')

cairo_pdf("figure2/beta_coef_box_fig2.pdf", width=5, height = 5)
beta_coeff_box
dev.off()

####calculating mean, median, CI for HMP2 nonIBD RC1 to add a line and CI
#### to the infant betas boxplots in Figure 3
hmp2_nonIBD_RC1 = for_boxplot[for_boxplot$metadata == 'HMP2 nonIBD Dietary pattern 1',]
hmp2_nonIBD_RC1$coef = abs(hmp2_nonIBD_RC1$coef)
l.model = lm(coef~1,hmp2_nonIBD_RC1)
confint(l.model,level=0.95)
mean(hmp2_nonIBD_RC1$coef)
median(hmp2_nonIBD_RC1$coef)
#### Figure 2D,E,F: specific species interaction with diet ##############
hmp2meta_df <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/HMP2_PCA_df.csv", header=T, check.names=T))
rownames(hmp2meta_df) = hmp2meta_df$V1
hmp2meta_df$V1 = NULL
hmp2meta_df$diagnosis <- gsub('a_nonIBD', 'nonIBD', hmp2meta_df$diagnosis)

scaleFUN <- function(x) sprintf("%.2f", x)

LOG <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(log2(y))
}

zeros_df = as.data.frame(aggregate(.~diagnosis,hmp2meta_df,FUN=function(x) sum(x==0)))
totals = as.data.frame(t(table(hmp2meta_df$diagnosis)))
zeros_df$totals = totals$Freq
zeros_df[,2:16] = zeros_df[,2:16]/zeros_df$totals

a_had_bar = ggplot(zeros_df, aes(x=diagnosis, y=zeros_df$Anaerostipes.hadrus, fill=diagnosis)) +
  geom_bar(stat="identity") +
  ylab("Zeros (%)") + 
  scale_fill_manual("diagnosis", 
                    breaks = c("nonIBD", "CD", "UC"),
                    values = c("azure4", "indianred3","burlywood3")) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=10)) +
  theme(legend.title=element_text(size=10)) +
  xlab("")
  

ahad_greater_zero = hmp2meta_df[hmp2meta_df$Anaerostipes.hadrus != 0, ]

smooth_CD = loess(Anaerostipes.hadrus ~ RC2, data=ahad_greater_zero[ahad_greater_zero$diagnosis=="CD",], span=0.10)
pred_CD = predict(smooth_CD)

scatter_2 = ggplot(hmp2meta_df, aes(x=RC2, y = LOG(Anaerostipes.hadrus), color = diagnosis)) +
  geom_point(alpha=0.35, size=3) +
  geom_smooth(aes(col=diagnosis), method="loess",se=TRUE, fullrange=TRUE,size=1.5,span=.9) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_colour_manual("Diagnosis", 
                      breaks = c("nonIBD", "CD", "UC"),
                      values = c("azure4", "indianred3","burlywood3")) +
  xlab("HMP2 dietary pattern RC2") +
  ylab("Anaerostipes hadrus") +
  theme(legend.position=c(1,.875)) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(axis.text.x = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=10)) +
  theme(legend.key=element_blank()) +
  theme(legend.title=element_text(size=10)) +
  theme(legend.position = 'none')

cairo_pdf("figure2/Ahadrus_smoothregression.pdf",width=5,height=4)
scatter_2
dev.off()

X<-split(hmp2meta_df, hmp2meta_df$diagnosis)
CD_df = as.data.frame(X[1])
UC_df = as.data.frame(X[3])
nonIBD_df = as.data.frame(X[2])

#bring in Hadza bugs and metadata
hadza_meta_df <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hadza_seasonal_metadata_processed.csv", header=T, check.names=T))
hadza_bug_df <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hadza_seasonal_species.csv", header=T, check.names=T))
had_df = cbind(hadza_bug_df,hadza_meta_df)
rownames(had_df) = had_df$V1
had_df$V1 = NULL

had_df$Season[had_df$Season == "2013-LD"] <- '2013 Late Dry'
had_df$Season[had_df$Season == "2014-ED"] <- '2014 Early Dry'
had_df$Season[had_df$Season == "2014-EW"] <- '2014 Early Wet'
had_df$Season[had_df$Season == "2014-LD"] <- '2014 Late Dry'
had_df$Season[had_df$Season == "2014-LW"] <- '2014 Late Wet'

hadza_box = ggplot(had_df, aes(x = Binary_season, y = had_df$g__Anaerostipes)) + 
  geom_boxplot() +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  ylab("g__Anaerostipes") +
  #theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(axis.text.x = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_text(size = 14)) 
  #annotate("text", x=.5, y=0.035, label= "\U03B2 coefficient: 4.69\nFDR: 3.30e-05",hjust = 0,size=5)

cairo_pdf("figure2/hadza_Ahadrus_boxplot.pdf", width=5, height = 4)
hadza_box
dev.off()

### MLVS scatter ####
MLVSmeta_df <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_combined_df.csv", header=T, check.names=T))
rownames(MLVSmeta_df) = MLVSmeta_df$V1
MLVSmeta_df$V1 = NULL
mlvs_scatter_RC1 = ggplot(MLVSmeta_df, aes(x=MLVSmeta_df$RC1, y = LOG(MLVSmeta_df$Anaerostipes.hadrus))) +
  geom_point(alpha=0.4, size=3,color="azure4") +
  geom_smooth(method="loess",se=TRUE, fullrange=TRUE,size=1.5,span=.9,color="darkgrey") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  xlab("MLVS dietary pattern RC1") +
  ylab("Anaerostipes.hadrus") +
  theme(legend.position=c(.9,.75)) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(axis.text.x = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12))
  #annotate("text", x=-3.5, y=0.15, label= "\U03B2 coefficient: 0.23\nFDR: 0.35",hjust = 0,size=5)


cairo_pdf("figure2/mlvs_Ahadrus_RC1_scatter.pdf", width=5, height = 4)
mlvs_scatter_RC1
dev.off()

vTert = quantile(MLVSmeta_df$RC1, c(0:3/3))

# classify values
MLVSmeta_df$tert = with(MLVSmeta_df, 
               cut(RC1, 
                   vTert, 
                   include.lowest = T, 
                   labels = c("Low", "Medium", "High")))

mlvs_scatter_RC1_ER = ggplot(MLVSmeta_df, aes(x=MLVSmeta_df$tert, y = MLVSmeta_df$Eubacterium.rectale)) +
  geom_boxplot(outlier.shape = NA,lwd=1) +
  geom_jitter(shape=16, position=position_jitter(0.3),alpha=.3,size=3) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  xlab("MLVS dietary pattern RC1") +
  ylab("Eubacterium rectale") + 
  theme(legend.position=c(.9,.75)) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(axis.text.x = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14,face="italic")) +
  theme(legend.text=element_text(size=12))

mlvs_scatter_RC1_ER_log = ggplot(MLVSmeta_df, aes(x=MLVSmeta_df$tert, y = LOG(MLVSmeta_df$Eubacterium.rectale))) +
  geom_boxplot(outlier.shape = NA,lwd=1) +
  geom_jitter(shape=16, position=position_jitter(0.3),alpha=.3,size=3) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  xlab("MLVS dietary pattern RC1") +
  ylab("Eubacterium rectale") + 
  theme(legend.position=c(.9,.75)) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(axis.text.x = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14,face="italic")) +
  theme(legend.text=element_text(size=12))

cairo_pdf("figure2/mlvs_EubactRect_RC1_scatter.pdf", width=5, height = 4)
mlvs_scatter_RC1_ER
dev.off()

cairo_pdf("figure2/mlvs_EubactRect_RC1_scatter_LOG.pdf", width=5, height = 4)
mlvs_scatter_RC1_ER_log
dev.off()

mlvs_scatter_RC2 = ggplot(MLVSmeta_df, aes(x=MLVSmeta_df$RC2, y = MLVSmeta_df$Anaerostipes.hadrus)) +
  geom_point(alpha=0.4, size=3,color="azure4") +
  geom_smooth(method = "lm", se=F, fullrange=TRUE, color="black") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  xlab("MLVS dietary pattern RC2") +
  ylab("Anaerostipes.hadrus") +
  theme(legend.position=c(.95,.75)) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(axis.text.x = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12))
  #annotate("text", x=-3.5, y=0.15, label= "\U03B2 coefficient: 0.20\nFDR: 0.45",hjust = 0,size=5)

cairo_pdf("figure2/mlvs_Ahadrus_RC2_scatter.pdf", width=5, height = 4)
mlvs_scatter_RC2
dev.off()

###########################################################
#Supplemental: PCs vs. raw mlvs and hmp2 maaslin#
###########################################################

# boxplot for hmp2 maaslin effect sizes - factor scores
hmp2_factor_results_CD = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_diseaseStrat_CD/all_results.tsv", header=T, check.names=T))
hmp2_factor_results_CD$metadata[hmp2_factor_results_CD$metadata == 'diagnosis2'] <- 'dysbiosis'

hmp2_factor_results_UC = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_diseaseStrat_UC/all_results.tsv", header=T, check.names=T))
hmp2_factor_results_UC$metadata[hmp2_factor_results_UC$metadata == 'diagnosis2'] <- 'dysbiosis'

hmp2_factor_results_nonIBD = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_PCs_diseaseStrat_a_nonIBD/all_results.tsv", header=T, check.names=T))
hmp2_factor_results_nonIBD$metadata[hmp2_factor_results_nonIBD$metadata == 'diagnosis2'] <- 'dysbiosis'

hmp2_factor_results = rbind(hmp2_factor_results_CD,hmp2_factor_results_UC,hmp2_factor_results_nonIBD)

filttest = apply(hmp2_factor_results, 1, function(r) any(r %in% c("antibiotic", "consent_age", "diagnosis", "dysbiosis")))
test_covariates = hmp2_factor_results[filttest,]
test_microbiome = hmp2_factor_results[!filttest,]

hmp2_factor_results$metadata[hmp2_factor_results$metadata == "RC1"] = "dietary pattern 1"
hmp2_factor_results$metadata[hmp2_factor_results$metadata == "RC2"] = "dietary pattern 2"
hmp2_factor_results$metadata[hmp2_factor_results$metadata == "consent_age"] = "consent age"

hmp2_factor = ggplot(hmp2_factor_results,aes(x=reorder(metadata,-abs(coef)), y=abs(coef))) +
  geom_boxplot() +
  labs(y=expression(beta~'coefficients'), x="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.text.x = element_text(size = 9)) +
  theme(axis.text.y = element_text(size = 12)) + 
  labs(title = "HMP2", subtitle = "Maaslin2 Effect Sizes - PC Scores")

#boxplot for hmp2 maaslin effect sizes - raw diet
hmp2_raw_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/hmp2_maaslin2_raw/all_results.tsv", header=T, check.names=T))
hmp2_raw_results$metadata[hmp2_raw_results$metadata == 'diagnosis2'] <- 'dysbiosis'

filttest = apply(hmp2_raw_results, 1, function(r) any(r %in% c("antibiotic", "consent_age", "diagnosis", "diagnosis2")))
test_covariates = hmp2_raw_results[filttest,]
test_microbiome = hmp2_raw_results[!filttest,]

hmp2_raw = ggplot(hmp2_raw_results,aes(x=reorder(metadata,-abs(coef)), y=abs(coef))) +
  geom_boxplot() +
  labs(y=expression(beta~'coefficients'), x="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.text.x = element_text(size = 7)) +
  theme(axis.text.y = element_text(size = 12)) + 
  labs(title = "HMP2", subtitle = "Maaslin2 Effect Sizes - Raw Diet") +
  scale_x_discrete(label=function(x) str_trunc(x, 15, "right"))

mlvs_factor_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_maaslin2_PCs/all_results.tsv", header=T, check.names=T))

filttest = apply(mlvs_factor_results, 1, function(r) any(r %in% c("antibiotic", "consent_age")))
test_covariates = mlvs_factor_results[filttest,]
test_microbiome = mlvs_factor_results[!filttest,]

mlvs_factor_results$metadata[mlvs_factor_results$metadata == "RC1"] = "dietary pattern 1"
mlvs_factor_results$metadata[mlvs_factor_results$metadata == "RC2"] = "dietary pattern 2"
mlvs_factor_results$metadata[mlvs_factor_results$metadata == "consent_age"] = "consent age"

mlvs_factor = ggplot(mlvs_factor_results,aes(x=reorder(metadata,-abs(coef)), y=abs(coef))) +
  geom_boxplot() +
  labs(y=expression(beta~'coefficients'), x="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.text.x = element_text(size = 9)) +
  theme(axis.text.y = element_text(size = 12)) + 
  labs(title = "MLVS", subtitle = "Maaslin2 Effect Sizes - PC Scores") +
  scale_y_continuous(limits = c(0, 3))

#boxplot for mlvs maaslin effect sizes - raw diet
mlvs_raw_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mlvs_maaslin2_raw/all_results.tsv", header=T, check.names=T))

filttest = apply(mlvs_raw_results, 1, function(r) any(r %in% c("antibiotic", "consent_age")))
test_covariates = mlvs_raw_results[filttest,]
test_microbiome = mlvs_raw_results[!filttest,]

mlvs_raw = ggplot(mlvs_raw_results,aes(x=reorder(metadata,-abs(coef)), y=abs(coef))) +
  geom_boxplot() +
  labs(y=expression(beta~'coefficients'), x="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.text.x = element_text(size = 7)) +
  theme(axis.text.y = element_text(size = 12)) + 
  labs(title = "MLVS", subtitle = "Maaslin2 Effect Sizes - Raw Diet") +
  scale_y_continuous(limits = c(0, 3))

### combining plots ###
westernplot = ggarrange(hmp2_factor,hmp2_raw,mlvs_factor,mlvs_raw, labels = c("A","B","C","D"),nrow = 2, ncol=2,align = "hv")
ggsave(
  "figure2/WesternAdult_PCVsRawBox.png",
  westernplot,
  width = 8,
  height = 8,
  dpi = 1200
)
