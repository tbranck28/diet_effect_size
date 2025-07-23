#########################
### Figure 4 plotting ###
#########################
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

#set directory paths
input_dir = file.path("/home","tobynb","diet_test","input","Faith_mice","Intermediate",fsep = "/")
output_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")

### Figure 4C: E. rectale vs diet (food items) ###
adonis_result = as.data.frame(data.table::fread(paste0(output_dir,"/mouse_complex_adonis.csv"), header=T, check.names=T))

complex_diet_adonis = ggplot(adonis_result, aes(x=reorder(diet_variable,R2), y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() +
  geom_point(data = adonis_result[adonis_result$p_adj < 0.1, ],aes(x=reorder(diet_variable,R2), R2 + 0.01), shape = "*", size=4.233, color="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(expression("R"^2)) +
  xlab("Diet variable")

cairo_pdf(file = paste0(output_dir,"/adonis_complex_diet.pdf"), width=8, height = 7)
complex_diet_adonis
dev.off()

### Figure 4C: E. rectale vs diet (food items) ###
#import taxonomy table, metadata, and maaslin beta coefficient results
filttest_complex = as.data.frame(data.table::fread(paste0(input_dir,"/mouse_complex_species.csv"), header=T, check.names=T))
metadata_input_complex = as.data.frame(data.table::fread(paste0(output_dir,"/mouse_complex_metadata_binned.csv"), header=T, check.names=T))
vec = as.data.frame(data.table::fread(paste0(output_dir,"/mouse_complexdiet_maaslin2_binned/all_results.tsv"), header=T, check.names=T))
vec = vec$coef
rownames(filttest_complex) = filttest_complex$V1
filttest_complex$V1 = NULL
rownames(metadata_input_complex) = metadata_input_complex$V1
metadata_input_complex$V1 = NULL
#merge into one data for ease of plotting
erec_df = merge(filttest_complex, metadata_input_complex, by=0, all=TRUE)

# E. rectale vs. apple sauce box plot
eub_applesauce = ggplot(erec_df, aes(x=factor(apple_sauce, levels=c("None","Mid","High")), y=Eubacterium.rectale)) + 
  geom_boxplot(outlier.shape = NA,lwd=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("applesauce (g/kg)") + ylab("")+
  theme_classic(base_size = 10) + 
  geom_jitter(shape=16, position=position_jitter(0.3), alpha=.3, size = 3)

# E. rectale vs. peaches box plot
eub_peaches = ggplot(erec_df, aes(x=factor(peaches, levels=c("None","Mid","High")), y=Eubacterium.rectale)) + 
  geom_boxplot(outlier.shape = NA,lwd=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("peaches (g/kg)") + ylab("") +
  theme_classic(base_size = 10) + 
  geom_jitter(shape=16, position=position_jitter(0.3), alpha=.3, size = 3)

# E. rectale vs. chicken box plot
eub_chicken = ggplot(erec_df, aes(x=factor(chicken, levels=c("None","Mid","High")), y=Eubacterium.rectale)) + 
  geom_boxplot(outlier.shape = NA,lwd=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("chicken (g/kg)") + ylab("") +
  theme_classic(base_size = 10) + 
  geom_jitter(shape=16, position=position_jitter(0.3), alpha=.3,size = 3)

vec = as.data.frame(vec)
vec_box = ggplot(abs(vec),aes(y=vec)) + 
  geom_boxplot() +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Mouse diet~microbiome") + ylab("Beta coefficients (absolute value)")

### Figure 4B: beta coefficient box plots for all hosts
#bring in all of the beta coefficient info two plot distribution of beta coefficients
for_boxplot = as.data.frame(data.table::fread(paste0(output_dir,"/beta_boxplot_df.csv"), header=T, check.names=T)) #made in the Figure2_PCs.R script
for_boxplot$V1 = NULL
nhp_sig_results = as.data.frame(data.table::fread(paste0(output_dir,"/nonhumanPrimates_maaslin2_output/all_results.tsv"), header=T, check.names=T))
nhp_sig_results$metadata[nhp_sig_results$metadata == 'Host_diet'] = 'NHP_Host_diet'
mouse_sig_results = as.data.frame(data.table::fread(paste0(output_dir,"/mouse_complexdiet_maaslin2_binned/all_results.tsv"), header=T, check.names=T))
mouse_sig_results$metadata <- "Mouse"

infant_sig_results = as.data.frame(data.table::fread(paste0(output_dir,"/infants_maaslin_combined.csv"), header=T, check.names=T))
infant_sig_results$V1 = NULL
infant_sig_results$Study = NULL
infant_sig_results$metadata[1:210] = 'infants_Backhed'
infant_sig_results$metadata[211:322] = 'infants_Murphy'

for_boxplot = rbind(for_boxplot,nhp_sig_results,mouse_sig_results,infant_sig_results)
for_boxplot$metadata = gsub("Dietary","dietary",for_boxplot$metadata)
beta_coeff_box = ggplot(for_boxplot,aes(x=reorder(metadata,-abs(coef)), y=abs(coef),fill=metadata)) +
  geom_boxplot() +
  labs(y=expression(beta~'coefficients'), x="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values=c("darkseagreen","indianred3","indianred3","azure3","azure3","lemonchiffon3","lemonchiffon3","darkorchid1","coral1","azure4","azure4","lightsalmon4","orange2")) +
  theme(legend.position=c(.85,.6)) +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 10)) +
  theme(legend.title=element_blank()) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(legend.text=element_text(size=10)) +
  theme(legend.key=element_blank()) +
  theme(legend.position="none")

cairo_pdf(paste0(output_dir,"/betas.pdf"),width=5,height=4)
beta_coeff_box
dev.off()

top = ggarrange(complex_diet_adonis,beta_coeff_box, labels=c("A","B"))
boxes = ggarrange(eub_applesauce,eub_peaches,eub_chicken,nrow=1,ncol=3,labels=c("C","",""))
boxes = annotate_figure(boxes, left = textGrob("Eubacterium rectale",rot = 90, vjust = 1, gp = gpar(fontsize = 8,cex = 1.3)))
pdf(paste0(output_dir,"/Figure_4.pdf"), width=8, height = 6)
ggarrange(top,boxes, nrow = 2,ncol=1,widths = c(6,4))
dev.off()

cairo_pdf(file = paste0(output_dir,"/Erectale_boxplots.pdf"), width=8, height = 7)
boxes
dev.off()