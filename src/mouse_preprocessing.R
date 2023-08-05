setwd("/Users/tobynbranck/Documents/Effect_size/data/")
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
library("pheatmap")
library("ape")
library("RColorBrewer")
library("ggpubr")
library(MicEco)

test = data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/mouse_faith/metaphlan_taxonomic_profiles.tsv", header=T, check.names = F)
names(test) <- gsub("# ", "", names(test), fixed = TRUE)
test = test[grep("s__", test$taxonomy),]
names(test) = gsub(pattern = "_.*", replacement = "", x = names(test))
test$taxonomy <- gsub('^.*?s__', '', test$taxonomy)
test$taxonomy = gsub("_", " ", test$taxonomy)

#plot to observe distribution of sum of abundance by sample - for samples with < 100, how to handle? why are some so low (in 80-98 range)??
qplot(colSums(test[,-1]), geom="histogram")
#sample SRR097225 has all zeros
test$SRR097225 = NULL

#divide by 100 to put on 0-1 scale
rownames(test) = test$taxonomy
test$taxonomy = NULL
bugnames = rownames(test)
test = test/100
rownames(test) = bugnames

#less strict filter for prevalence / abundance; 10^-5 in at least 25% of the samples
test$taxonomy = bugnames
filttest = test[rowSums(test > 0.00001) >= dim(test)[2]*.25, ]
filttest = data.frame(filttest)
rownames(filttest) = filttest$taxonomy
filttest$taxonomy = NULL
names(filttest) <- sub("^X", "", names(filttest))

filttest = as.data.frame(t(filttest))
samplist = rownames(filttest)
write.table(filttest, "/Users/tobynbranck/Documents/Effect_size/analysis/mouse_species.csv", sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

#stricter filter for prevalence / abundance; 10^-5 in at least 60% of the samples
test$taxonomy = bugnames
filttest_strict = test[rowSums(test > 0.00001) >= dim(test)[2]*.25, ]
filttest_strict = data.frame(filttest_strict)
rownames(filttest_strict) = filttest_strict$taxonomy
filttest_strict$taxonomy = NULL
names(filttest_strict) <- sub("^X", "", names(filttest_strict))

filttest_strict = as.data.frame(t(filttest_strict))
samplist = rownames(filttest_strict)
#write.table(filttest_strict, "/Users/tobynbranck/Documents/Effect_size/analysis/mouse_species_strict.csv", sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

######### sample download check --> all samples I need for all three sub-experiments are present and have 
######### been processed through bb3 (they are included in the taxonomic profiles). NOTE: simple_diet_1, 
######### simple_diet_2, complex_diet samples to make lists were taken from the ipynb "mouse_dataset_samples.ipynb"
samps_ran_bb3 = colnames(test[,-402])
meta = data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/mouse_faith/mouse_metadata.csv", header=T, check.names = F)
#keep samples with at least 500,000 depth
kneadd = data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/mouse_faith/kneaddata_read_count_table.tsv", header=T, check.names = F)
kneadd = kneadd[kneadd$`final single`>=500000,]
meta = meta[meta$Run %in% kneadd$Sample,]

complex_samps_df = dplyr::filter(meta, grepl("puree",diet))
complex_samps = complex_samps_df$Run

simple_diet_1_samples = 

simple_diet_1_samples = c('SRR097152', 'SRR097140', 'SRR097029', 'SRR097052', 'SRR096959',
'SRR097037', 'SRR097100', 'SRR097189', 'SRR097123', 'SRR097007',
'SRR097098', 'SRR097186', 'SRR096987', 'SRR097121', 'SRR097125',
'SRR097218', 'SRR097017', 'SRR097212', 'SRR096964', 'SRR097040',
'SRR096986', 'SRR097118', 'SRR097126', 'SRR097185', 'SRR097021',
'SRR096998', 'SRR096984', 'SRR097055', 'SRR097171',
'SRR097159', 'SRR097164', 'SRR096999', 'SRR097117', 'SRR097051',
'SRR097104', 'SRR097206', 'SRR097004', 'SRR097168', 'SRR097011',
'SRR097173', 'SRR097128', 'SRR097054', 'SRR097063', 'SRR097012',
'SRR097046', 'SRR097081', 'SRR097131', 'SRR097088', 'SRR097127',
'SRR097224', 'SRR097187', 'SRR097103', 'SRR097036', 'SRR097072',
'SRR097001', 'SRR097084', 'SRR097060', 'SRR097149', 'SRR097077',
'SRR097120', 'SRR096966', 'SRR097157', 'SRR097153', 'SRR097167',
'SRR096977', 'SRR097064', 'SRR097018', 'SRR097003', 'SRR097135',
'SRR097044', 'SRR096973', 'SRR097034', 'SRR096980', 'SRR097090',
'SRR097095', 'SRR097217', 'SRR097151', 'SRR097156', 'SRR097142')
simple_diet_2_samples = c('SRR097137', 'SRR097107', 'SRR097022', 'SRR097089', 'SRR097058',
                  'SRR097035', 'SRR097148', 'SRR096972', 'SRR097030', 'SRR097216',
                  'SRR097170', 'SRR097028', 'SRR097009', 'SRR096960', 'SRR097075',
                  'SRR097122', 'SRR097079', 'SRR096975', 'SRR097066', 'SRR097010',
                  'SRR096965', 'SRR097027', 'SRR096952', 'SRR097130', 'SRR097014',
                  'SRR096968', 'SRR097039', 'SRR097210', 'SRR097158', 'SRR097042',
                  'SRR097177', 'SRR096982', 'SRR097183', 'SRR097190', 'SRR097019',
                  'SRR097144', 'SRR097002', 'SRR097094', 'SRR097179', 'SRR097143',
                  'SRR097209', 'SRR097182', 'SRR097053', 'SRR096978', 'SRR097091',
                  'SRR097133', 'SRR097047', 'SRR097188', 'SRR097023', 'SRR097099',
                  'SRR097065', 'SRR097146', 'SRR096985', 'SRR097221', 'SRR097163',
                  'SRR097176', 'SRR097191', 'SRR097102', 'SRR097175', 'SRR097222',
                  'SRR097083', 'SRR097067', 'SRR097180', 'SRR096981', 'SRR097154',
                  'SRR096969', 'SRR096979', 'SRR097106', 'SRR097124', 'SRR097033',
                  'SRR097038', 'SRR097174', 'SRR097172', 'SRR097097', 'SRR097005',
                  'SRR097031', 'SRR097049', 'SRR097207', 'SRR097076', 'SRR097043',
                  'SRR096956', 'SRR096958', 'SRR097000', 'SRR096974', 'SRR097068',
                  'SRR097093')
complex_diet_samples = c('SRR097276', 'SRR097288', 'SRR097301', 'SRR097295', 'SRR097269',
                 'SRR097347', 'SRR097275', 'SRR097316', 'SRR097320', 'SRR097389',
                 'SRR097317', 'SRR097279', 'SRR097324', 'SRR097354', 'SRR097250',
                 'SRR097322', 'SRR097314', 'SRR097311', 'SRR097266', 'SRR097310',
                 'SRR097345', 'SRR097390', 'SRR097293', 'SRR097306', 'SRR097248',
                 'SRR097303', 'SRR097277', 'SRR097304', 'SRR097308', 'SRR097376',
                 'SRR097330', 'SRR097391', 'SRR097366', 'SRR097361', 'SRR097388',
                 'SRR097312', 'SRR097290', 'SRR097365', 'SRR097397', 'SRR097319',
                 'SRR097369', 'SRR097278', 'SRR097328', 'SRR097318', 'SRR097283',
                 'SRR097298', 'SRR097291', 'SRR097254', 'SRR097326', 'SRR097331',
                 'SRR097380', 'SRR097273', 'SRR097378', 'SRR097370', 'SRR097329',
                 'SRR097350', 'SRR097244', 'SRR097399', 'SRR097287', 'SRR097313',
                 'SRR097252', 'SRR097302')

###########################################
###### dietary metadata curation ##########
###########################################
simple_diet_exp1 = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mouse_simplediet_exp1.csv", header = TRUE, check.names = F))
simple_diet_exp2 = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mouse_simplediet_exp2.csv", header = TRUE, check.names = F))
complex_diet = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mouse_complexdiet.csv", header = TRUE, check.names = F))
simple_diets = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/mouse_faith/mouse_simple_diets.txt", sep = '\t', header = TRUE, check.names = F))
complex_diets = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/mouse_faith/mouse_complex_diet.txt", sep = '\t', header = TRUE, check.names = F))

melted = melt(complex_diets, id=c("V1","g/kg"))

melted = melted[1:232,]
melted$'apple sauce' <- ifelse(melted$value == 'apple sauce', 1, 0)
melted$'peaches' <- ifelse(melted$value == 'peaches', 1, 0)
melted$'chicken' <- ifelse(melted$value == 'chicken', 1, 0)
melted$'sweet potatoes' <- ifelse(melted$value == 'sweet potatoes', 1, 0)
melted$'oatmeal' <- ifelse(melted$value == 'oatmeal', 1, 0)
melted$'peas' <- ifelse(melted$value == 'peas', 1, 0)
melted$'rice' <- ifelse(melted$value == 'rice', 1, 0)
melted$'beef' <- ifelse(melted$value == 'beef', 1, 0)
melted$'peas' <- ifelse(melted$value == 'peas', 1, 0)
melted$value = NULL

melted = dcast(melted, variable + V1 ~ value)
#remove even number weeks (these are the "reset" weeks where mice were given the refined diet to stabilize 
#the microbiome)
melted = melted[!grepl("week2", melted$V1),]
melted = melted[!grepl("week4", melted$V1),]
melted = melted[!grepl("week6", melted$V1),]
melted = melted[!grepl("week8", melted$V1),]
melted = melted[!grepl("week10", melted$V1),]
melted$value = NULL

# replace 1 in the row by the g/kg value
melted = melted %>% 
  dplyr::mutate_if(. %in% melted[,4:11], ~ dplyr::recode(., "1" = melted$`g/kg`, "0" = '0'))
melted[,4:11] <- sapply(melted[,4:11],as.numeric)

#collapse so that there is 1 row per week, per mouse
complex_diets = aggregate(cbind(melted$`apple sauce`,melted$peaches,melted$chicken,melted$`sweet potatoes`,melted$oatmeal,melted$peas,melted$rice,melted$beef) ~ V1 + variable, data = melted, FUN = sum, na.rm = TRUE)
names(complex_diets) = c("Week","Mouse","apple sauce","peaches","chicken","sweet potatoes","oatmeal","peas","rice","beef")

#now add sample name
complex_diet$diet_period[complex_diet$diet_period == 'E9'] <- 'week1'
complex_diet$diet_period[complex_diet$diet_period == 'E11'] <- 'week3'
complex_diet$diet_period[complex_diet$diet_period == 'E13'] <- 'week5'
complex_diet$diet_period[complex_diet$diet_period == 'E15'] <- 'week7'
complex_diet$diet_period[complex_diet$diet_period == 'E17'] <- 'week9'
complex_diet$diet_period[complex_diet$diet_period == 'E19'] <- 'week11'

complex_diet$Individual[complex_diet$Individual == 'm1'] <- 'mouse1'
complex_diet$Individual[complex_diet$Individual == 'm2'] <- 'mouse2'
complex_diet$Individual[complex_diet$Individual == 'm3'] <- 'mouse3'
complex_diet$Individual[complex_diet$Individual == 'm4'] <- 'mouse4'
complex_diet$Individual[complex_diet$Individual == 'm5'] <- 'mouse5'
complex_diet$Individual[complex_diet$Individual == 'm6'] <- 'mouse6'
complex_diet$Individual[complex_diet$Individual == 'm7'] <- 'mouse7'
complex_diet$Individual[complex_diet$Individual == 'm8'] <- 'mouse8'

# merge individual and diet week columns to be able to get sample name 
complex_diet$identifier <- paste(complex_diet$Individual,complex_diet$diet_period)
complex_diets$identifier = paste(complex_diets$Mouse, complex_diets$Week)

#complex_diet$Run<-complex_diets$Run[match(complex_diet$identifier,complex_diets$identifier)]
complex_diet$'apple sauce'<-complex_diets$`apple sauce`[match(complex_diet$identifier,complex_diets$identifier)]
complex_diet$'peaches'<-complex_diets$peaches[match(complex_diet$identifier,complex_diets$identifier)]
complex_diet$'chicken'<-complex_diets$chicken[match(complex_diet$identifier,complex_diets$identifier)]
complex_diet$'sweet potatoes'<-complex_diets$`sweet potatoes`[match(complex_diet$identifier,complex_diets$identifier)]
complex_diet$'oatmeal'<-complex_diets$oatmeal[match(complex_diet$identifier,complex_diets$identifier)]
complex_diet$'peas'<-complex_diets$peas[match(complex_diet$identifier,complex_diets$identifier)]
complex_diet$'rice'<-complex_diets$rice[match(complex_diet$identifier,complex_diets$identifier)]
complex_diet$'beef'<-complex_diets$beef[match(complex_diet$identifier,complex_diets$identifier)]

complex_diet_metadata = complex_diet[,c("Run","Individual","sampling_day","diet_period","apple sauce","peaches","chicken","sweet potatoes","oatmeal","peas","rice","beef")]
rownames(complex_diet_metadata) = complex_diet_metadata$Run
complex_diet_metadata$Run = NULL

#get common samples in filttest and complex_diet_metadata
filttest_complex = filttest[match(row.names(complex_diet_metadata), row.names(filttest)), ]
filttest_strict_complex = filttest_strict[match(row.names(complex_diet_metadata), row.names(filttest_strict)), ]
setdiff(rownames(filttest_complex), rownames(complex_diet_metadata))
#### plot diet bar plots ####
for_bar = melt(complex_diet_metadata, id=c("Individual","sampling_day","diet_period"))
for_bar$concat = paste(for_bar$Individual,for_bar$diet_period,for_bar$variable)
for_bar = for_bar[order(for_bar[,'Individual']),]
for_bar = for_bar[!duplicated(for_bar$concat),]
for_bar$sampling_day = NULL
for_bar$concat = paste(for_bar$Individual,for_bar$diet_period)
colors = c( "#A54657",  "#582630", "#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A","#86BBD8")
pdf(file = "/Users/tobynbranck/Documents/Effect_size/analysis/mouse_complex_diet_barplots.pdf",width=12,height=9)
ggplot(for_bar, aes(x = factor(diet_period, levels=c("week1","week3","week5","week7","week9","week11")), fill = variable, y = value)) + 
  geom_bar(stat = "identity", colour = "black") + 
  theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1, face= "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold", colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "g/kg", fill = "diet variable") + 
  facet_wrap(~Individual) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

#### plot heatmap of complex diet variables ####
mat_col = data.frame(Individual = complex_diet_metadata$Individual)
rownames(mat_col) = rownames(complex_diet_metadata)

mat_colors = list(Individual = brewer.pal(8,"Set3"))
names(mat_colors$Individual) = unique(complex_diet_metadata$Individual)
write.csv(complex_diet_metadata,"mouse_complex_diet_forSimul.csv")
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
z_diets_complex = as.data.frame(scale(complex_diet_metadata[,4:11]))
z_diets_complex = do.call(data.frame,lapply(z_diets_complex, function(x) replace(x, is.infinite(x),0)))
rownames(z_diets_complex) = rownames(complex_diet_metadata)
## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))
cairo_pdf(file = "/Users/tobynbranck/Documents/Effect_size/analysis/figure4/complex_diet_heatmap.pdf", width=8, height = 7)
pheatmap(z_diets_complex, annotation_row = mat_col, annotation_colors = mat_colors, border_color = NA,show_rownames = FALSE,annotation_names_row=FALSE, color = brewer.pal(6,"Purples"))
dev.off()

###plotting##################################################
#pcoa
mouse_dist_matrix = vegdist(filttest_complex, method="bray",diag=TRUE)
mouse_pcoa = pcoa(mouse_dist_matrix)
mouse_pcoa_df = as.data.frame(mouse_pcoa$vectors[,1:2])
#overlaying ordinations by covariates
for (col in seq_along(complex_diet_metadata)){
  print(col)
  var = names(complex_diet_metadata[col])
  print(var)
  plot = ggplot(mouse_pcoa_df, aes(x = Axis.1, y = Axis.2, color=factor(complex_diet_metadata[[col]]))) + 
    geom_point(alpha = 5/10, size=3) +
    labs(colour = "g/kg") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(x=paste0("PC1: ",round(mouse_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
         y=paste0("PC2: ",round(mouse_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
    ggtitle(paste("complex diet; colored by",var)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = c(.9, .25)) +
    theme(legend.key.size = unit(0.4, "cm")) +
    theme(legend.title = element_text(size = 8), 
          legend.text = element_text(size = 8))
  ggsave(plot, file=paste0("/Users/tobynbranck/Documents/Effect_size/analysis/mouse_ords_overlayed_covariates/mouse_complex_ord_", var,".png"), width = 14, height = 10, units = "cm")
}

pdf(file = "/Users/tobynbranck/Documents/Effect_size/analysis/mouse_sampsummary_ord.pdf", width=10, height = 12)
ggplot(mouse_pcoa_df, aes(x = Axis.1, y = Axis.2, color=as.factor(complex_diet_metadata$Individual))) + 
  #geom_point(alpha = 5/10) +
  #labs(colour = "Individual") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Mouse cohort sample summary") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(.9, .25)) + 
  geom_line(aes(group = interaction(complex_diet_metadata$Individual, complex_diet_metadata$diet_period)),color="black",size=.3) +
  geom_point(aes(shape = factor(complex_diet_metadata$sampling_day)), size = 3) +
  guides(col=guide_legend("Individual"),
         shape=guide_legend("Sampling Day")) +
  labs(x=paste0("PC1: ",round(mouse_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(mouse_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
  geom_polygon(aes(group = interaction(complex_diet_metadata$Individual, complex_diet_metadata$diet_period)),alpha=0.2) + theme(legend.title = element_text(size = 7),                                                                                                                                 legend.text = element_text(size = 7))
dev.off()  

#to show distribution of individuals' samples across diets and days after diet switch
ggplot(complex_diet_metadata, aes(y=interaction(Individual,diet_period), x=as.factor(sampling_day)), color=as.factor(complex_diet_metadata$Individual)) + 
  geom_point(stat="identity", aes(colour = factor(Individual))) +
  scale_fill_viridis(discrete = T, option = "E") +
  ggtitle("") +
  theme(legend.position="none") +
  xlab("")

numSamps_perInd = complex_diet_metadata %>% count("Individual")
mean(numSamps_perInd$freq)
sum(numSamps_perInd$freq)
###################################################
######complex diet analysis #######################
#adonis
adonis_output = vector()
adonis_res_rsq = vector()
adonis_res_pval = vector()
#complex_diet_metadata[,4:11] <- lapply(complex_diet_metadata[,4:11], as.character)
complex_diet_metadata_adon = complex_diet_metadata[,4:11]
for (col in seq_along(complex_diet_metadata_adon)){
  adonis.univ = adonis(filttest_complex ~ complex_diet_metadata_adon[[col]], data = complex_diet_metadata_adon, method = "bray", strata = complex_diet_metadata$Individual)
  print(adonis_OmegaSq(adonis.univ, partial = TRUE))
  #adonis.univ = adonis(filttest_complex ~ complex_diet_metadata$Individual + complex_diet_metadata_adon[[col]], data = complex_diet_metadata_adon, method = "bray")
  print(adonis.univ)
  adonis_output[col] = adonis.univ
  adonis_res_rsq[col] = adonis.univ$aov.tab[1,]$R2
  adonis_res_pval[col] = adonis.univ$aov.tab[1,]$`Pr(>F)`}
univar_res_rgnav = rbind(adonis_res_pval, adonis_res_rsq)
univar_res_rgnav = as.data.frame(t(univar_res_rgnav))
names(univar_res_rgnav) = c("P-Value", "R2")
univar_res_rgnav$`P-Value` = as.numeric(univar_res_rgnav$`P-Value`)
univar_res_rgnav$R2 = as.numeric(univar_res_rgnav$R2)
univar_res_rgnav$p_adj = p.adjust(univar_res_rgnav$`P-Value`, "fdr")
univar_res_rgnav$stars = cut(univar_res_rgnav$p_adj, c(0, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "-", ""))
univar_res_rgnav$diet_variable = colnames(complex_diet_metadata[,4:11])
univar_res_rgnav = univar_res_rgnav[order(univar_res_rgnav$R2),]

complex_diet_adonis = ggplot(univar_res_rgnav, aes(x=reorder(diet_variable,R2), y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() +
  geom_point(data = univar_res_rgnav[univar_res_rgnav$p_adj < 0.05, ],aes(x=reorder(diet_variable,R2), R2 + 0.01), shape = "*", size=4.233, color="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  #ggtitle("Mouse cohort adonis") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(expression("R"^2)) +
  xlab("Diet variable")

#saving complex_mouse metadata
write.table(complex_diet_metadata,"/Users/tobynbranck/Documents/Effect_size/analysis/mouse_complex_metadata.csv", sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)
#write.table(filttest_complex,"/Users/tobynbranck/Documents/Effect_size/analysis/mouse_complex_species.csv", sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

#maaslin2
metadata_input_complex = complex_diet_metadata
metadata_input_complex$sampling_day = NULL
metadata_input_complex$diet_period = NULL
colnames(metadata_input_complex)[2] <- "apple_sauce"
colnames(metadata_input_complex)[5] <- "sweet_potatoes"

Maaslin2(input_data = filttest_strict_complex, input_metadata = metadata_input_complex, output = "/Users/tobynbranck/Documents/Effect_size/analysis/mouse_complexdiet_maaslin2", random_effects = c("Individual"))

metadata_input_complex[] <- lapply(metadata_input_complex, as.character)
Maaslin2(input_data = filttest_strict_complex, input_metadata = metadata_input_complex, output = "/Users/tobynbranck/Documents/Effect_size/analysis/mouse_complexdiet_maaslin2_categorical", random_effects = c("Individual"))

metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "666.7", replacement = "High", fixed = TRUE)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "421.1", replacement = "High", fixed = TRUE)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "250", replacement = "Mid", fixed = TRUE)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "222.2", replacement = "Mid", fixed = TRUE)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "105.3", replacement = "Mid", fixed = TRUE)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "55.6", replacement = "Mid", fixed = TRUE)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "52.6", replacement = "Mid", fixed = TRUE)
metadata_input_complex[] <- lapply(metadata_input_complex, gsub, pattern = "0", replacement = "None", fixed = TRUE)
Maaslin2(input_data = filttest_strict_complex, input_metadata = metadata_input_complex, output = "/Users/tobynbranck/Documents/Effect_size/analysis/mouse_complexdiet_maaslin2_binned", random_effects = c("Individual"), reference=c("apple_sauce,None","peaches,None","chicken,None","sweet_potatoes,None","oatmeal,None","peas,None","rice,None","beef,None"))

adonis_output = vector()
adonis_res_rsq = vector()
adonis_res_pval = vector()
metadata_input_complex_adon = metadata_input_complex[,2:9]
for (col in seq_along(metadata_input_complex_adon)){
  print(metadata_input_complex_adon[col])
  adonis.univ = adonis(filttest_complex ~ metadata_input_complex_adon[[col]], data = metadata_input_complex_adon, method = "bray", strata = metadata_input_complex$Individual)
  print(adonis_OmegaSq(adonis.univ, partial = TRUE))
  adonis_output[col] = adonis.univ
  adonis_res_rsq[col] = adonis.univ$aov.tab[1,]$R2
  adonis_res_pval[col] = adonis.univ$aov.tab[1,]$`Pr(>F)`}
univar_res_rgnav = rbind(adonis_res_pval, adonis_res_rsq)
univar_res_rgnav = as.data.frame(t(univar_res_rgnav))
names(univar_res_rgnav) = c("P-Value", "R2")
univar_res_rgnav$`P-Value` = as.numeric(univar_res_rgnav$`P-Value`)
univar_res_rgnav$R2 = as.numeric(univar_res_rgnav$R2)
univar_res_rgnav$p_adj = p.adjust(univar_res_rgnav$`P-Value`, "fdr")
univar_res_rgnav$stars = cut(univar_res_rgnav$p_adj, c(0, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "-", ""))
univar_res_rgnav$diet_variable = colnames(metadata_input_complex[,2:9])
univar_res_rgnav = univar_res_rgnav[order(univar_res_rgnav$R2),]

complex_diet_adonis = ggplot(univar_res_rgnav, aes(x=reorder(diet_variable,R2), y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() +
  geom_point(data = univar_res_rgnav[univar_res_rgnav$p_adj < 0.1, ],aes(x=reorder(diet_variable,R2), R2 + 0.01), shape = "*", size=4.233, color="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  #ggtitle("Mouse cohort adonis") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(expression("R"^2)) +
  xlab("Diet variable")

cairo_pdf(file = "/Users/tobynbranck/Documents/Effect_size/analysis/figure4/adonis_complex_diet.pdf", width=8, height = 7)
ggplot(univar_res_rgnav, aes(x=reorder(diet_variable,R2), y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() +
  geom_point(data = univar_res_rgnav[univar_res_rgnav$p_adj < 0.1, ],aes(x=reorder(diet_variable,R2), R2 + 0.01), shape = "*", size=4.233, color="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  #ggtitle("Mouse cohort adonis") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(expression("R"^2)) +
  xlab("Diet variable")
dev.off()

### E. rectale plots ###
vec = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mouse_complexdiet_maaslin2_binned/all_results.tsv", header=T, check.names=T))
vec = vec$coef
erec_df = merge(filttest_strict_complex, metadata_input_complex, by=0, all=TRUE)

eub_applesauce = ggplot(erec_df, aes(x=factor(erec_df$apple_sauce, levels=c("None","Mid","High")), y=erec_df$`Eubacterium rectale`)) + 
  geom_boxplot(outlier.shape = NA,lwd=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("applesauce (g/kg)") + ylab("")+
  theme_classic(base_size = 10) + 
  geom_jitter(shape=16, position=position_jitter(0.3), alpha=.3, size = 3)

eub_peaches = ggplot(erec_df, aes(x=factor(erec_df$peaches, levels=c("None","Mid","High")), y=erec_df$`Eubacterium rectale`)) + 
  geom_boxplot(outlier.shape = NA,lwd=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("peaches (g/kg)") + ylab("") +
  theme_classic(base_size = 10) + 
  geom_jitter(shape=16, position=position_jitter(0.3), alpha=.3, size = 3)

eub_chicken = ggplot(erec_df, aes(x=factor(erec_df$chicken, levels=c("None","Mid","High")), y=erec_df$`Eubacterium rectale`)) + 
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

for_boxplot = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/beta_boxplot_df.csv", header=T, check.names=T))
for_boxplot$V1 = NULL
nhp_sig_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/dada2_test/nonhumanPrimates_maaslin2_output/all_results.tsv", header=T, check.names=T))
nhp_sig_results$metadata[nhp_sig_results$metadata == 'Host_diet'] = 'NHP_Host_diet'
mouse_sig_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mouse_complexdiet_maaslin2_binned/all_results.tsv", header=T, check.names=T))
mouse_sig_results$metadata <- "Mouse"

infant_sig_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/infants_maaslin_combined.csv", header=T, check.names=T))
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

cairo_pdf("/Users/tobynbranck/Documents/Effect_size/analysis/figure4/betas.pdf",width=5,height=4)
beta_coeff_box
dev.off()

top = ggarrange(complex_diet_adonis,beta_coeff_box, labels=c("A","B"))
boxes = ggarrange(eub_applesauce,eub_peaches,eub_chicken,nrow=1,ncol=3,labels=c("C","",""))
boxes = annotate_figure(boxes, left = textGrob("Eubacterium rectale",rot = 90, vjust = 1, gp = gpar(fontsize = 8,cex = 1.3)))
pdf("/Users/tobynbranck/Documents/Effect_size/analysis/figure4/Figure_4.pdf", width=8, height = 6)
ggarrange(top,boxes, nrow = 2,ncol=1,widths = c(6,4))
dev.off()

cairo_pdf(file = "/Users/tobynbranck/Documents/Effect_size/analysis/figure4/Erectale_boxplots.pdf", width=8, height = 7)
boxes
dev.off()
##############################################################################
####### Below is the code used for the refined diets shown in supplement #####
##############################################################################
simple_diets = as.data.frame(t(simple_diets))
colnames(simple_diets) = simple_diets[1,]
simple_diets = simple_diets[-1,]
simple_diets = simple_diets[,-1]
simple_diets$'diet_type' = rownames(simple_diets)

simple_diet1 = merge(simple_diet_exp1, simple_diets, by = 'diet_type')
simple_diet1_micro = simple_diet1[,1:18]
simple_diet1_macro = cbind(simple_diet1[,1:5],simple_diet1[,19:25])
rownames(simple_diet1_micro) = simple_diet1_micro$Run
#remove "SRR097225"#
simple_diet1_micro = simple_diet1_micro[row.names(simple_diet1_micro) != "SRR097225", , drop = FALSE]

simple_diet1_micro$Run = NULL
rownames(simple_diet1_macro) = simple_diet1_macro$Run
simple_diet1_macro$Run = NULL

simple_diet2 = merge(simple_diet_exp2, simple_diets, by = 'diet_type')
simple_diet2_micro = simple_diet2[,1:18]
simple_diet2_macro = cbind(simple_diet2[,1:5],simple_diet2[,19:25])
rownames(simple_diet2_micro) = simple_diet2_micro$Run
simple_diet2_micro$Run = NULL
rownames(simple_diet2_macro) = simple_diet2_macro$Run
simple_diet2_macro$Run = NULL

###########pcoa testing###############
simple_diets_combined = as.data.frame(c(rownames(simple_diet1_micro), rownames(simple_diet2_micro)))
rownames(simple_diets_combined) = simple_diets_combined$`c(rownames(simple_diet1_micro), rownames(simple_diet2_micro))`
simple_diets_combined$`c(rownames(simple_diet1_micro), rownames(simple_diet2_micro))` = NULL
simple_diets_combined$Exp = NA
simple_diets_combined$Exp[1:80] = "exper1"
simple_diets_combined$Exp[80:165] = "exper2"
individual = c(simple_diet1_micro$Individual,simple_diet2_micro$Individual)
simple_diets_combined$individual = individual

spec = filttest[match(row.names(simple_diets_combined), row.names(filttest)), ]

spec = spec[row.names(spec) != "SRR097225", , drop = FALSE]
simple_diets_combined = simple_diets_combined[row.names(simple_diets_combined) != "SRR097225", , drop = FALSE]
spec = filttest[match(row.names(simple_diets_combined), row.names(filttest)), ]

spec[is.na(spec)] <- 0
dist_matrix = vegdist(spec, method="bray",diag=TRUE)
pcoa = pcoa(dist_matrix)
pcoa_df = as.data.frame(pcoa$vectors[,1:2])

pdf(file = "/Users/tobynbranck/Documents/Effect_size/analysis/figure4/simpleExperCombined.pdf", width=6, height = 6)
ggplot(as.data.frame(pcoa_df), aes(x = Axis.1, y = Axis.2, color=simple_diets_combined$Exp)) + 
  geom_point(alpha = 4/10, size=3) +
  #geom_text(aes(label=rownames(simple_diets_combined)),hjust=0, vjust=0) +
  scale_color_manual(values = c('#E69F00','#56B4E9')) +
  labs(colour = "Sub-experiment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("PCoA - simple experiments") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(.9, .2)) +
  labs(x=paste0("PC1: ",round(pcoa$values$Relative_eig[1],digits=3)*100,"%"),
       y=paste0("PC2: ",round(pcoa$values$Relative_eig[2],digits=3)*100,"%"))
dev.off()

adonis_obj_subexps = adonis(spec ~ simple_diets_combined$Exp, data = simple_diets_combined, method = "bray",na.rm=TRUE, strata = simple_diets_combined$individual)

diets_combined = rbind(simple_diet1_micro,simple_diet2_micro)

diets_combined[,5:17] <- lapply(diets_combined[,5:17], as.numeric)
adonis_res_rsq_combined = vector()
adonis_res_pval_combined = vector()

diets_combined$`Maltodextrin Lo-Dex 10` = NULL
diets_combined$`Cellulose (Fiber)` = NULL
diets_combined$`Ethoxyquin (Liquid)` = NULL
diets_combined_adon = diets_combined[,5:14]

#### plot heatmap of simple combined diet variables #### - maybe split the 
#### map into macro vs. micro so that we can better tell the concentration
#### difference between the two, since the macro is masking the micro [ ]'s
mat_col = data.frame(Individual = diets_combined$Individual)
rownames(mat_col) = rownames(diets_combined)

mat_colors = brewer.pal(9,"RdPu")
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
z_diets_combined = as.data.frame(scale(diets_combined[,5:14]))
z_diets_combined = do.call(data.frame,lapply(z_diets_combined, function(x) replace(x, is.infinite(x),0)))
rownames(z_diets_combined) = rownames(diets_combined)
## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

cairo_pdf(file = "/Users/tobynbranck/Documents/Effect_size/analysis/figure4/simple_diet_heatmap.pdf", width=8, height = 7)
pheatmap(z_diets_combined, annotation_row = mat_col, border_color = NA,show_rownames = FALSE,annotation_names_row=FALSE, color = brewer.pal(9,"RdPu"))
dev.off()
#### Binning diets as opposed to running with raw values ###########
binned_diets = diets_combined
binned_diets$Casein[binned_diets$Casein<139] <- "low"
binned_diets$Casein[binned_diets$Casein>138 & binned_diets$Casein<600] <- "medium"
binned_diets$Casein[binned_diets$Casein!="low" & binned_diets$Casein!="medium"] <- "high"

binned_diets$`L-Cystine`[binned_diets$`L-Cystine`<2] <- "low"
binned_diets$`L-Cystine`[binned_diets$`L-Cystine`>2 & binned_diets$`L-Cystine`<40] <- "medium"
binned_diets$`L-Cystine`[binned_diets$`L-Cystine`!="low" & binned_diets$`L-Cystine`!="medium"] <- "high"

binned_diets$Sucrose[binned_diets$Sucrose<200] <- "low"
binned_diets$Sucrose[binned_diets$Sucrose>200 & binned_diets$Sucrose<500] <- "medium"
binned_diets$Sucrose[binned_diets$Sucrose!="low" & binned_diets$Sucrose!="medium"] <- "high"

binned_diets$`Corn Starch`[binned_diets$`Corn Starch`==0] <- "low"
binned_diets$`Corn Starch`[binned_diets$`Corn Starch`==100] <- "medium"
binned_diets$`Corn Starch`[binned_diets$`Corn Starch`==400] <- "high"

binned_diets$`Corn Oil`[binned_diets$`Corn Oil`<100] <- "low"
binned_diets$`Corn Oil`[binned_diets$`Corn Oil`!="low"] <- "high"

binned_diets$`79055 Mineral mix Ca-P Deficient`[binned_diets$`79055 Mineral mix Ca-P Deficient`<14] <- "low"
binned_diets$`79055 Mineral mix Ca-P Deficient`[binned_diets$`79055 Mineral mix Ca-P Deficient`==14.74] <- "medium"
binned_diets$`79055 Mineral mix Ca-P Deficient`[binned_diets$`79055 Mineral mix Ca-P Deficient`!="low" & binned_diets$`79055 Mineral mix Ca-P Deficient`!="medium"] <- "high"

binned_diets$`Calcium Carbonate`[binned_diets$`Calcium Carbonate`<6] <- "low"
binned_diets$`Calcium Carbonate`[binned_diets$`Calcium Carbonate`=="11.3"] <- "high"
binned_diets$`Calcium Carbonate`[binned_diets$`Calcium Carbonate`=="11.6"] <- "high"
binned_diets$`Calcium Carbonate`[binned_diets$`Calcium Carbonate`=="12.25"] <- "high"
binned_diets$`Calcium Carbonate`[binned_diets$`Calcium Carbonate`!="low" & binned_diets$`Calcium Carbonate`!= "high"] <- "medium"

binned_diets$`Calcium Phosphate`[binned_diets$`Calcium Phosphate`<4] <- "low"
binned_diets$`Calcium Phosphate`[binned_diets$`Calcium Phosphate`>4 & binned_diets$`Calcium Phosphate`<9] = "medium"
binned_diets$`Calcium Phosphate`[binned_diets$`Calcium Phosphate`!="low" & binned_diets$`Calcium Phosphate`!="medium"] <- "high"

binned_diets$`40077 Vitamin mix`[binned_diets$`40077 Vitamin mix`<16] <- "low"
binned_diets$`40077 Vitamin mix`[binned_diets$`40077 Vitamin mix`!="low"] <- "high"

binned_diets$`Choline Bitartrate`[binned_diets$`Choline Bitartrate`<2.4] <- "low"
binned_diets$`Choline Bitartrate`[binned_diets$`Choline Bitartrate`==2.40] <- "medium"
binned_diets$`Choline Bitartrate`[binned_diets$`Choline Bitartrate`==2.60] <- "high"

binned_diets_adon = binned_diets[,5:14]

for (col in seq_along(binned_diets_adon)){
  adonis_obj_subscombined = adonis(spec ~ binned_diets_adon[[col]], data = binned_diets_adon, method = "bray", strata = binned_diets$Individual,na.rm=TRUE)
  adonis_res_rsq_combined[col] = adonis_obj_subscombined$aov.tab[1,]$R2
  adonis_res_pval_combined[col] = adonis_obj_subscombined$aov.tab[1,]$`Pr(>F)`}
univar_res_rgnav_combined = rbind(adonis_res_pval_combined, adonis_res_rsq_combined)
univar_res_rgnav_combined = as.data.frame(t(univar_res_rgnav_combined))
names(univar_res_rgnav_combined) = c("P-Value", "R2")
univar_res_rgnav_combined$`P-Value` = as.numeric(univar_res_rgnav_combined$`P-Value`)
univar_res_rgnav_combined$R2 = as.numeric(univar_res_rgnav_combined$R2)
univar_res_rgnav_combined$p_adj = p.adjust(univar_res_rgnav_combined$`P-Value`, "fdr")
univar_res_rgnav_combined$stars = cut(univar_res_rgnav_combined$p_adj, c(0, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "-", ""))
univar_res_rgnav_combined$diet_variable = colnames(binned_diets_adon)
univar_res_rgnav_combined = univar_res_rgnav_combined[order(univar_res_rgnav_combined$R2),]

png("/Users/tobynbranck/Documents/Effect_size/analysis/simple_diets_combined_Binned_Analysis.png")
ggplot(univar_res_rgnav_combined, aes(x=reorder(diet_variable,R2), y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() +
  geom_point(data = univar_res_rgnav_combined[univar_res_rgnav_combined$p_adj < 0.05, ],aes(x=reorder(diet_variable,R2), R2 + 0.01), shape = "*", size=4.233, color="black") +
  geom_point(data = univar_res_rgnav_combined[univar_res_rgnav_combined$p_adj < 0.1, ],aes(x=reorder(diet_variable,R2), R2 + 0.01), shape = "-", size=4.233, color="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  #ggtitle("Mouse Refined diets combined (binned)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(expression("R^2")) +
  xlab("Diet variable") +
  theme(plot.title = element_text(size=11))
dev.off()

simp_comb_adonis_binned = ggplot(univar_res_rgnav_combined, aes(x=reorder(diet_variable,R2), y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() +
  geom_point(data = univar_res_rgnav_combined[univar_res_rgnav_combined$p_adj < 0.1, ],aes(x=reorder(diet_variable,R2), R2 + 0.01), shape = "*", size=4.233, color="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  #ggtitle("Mouse Refined diets combined (binned)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(expression("R^2")) +
  xlab("") +
  theme(plot.title = element_text(size=11))

diets_combined = binned_diets
spec_strict = spec
spec_strict$`Blautia coccoides` = NULL
spec_strict$`Lactococcus lactis` = NULL
spec_strict$`Desulfovibrio piger` = NULL
diets_combined$diet_type = NULL
diets_combined$sampling_day = NULL
diets_combined$diet_period = NULL
names(diets_combined)<-make.names(names(diets_combined),unique = TRUE)

diets_combined[,2:11] <- lapply(diets_combined[,2:11], as.character)
Maaslin2(input_data = spec_strict, input_metadata = diets_combined, output = "/Users/tobynbranck/Documents/Effect_size/analysis/mouse_simple_diets_combined_character_binned", random_effects = c("Individual"), reference=c("Casein,low","L.Cystine,low","Sucrose,low","Corn.Starch,low","Corn.Oil,low","X79055.Mineral.mix.Ca.P.Deficient,low","Calcium.Carbonate,low","Calcium.Phosphate,low","40077.Vitamin.mix,low","Choline.Bitartrate,low"))

#put the two combine simple diet figures together for the Supplement
simple_combine_binned_EffectS <- as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mouse_simple_diets_combined_character_binned/all_results.tsv", header=T, check.names=T))
complex_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/mouse_complexdiet_maaslin2_binned/all_results.tsv", header=T, check.names=T))
simple_vec = simple_combine_binned_EffectS$coef
comp_vec = complex_results$coef
simple_forplot = data.frame(group = "Refined diet", value = simple_vec)
complex_forplot = data.frame(group = "Complex diet", value = comp_vec)
plot.data = rbind(simple_forplot,complex_forplot)
box = ggplot(plot.data, aes(x=group, y=abs(value))) + geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y=expression(beta~'coefficients'), x="")

t.test(abs(simple_forplot$value), abs(complex_forplot$value))
  
cairo_pdf("/Users/tobynbranck/Documents/Effect_size/analysis/figure4/mouse_refined_combined_supplement.pdf",width = 10, height = 5)
ggarrange(simp_comb_adonis_binned, box, ncol = 2, nrow = 1,labels=c("A","B"))
dev.off()

########### Extra code #############################
########### Random Forests testing #################
#RF for combined refined diets
names(spec_strict)<-make.names(names(spec_strict),unique = TRUE)
diets_combined_rf = diets_combined[2:11]
for (col in seq_along(spec_strict)){
  bug = names(spec_strict[col])
  testdf = cbind(spec_strict[[col]],diets_combined_rf)
  names(testdf)[1] = bug
  names(testdf)[2:11] = names(diets_combined)[2:11]
  rf.object = randomForest(as.formula(paste(bug, "~ .")), data = testdf, ntree = 50, importance=TRUE)
  print(rf.object)}

########## simple diet analysis when the two simple diets are ran separately (not combined) ##########
#get common samples in filttest and simple diets
filttest_simple1 = filttest[match(row.names(simple_diet1_micro), row.names(filttest)), ]
simple_diet1_micro = simple_diet1_micro[match(row.names(filttest_simple1), row.names(simple_diet1_micro)), ]
filttest_simple1 = filttest_simple1[rowSums(is.na(filttest_simple1)) == 0, ]
simple_diet1_micro = simple_diet1_micro[rowSums(is.na(simple_diet1_micro)) == 0, ]

filttest_strict_simple1 = filttest_strict[match(row.names(simple_diet1_micro), row.names(filttest_strict)), ]

setdiff(rownames(filttest_simple1), rownames(simple_diet1_micro))
setdiff(rownames(filttest_strict_simple1), rownames(simple_diet1_micro))

filttest_simple2 = filttest[match(row.names(simple_diet2_micro), row.names(filttest)), ]
filttest_strict_simple2 = filttest_strict[match(row.names(simple_diet2_micro), row.names(filttest_strict)), ]
setdiff(rownames(filttest_simple2), rownames(simple_diet2_micro))
setdiff(rownames(filttest_strict_simple2), rownames(simple_diet2_micro))
#removing Maltodextrin, Cellulose, Ethoxyquin because each of these variables have only one value (no variation)
simple_diet1_micro$`Maltodextrin Lo-Dex 10` = NULL
simple_diet1_micro$`Cellulose (Fiber)` = NULL
simple_diet1_micro$`Ethoxyquin (Liquid)` = NULL

##### overlaying ordinations by covariates for simple diet 1 #####
simple1_dist_matrix = vegdist(filttest_simple1, method="bray",diag=TRUE)
simple1_pcoa = pcoa(simple1_dist_matrix)
simple1_pcoa_df = as.data.frame(simple1_pcoa$vectors[,1:2])
for (col in seq_along(simple_diet1_micro)){
  print(col)
  var = names(simple_diet1_micro[col])
  print(var)
  plot = ggplot(simple1_pcoa_df, aes(x = Axis.1, y = Axis.2, color=factor(simple_diet1_micro[[col]]))) + 
    geom_point(alpha = 5/10, size=3) +
    labs(colour = "g/kg") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(x=paste0("PC1: ",round(simple1_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
         y=paste0("PC2: ",round(simple1_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
    ggtitle(paste("simple diet 1; colored by",var)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = c(.9, .25)) +
    theme(legend.key.size = unit(0.4, "cm")) +
    theme(legend.title = element_text(size = 8), 
          legend.text = element_text(size = 8))
  ggsave(plot, file=paste0("/Users/tobynbranck/Documents/Effect_size/analysis/mouse_ords_overlayed_covariates/simple1_ordination_", var,".png"), width = 14, height = 10, units = "cm")
}

##### overlaying ordinations by covariates for simple diet 2 #####
simple2_dist_matrix = vegdist(filttest_simple2, method="bray",diag=TRUE)
simple2_pcoa = pcoa(simple2_dist_matrix)
simple2_pcoa_df = as.data.frame(simple2_pcoa$vectors[,1:2])
for (col in seq_along(simple_diet2_micro)){
  print(col)
  var = names(simple_diet2_micro[col])
  print(var)
  plot = ggplot(simple2_pcoa_df, aes(x = Axis.1, y = Axis.2, color=factor(simple_diet2_micro[[col]]))) + 
    geom_point(alpha = 5/10, size=3) +
    labs(colour = "g/kg") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(x=paste0("PC1: ",round(simple2_pcoa$values$Relative_eig[1],digits=3)*100,"%"),
         y=paste0("PC2: ",round(simple2_pcoa$values$Relative_eig[2],digits=3)*100,"%")) +
    ggtitle(paste("simple diet 2; colored by",var)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = c(.9, .25)) +
    theme(legend.key.size = unit(0.4, "cm")) +
    theme(legend.title = element_text(size = 8), 
          legend.text = element_text(size = 8))
  ggsave(plot, file=paste0("/Users/tobynbranck/Documents/Effect_size/analysis/mouse_ords_overlayed_covariates/simple2_ordination_", var,".png"), width = 14, height = 10, units = "cm")
}


##### simple 1 analysis #####
adonis_res_rsq_simple1 = vector()
adonis_res_pval_simple1 = vector()
simple_diet1_micro[,5:14] <- lapply(simple_diet1_micro[,5:14], as.numeric)
simple_diet1_micro_adon = simple_diet1_micro[,5:14]
for (col in seq_along(simple_diet1_micro_adon)){
  print (simple_diet1_micro_adon[[col]])
  adonis.univ_simple1 = adonis(filttest_simple1 ~ simple_diet1_micro_adon[[col]], data = simple_diet1_micro_adon, method = "bray", strata = simple_diet1_micro$Individual,na.rm=TRUE)
  print (adonis.univ_simple1)
  adonis_res_rsq_simple1[col] = adonis.univ_simple1$aov.tab[1,]$R2
  adonis_res_pval_simple1[col] = adonis.univ_simple1$aov.tab[1,]$`Pr(>F)`}
univar_res_rgnav_simple1 = rbind(adonis_res_pval_simple1, adonis_res_rsq_simple1)
univar_res_rgnav_simple1 = as.data.frame(t(univar_res_rgnav_simple1))
names(univar_res_rgnav_simple1) = c("P-Value", "R2")
univar_res_rgnav_simple1$`P-Value` = as.numeric(univar_res_rgnav_simple1$`P-Value`)
univar_res_rgnav_simple1$R2 = as.numeric(univar_res_rgnav_simple1$R2)
univar_res_rgnav_simple1$p_adj = p.adjust(univar_res_rgnav_simple1$`P-Value`, "fdr")
univar_res_rgnav_simple1$stars = cut(univar_res_rgnav_simple1$p_adj, c(0, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "-", ""))
univar_res_rgnav_simple1$diet_variable = colnames(simple_diet1_micro_adon)
univar_res_rgnav_simple1 = univar_res_rgnav_simple1[order(univar_res_rgnav_simple1$R2),]

ggplot(univar_res_rgnav_simple1, aes(x=reorder(diet_variable,R2), y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() +
  geom_point(data = univar_res_rgnav_simple1[univar_res_rgnav_simple1$p_adj < 0.05, ],aes(x=reorder(diet_variable,R2), R2 + 0.01), shape = "*", size=4.233, color="black") +
  geom_point(data = univar_res_rgnav_simple1[univar_res_rgnav_simple1$p_adj < 0.1, ],aes(x=reorder(diet_variable,R2), R2 + 0.01), shape = "-", size=4.233, color="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Mouse cohort adonis, simple experiment 1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(expression("R^2")) +
  xlab("Diet variable")
#maaslin2
metadata_input_simple_diet1_micro = simple_diet1_micro
metadata_input_simple_diet1_micro$sampling_day = NULL
metadata_input_simple_diet1_micro$diet_period = NULL
metadata_input_simple_diet1_micro$diet_type = NULL
names(metadata_input_simple_diet1_micro) <- gsub(" ", "_", names(metadata_input_simple_diet1_micro))
names(metadata_input_simple_diet1_micro) <- gsub("-", "_", names(metadata_input_simple_diet1_micro))
names(metadata_input_simple_diet1_micro) <- gsub("_", ".", names(metadata_input_simple_diet1_micro))

colnames(metadata_input_simple_diet1_micro)[7] <- "Mineral_mix"
colnames(metadata_input_simple_diet1_micro)[10] <- "Vitamin_mix"
individuals = metadata_input_simple_diet1_micro[,1]
metadata_input_simple_diet1_micro[,2:11] <- lapply(metadata_input_simple_diet1_micro[,2:11], as.numeric)

metadata_input_simple_diet1_micro %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()
#Maaslin2(input_data = filttest_strict_simple1, input_metadata = metadata_input_simple_diet1_micro, output = "/Users/tobynbranck/Documents/Effect_size/analysis/mouse_simple1_diet_maaslin2", random_effects = c("Individual"))

metadata_input_simple_diet1_micro[] <- lapply(metadata_input_simple_diet1_micro, as.character)
#Maaslin2(input_data = filttest_strict_simple1, input_metadata = metadata_input_simple_diet1_micro, output = "/Users/tobynbranck/Documents/Effect_size/analysis/mouse_simple1_diet_maaslin2_categorical", random_effects = c("Individual"), reference=c("Sucrose","81.38"))

####simple 2 analysis
simple_diet2_micro$`Maltodextrin Lo-Dex 10` = NULL
simple_diet2_micro$`Cellulose (Fiber)` = NULL
simple_diet2_micro$`Ethoxyquin (Liquid)` = NULL

adonis_output_simple2 = vector()
adonis_res_rsq_simple2 = vector()
adonis_res_pval_simple2 = vector()
simple_diet2_micro[,5:14] <- lapply(simple_diet2_micro[,5:14], as.numeric)
#simple_diet2_micro[,5:14] <- lapply(simple_diet2_micro[,5:14], as.character)
simple_diet2_micro_adon = simple_diet2_micro[,5:14]
for (col in seq_along(simple_diet2_micro_adon)){
  print (simple_diet2_micro_adon[[col]])
  adonis.univ_simple2 = adonis(filttest_simple2 ~ simple_diet2_micro_adon[[col]], data = simple_diet2_micro_adon, method = "bray", strata = simple_diet2_micro$Individual,na.rm=TRUE)
  adonis_res_rsq_simple2[col] = adonis.univ_simple2$aov.tab[1,]$R2
  adonis_res_pval_simple2[col] = adonis.univ_simple2$aov.tab[1,]$`Pr(>F)`}
univar_res_rgnav_simple2 = rbind(adonis_res_pval_simple2, adonis_res_rsq_simple2)
univar_res_rgnav_simple2 = as.data.frame(t(univar_res_rgnav_simple2))
names(univar_res_rgnav_simple2) = c("P-Value", "R2")
univar_res_rgnav_simple2$`P-Value` = as.numeric(univar_res_rgnav_simple2$`P-Value`)
univar_res_rgnav_simple2$R2 = as.numeric(univar_res_rgnav_simple2$R2)
univar_res_rgnav_simple2$p_adj = p.adjust(univar_res_rgnav_simple2$`P-Value`, "fdr")
univar_res_rgnav_simple2$stars = cut(univar_res_rgnav_simple2$p_adj, c(0, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "-", ""))
univar_res_rgnav_simple2$diet_variable = colnames(simple_diet2_micro_adon)
univar_res_rgnav_simple2 = univar_res_rgnav_simple2[order(univar_res_rgnav_simple2$R2),]

ggplot(univar_res_rgnav_simple2, aes(x=reorder(diet_variable,R2), y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() +
  geom_point(data = univar_res_rgnav_simple2[univar_res_rgnav_simple2$p_adj < 0.05, ],aes(x=reorder(diet_variable,R2), R2 + 0.01), shape = "*", size=4.233, color="black") +
  geom_point(data = univar_res_rgnav_simple2[univar_res_rgnav_simple2$p_adj < 0.1, ],aes(x=reorder(diet_variable,R2), R2 + 0.01), shape = "-", size=4.233, color="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Mouse cohort adonis, simple experiment 2") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(expression("R^2")) +
  xlab("Diet variable")

#maaslin2
metadata_input_simple_diet2_micro = simple_diet2_micro
metadata_input_simple_diet2_micro$sampling_day = NULL
metadata_input_simple_diet2_micro$diet_period = NULL
metadata_input_simple_diet2_micro$diet_type = NULL
names(metadata_input_simple_diet2_micro) <- gsub(" ", "_", names(metadata_input_simple_diet2_micro))
names(metadata_input_simple_diet2_micro) <- gsub("-", "_", names(metadata_input_simple_diet2_micro))
names(metadata_input_simple_diet2_micro) <- gsub("_", ".", names(metadata_input_simple_diet2_micro))

#colnames(metadata_input_simple_diet2_micro)[2] <- "apple_sauce"
colnames(metadata_input_simple_diet2_micro)[7] <- "Mineral_mix"
colnames(metadata_input_simple_diet2_micro)[10] <- "Vitamin_mix"
metadata_input_simple_diet2_micro[,2:11] <- lapply(metadata_input_simple_diet2_micro[,2:11], as.numeric)
metadata_input_simple_diet2_micro %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()
Maaslin2(input_data = filttest_strict_simple2, input_metadata = metadata_input_simple_diet2_micro, output = "/Users/tobynbranck/Documents/Effect_size/analysis/mouse_simple2_diet_maaslin2", random_effects = c("Individual"))

metadata_input_simple_diet2_micro[] <- lapply(metadata_input_simple_diet2_micro, as.character)
Maaslin2(input_data = filttest_strict_simple2, input_metadata = metadata_input_simple_diet2_micro, output = "/Users/tobynbranck/Documents/Effect_size/analysis/mouse_simple2_diet_maaslin2_categorical", random_effects = c("Individual"))
