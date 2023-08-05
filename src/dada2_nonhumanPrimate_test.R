setwd("/Users/tobynbranck/Documents/Effect_size/analysis/dada2_test")
library(dada2); packageVersion("dada2")
library(phyloseq)
library(Maaslin2)
library(vegan)

path <- "/Users/tobynbranck/Documents/Effect_size/analysis/dada2_test/nonhumanPrimates"
list.files(path)

#### string manipulation to parse sample names (and to sort forward vs. reverse reads but ***not relevent
#### in my case)
fnFs <- list.files(path, pattern="SRR", full.names = TRUE)
sample.names = sapply(strsplit(basename(fnFs), "[.]"), `[`, 1)

#### plots quality of reads in first two samples
pdf(file = "nonHumanPrimReads_quality_profile.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4)
plotQualityProfile(fnFs[1:2])
dev.off()

#### creates directory for filtered reads, assigns name to files containing post-filtered/trimmed reads,
#### performs filter/trim function with parameters specified below
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, truncLen=c(50),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

#### models and plots post-filtered read error rates
errF <- learnErrors(filtFs, multithread=TRUE)
pdf(file = "nonHumanPrim_read_errors.pdf",
         width = 4,
         height = 4)
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf(file = "nonHumanPrimReads_quality_profile_PostFilter73_74.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)
plotQualityProfile(filtFs[1:3])
dev.off()

#### filtering out samples with low quality reads #########
filt_qual_plot = plotQualityProfile(filtFs[1:74])
filt_scores = ggplot_build(filt_qual_plot)
mean_qual_scors = as.data.frame(filt_scores$data[[2]])
qualscores = as.vector(mean_qual_scors$y)
means = .colMeans(qualscores, 50, length(qualscores) / 50)
indx_keep = which(means >= 32)
filtFs = filtFs[indx_keep]

cairo_pdf("NHP_quality_scores.pdf",width=12,height=12)
filt_scores
dev.off()

#### sample inference on filtered reads
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaFs[[1]]
# ****skipping "merge" function since only working with single-end reads

#### constructs ASV table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#### removes chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#### tracks number of reads that make it through each step of pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

### assign taxonomy function
taxa <- assignTaxonomy(seqtab.nochim, "./silva_nr_v132_train_set.fa.gz", tryRC=TRUE, multithread=TRUE)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
unname(taxa)

###from metaphlan ###
metaphlan_profiles = as.data.frame(data.table::fread("/Users/tobynbranck/Downloads/all_samples_taxonomy_closed_reference_withseqs.tsv", header=T, check.names=T))
metaphlan_profiles$taxonomy = NULL
metaphlan_profiles = as.data.frame(t(metaphlan_profiles))
names(metaphlan_profiles) <- metaphlan_profiles[1,]
metaphlan_profiles = metaphlan_profiles[-1,]
rownames(metaphlan_profiles) = gsub(".fastq","",rownames(metaphlan_profiles))
metaphlan_profiles = mutate_all(metaphlan_profiles, function(x) as.numeric(as.character(x)))

metaphlan_profiles = metaphlan_profiles[rownames(metaphlan_profiles) %in% rownames(seqtab.nochim),]

phyloseq::otu_table(metaphlan_profiles,taxa_are_rows = FALSE)

se = readRDS("/Users/tobynbranck/Downloads/seqtab_final.rds")
test = assignTaxonomy(se, "./silva_nr_v132_train_set.fa.gz", tryRC=TRUE, multithread=TRUE)
taxa.print <- test
rownames(taxa.print) <- NULL
head(taxa.print)
unname(test)


#########################################
######### phyloseq portion ##############
#########################################

#read in metadata file
samples_df = data.table::fread("/Users/tobynbranck/Documents/Effect_size/data/nonhumanPrimate_metadata.txt", header=T, check.names=F)
#make "Run" column the into the rownames
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Run")

#create phyloseq object
ps <- phyloseq(otu_table(metaphlan_profiles, taxa_are_rows=FALSE), sample_data(samples_df), tax_table(test))

#create phyloseq object
#ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samples_df), tax_table(taxa))

#set up container for asv sequences using the DNAString function (for easier string manipulation) ?
dna <- Biostrings::DNAStringSet(taxa_names(ps))
#assigning name to each seq - not sure why we need to assign each seq's name to be the actual sequence
names(dna) <- taxa_names(ps)
#merging the phyloseq object with the dna seqs object
ps <- merge_phyloseq(ps, dna)
#assigning taxa name as ASVXXX to each sequence
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
#transforms otus into relative abundance
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
### make dataframe out of otu_table(ps) 
ps.prop = tax_glom(ps.prop, taxrank="Genus", NArm=TRUE)
df_otutable = as.data.frame(otu_table(ps.prop))
df_taxa = as.data.frame(tax_table(ps.prop))
df_taxa = df_taxa["Genus"]
names(df_otutable) = df_taxa$Genus
##### separate the df_otutable by phenotype, high fat diet vs. low fat diet, so that we can
##### apply filtering strategy by phenotype
#split df_otutable:
samples_df_h = samples_df[samples_df$Host_diet == 'high fat',]
samples_df_l = samples_df[samples_df$Host_diet == 'herbivore',]
df_otutable_h = subset(df_otutable, row.names(df_otutable) %in% row.names(samples_df_h))
df_otutable_l = subset(df_otutable, row.names(df_otutable) %in% row.names(samples_df_l))

df_otutable_h = as.data.frame(t(df_otutable_h))
df_otutable_l = as.data.frame(t(df_otutable_l))

filtered_h = df_otutable_h[rowSums(df_otutable_h > 0.001) >= dim(df_otutable_h)[2]*.25, ]
filtered_h = data.frame(filtered_h)
filtered_h = as.data.frame(t(filtered_h))

filtered_l = df_otutable_l[rowSums(df_otutable_l > 0.001) >= dim(df_otutable_l)[2]*.25, ]
filtered_l = data.frame(filtered_l)
filtered_l = as.data.frame(t(filtered_l))

filtered_strict_h = df_otutable_h[rowSums(df_otutable_h > 0.001) >= dim(df_otutable_h)[2]*.25, ]
filtered_strict_h = data.frame(filtered_strict_h)
filtered_strict_h = as.data.frame(t(filtered_strict_h))

filtered_strict_l = df_otutable_l[rowSums(df_otutable_l > 0.001) >= dim(df_otutable_l)[2]*.25, ]
filtered_strict_l = data.frame(filtered_strict_l)
filtered_strict_l = as.data.frame(t(filtered_strict_l))

# find union of ASVs after apply prevalence/abundance filtering to the subsetted-by-phenotype samples -
# two groups - high fat and herbivore; merge filtered_h and filtered_l into filtered
h_asvs = colnames(filtered_h)
l_asvs = colnames(filtered_l)
union_filt = union(h_asvs,l_asvs)
filtered = df_otutable[,union_filt]
filtered = filtered[,mixedsort(colnames(filtered))]
# same for the strict filtering; merge filtered_strict_h and filtered_strict_l into filtered_strict
h_asvs_strict = colnames(filtered_strict_h)
l_asvs_strict = colnames(filtered_strict_l)
union_filt_strict = union(h_asvs_strict,l_asvs_strict)
filtered_strict = df_otutable[,union_filt_strict]
filtered_strict = filtered_strict[,mixedsort(colnames(filtered_strict))]

##### subset sample metadata matrix to match QC'd samples in filtered ASV matrix #####
samples_df = samples_df['Host_diet']
filtered = filtered[rownames(filtered) %in% rownames(samples_df),]
samples_df = samples_df[rownames(samples_df) %in% rownames(filtered),]
samples_df = as.data.frame(samples_df)
samples_df$Host_diet = samples_df$samples_df
rownames(samples_df) = rownames(filtered)
samples_df$samples_df = NULL

identical(row.names(filtered), row.names(samples_df))

write.csv(filtered, "nonHumanPrim_Phenotp_.001abun.25prev_table.csv", row.names = TRUE)
write.csv(samples_df, "nonHumanPrim_processedmetadata.csv", row.names = TRUE)


######## adonis section ##########
identical(row.names(filtered), row.names(samples_df))
adonis_obj = adonis(filtered ~ samples_df$Host_diet, data = samples_df, method = "bray")

######## maaslin2 section, need more stringent prev filtering ########
samples_df = samples_df['Host_diet']
#filtered_strict = filtered_strict[rownames(filtered_strict) %in% rownames(samples_df),]
samples_df = samples_df[rownames(samples_df) %in% rownames(filtered),]
samples_df = as.data.frame(samples_df)
samples_df$Host_diet = samples_df$samples_df
rownames(samples_df) = rownames(filtered)
samples_df$samples_df = NULL

Maaslin2(input_data = filtered, input_metadata = samples_df, output = "nonhumanPrimates_maaslin2_output",max_significance = 1)

### making plot for NHP bugs effect size vs. qvalue ###
NHP_results = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/Effect_size/analysis/dada2_test/nonhumanPrimates_maaslin2_output/all_results.tsv", header=T, check.names=T))
ggplot(NHP_results,aes(x=qval,y=coef)) +
  geom_point()

library(ggrepel)
cairo_pdf("NHP_maaslin_plot.pdf",width=10,height=8)
ggplot(NHP_results,aes(x=qval,y=coef,label=ifelse(NHP_results$coef<2 & NHP_results$coef>-2,"",NHP_results$feature))) +
    geom_point(size=2,color=ifelse(NHP_results$coef<2 & NHP_results$coef>-2, 'grey', 'black')) +
    xlab("q-value") +
    ylab("Effect size (b-coefficient)") +
    geom_vline(xintercept=0.25, linetype="dashed", color = "grey") +
    geom_hline(yintercept=2, linetype="dashed", color = "grey") +
    geom_hline(yintercept=-2, linetype="dashed", color = "grey") +
    theme_classic() +
    scale_y_continuous(breaks=c(-4,-3,-2,-1,0,1,2,3,4)) +
    geom_label_repel(size=3,box.padding = unit(0.2, "lines"))
dev.off()

cairo_pdf("NHP_maaslin_plot_volcano.pdf",width=8,height=8)
ggplot(data=NHP_results, aes(x=coef, y=-log10(qval),label=ifelse(NHP_results$coef<2 & NHP_results$coef>-2,"",NHP_results$feature))) + 
  geom_point(size=2,color=ifelse(NHP_results$coef<2 & NHP_results$coef>-2, 'grey', 'black')) +
  theme_classic() +
  xlab("Effect Size") +
  ylab("-log10(q-value)") +
  geom_label_repel(size=3,box.padding = unit(0.2, "lines")) +
  theme(axis.text=element_text(size=10))
dev.off() 
  
###### Extra code #############################################
###### saving taxa table (taxonomic annotations of ASVs) ######
write.csv(as.data.frame(tax_table(ps.prop)), "taxa_table_ps.prop.csv", row.names = TRUE)
taxonomy_df = as.data.frame(tax_table(ps.prop))
filtered_strict_t = as.data.frame(t(filtered_strict))
ASVs_postfilter = taxonomy_df[rownames(taxonomy_df) %in% rownames(filtered_strict_t),]
write.csv(ASVs_postfilter, "taxa_after_filter_by_pt.csv", row.names = TRUE)


##############     Extra phyloseq functions     ###############
#makes ordination object
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
#plots ordination
ordination = plot_ordination(ps.prop, ord.nmds.bray, color="Host_diet", title="Bray NMDS")

##subsetting only the Genus ASVs
ps_genus <- tax_glom(ps.prop, "Genus")
#annotating otu_table with genus names (replacing ASVXXX with genus name)
taxa_names(ps_genus) <- tax_table(ps_genus)[, 6]
#transforming otu_table into dataframe and saving as csv
otu_table_genus = as.data.frame(otu_table(ps_genus))
write.csv(otu_table_genus, "nonHumanPrim_abundance_table.csv", row.names = TRUE)

otu_table_genus = t(otu_table_genus)
otu_table_genus = as.data.frame(otu_table_genus)
filtered = otu_table_genus[rowSums(otu_table_genus > 0.0001) >= dim(otu_table_genus)[2]*.25, ]
filtered = data.frame(filtered)

################### extra code for visualization ###########################################
#for taxonomic bar plot: selects top 20 most abundant ASVs, re-calculates relative abundance
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:25]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
barplot = plot_bar(ps.top20, fill="Family") + facet_wrap(~Host_diet, scales="free_x")
