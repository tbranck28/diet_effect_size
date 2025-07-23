#This script analyzes the non-human primate microbiome dataset to
#understand the effects of diet on microbiome in a non-human gut environment.
#adapted from the tutorial found here: https://benjjneb.github.io/dada2/tutorial.html
#resulting number of samples after pairing taxonomic profiles with dietary metadata: 34 (13 high fat diet; 21 low fat)
#resulting number of individuals: 34

library(dada2); packageVersion("dada2")
library(phyloseq)
library(Maaslin2)
library(vegan)

#set directory paths
input_dir = file.path("/home","tobynb","diet_test","input","nonhuman_primates","Initial",fsep = "/")
output_dir = file.path("/home","tobynb","diet_test","output_test",fsep = "/")
seqs_path <- "/home/tobynb/diet_test/input/nonhuman_primates/Initial/NHP_seqs"
list.files(seqs_path)

### manually performing QC of reads (meaning outside of biobakery workflow) to get list of samples that have sufficient read quality
#### string manipulation to parse sample names (and to sort forward vs. reverse reads but ***not relavent
#### in this case)
fnFs <- list.files(seqs_path, pattern="SRR", full.names = TRUE)
sample.names = sapply(strsplit(basename(fnFs), "[.]"), `[`, 1)

#### plots quality of reads in first two samples
plotQualityProfile(fnFs[1:2])

#### creates directory for filtered reads, assigns name to files containing post-filtered/trimmed reads,
#### performs filter/trim function with parameters specified below
filtFs <- file.path(seqs_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, truncLen=c(50),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

#### models and plots post-filtered read error rates
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

plotQualityProfile(filtFs[1:3])

#### getting list of samples with low quality reads #########
filt_qual_plot = plotQualityProfile(filtFs[1:74])
filt_scores = ggplot_build(filt_qual_plot)
mean_qual_scors = as.data.frame(filt_scores$data[[2]])
qualscores = as.vector(mean_qual_scors$y)
means = .colMeans(qualscores, 50, length(qualscores) / 50)
indx_keep = which(means >= 32)
filtFs = filtFs[indx_keep]
filt_scores

#### sample inference on filtered reads
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaFs[[1]]
# ****skipping "merge" function since only working with single-end reads

#### constructs ASV table - this will be used for filtering samples only
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

#### removes chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# #### tracks number of reads that make it through each step of pipeline --below two sections are commented out and not needed, since the biobakery employment
#of dada2 handles this. 
# getN <- function(x) sum(getUniques(x))
# track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
# rownames(track) <- sample.names
# head(track)
# 
# ### assign taxonomy function
# taxa <- assignTaxonomy(seqtab.nochim, paste0(input_dir,"/silva_nr_v132_train_set.fa.gz"), tryRC=TRUE, multithread=TRUE)
# taxa.print <- taxa
# rownames(taxa.print) <- NULL
# head(taxa.print)
# unname(taxa)

#bring in metaphlan table generated from biobakery (which employed dada2) and format table
metaphlan_profiles = as.data.frame(data.table::fread(paste0(input_dir,"/all_samples_taxonomy_closed_reference_withseqs.tsv"), header=T, check.names=T))
metaphlan_profiles$taxonomy = NULL
metaphlan_profiles = as.data.frame(t(metaphlan_profiles))
names(metaphlan_profiles) <- metaphlan_profiles[1,]
metaphlan_profiles = metaphlan_profiles[-1,]
rownames(metaphlan_profiles) = gsub(".fastq","",rownames(metaphlan_profiles))
metaphlan_profiles = mutate_all(metaphlan_profiles, function(x) as.numeric(as.character(x)))

#keep samples that are in the seqtab.nochim table, recall we filtered samples to remove those with low quality reads
metaphlan_profiles = metaphlan_profiles[rownames(metaphlan_profiles) %in% rownames(seqtab.nochim),]

phyloseq::otu_table(metaphlan_profiles,taxa_are_rows = FALSE)

#bring in ASV table from biobakery workflow which employed dada2
se = readRDS(paste0(input_dir,"/seqtab_final.rds"))
test = assignTaxonomy(se, paste0(input_dir,"/silva_nr_v132_train_set.fa.gz"), tryRC=TRUE, multithread=TRUE)
taxa.print <- test
rownames(taxa.print) <- NULL
head(taxa.print)
unname(test)

#########################################
######### phyloseq portion ##############
#########################################

#read in metadata file from SRA
samples_df = data.table::fread(paste0(input_dir,"/nonhumanPrimate_metadata.txt"), header=T, check.names=F)

#make "Run" column the rownames
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Run")

#create phyloseq object
ps <- phyloseq(otu_table(metaphlan_profiles, taxa_are_rows=FALSE), sample_data(samples_df), tax_table(test))

#set up container for asv sequences using the DNAString function (for easier string manipulation)
dna <- Biostrings::DNAStringSet(taxa_names(ps))

#assigning name to each seq
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

### separate the df_otutable by diet group, high fat diet vs. low fat diet, so that we can apply filtering strategy by diet group
samples_df_h = samples_df[samples_df$Host_diet == 'high fat',]
samples_df_l = samples_df[samples_df$Host_diet == 'herbivore',]
df_otutable_h = subset(df_otutable, row.names(df_otutable) %in% row.names(samples_df_h))
df_otutable_l = subset(df_otutable, row.names(df_otutable) %in% row.names(samples_df_l))

df_otutable_h = as.data.frame(t(df_otutable_h))
df_otutable_l = as.data.frame(t(df_otutable_l))

filtered_h = df_otutable_h[rowSums(df_otutable_h > 0.001) >= dim(df_otutable_h)[2]*.25, ]
filtered_h = as.data.frame(t(filtered_h))

filtered_l = df_otutable_l[rowSums(df_otutable_l > 0.001) >= dim(df_otutable_l)[2]*.25, ]
filtered_l = as.data.frame(t(filtered_l))

# find union of ASVs after apply prevalence/abundance filtering to the subsetted-by-diet samples; merge filtered_h and filtered_l into filtered
h_asvs = colnames(filtered_h)
l_asvs = colnames(filtered_l)
union_filt = union(h_asvs,l_asvs)
filtered = df_otutable[,union_filt]
filtered = filtered[,mixedsort(colnames(filtered))]

##### subset sample metadata matrix to match QC'd samples in filtered ASV matrix #####
samples_df = samples_df['Host_diet']
samples_df$samples = rownames(samples_df)
filtered = filtered[rownames(filtered) %in% rownames(samples_df),]
samples_df = samples_df[rownames(samples_df) %in% rownames(filtered),]
samples_df$samples = NULL
identical(row.names(filtered), row.names(samples_df))

write.csv(filtered, paste0(output_dir,"/nonHumanPrim_Phenotp_.001abun.25prev_table.csv"), row.names = TRUE)
write.csv(samples_df, paste0(output_dir,"/nonHumanPrim_processedmetadata.csv"), row.names = TRUE)

######## adonis section ##########
identical(row.names(filtered), row.names(samples_df))
adonis_obj = adonis2(filtered ~ samples_df$Host_diet, data = samples_df, method = "bray",by="margin")

######## maaslin2 section ########
Maaslin2(input_data = filtered, input_metadata = samples_df, output = paste0(output_dir,"/nonhumanPrimates_maaslin2_output"),max_significance = 1)

### making plot for NHP bugs effect size vs. qvalue ###
NHP_results = as.data.frame(data.table::fread(paste0(output_dir,"/nonhumanPrimates_maaslin2_output/all_results.tsv"), header=T, check.names=T))
cairo_pdf(paste0(output_dir,"/NHP_maaslin_plot_volcano.pdf"),width=8,height=8)
ggplot(data=NHP_results, aes(x=coef, y=-log10(qval),label=ifelse(NHP_results$coef<2 & NHP_results$coef>-2,"",NHP_results$feature))) + 
  geom_point(size=2,color=ifelse(NHP_results$coef<2 & NHP_results$coef>-2, 'grey', 'black')) +
  theme_classic() +
  xlab("Effect Size") +
  ylab("-log10(q-value)") +
  ggrepel::geom_label_repel(size=3,box.padding = unit(0.2, "lines")) +
  theme(axis.text=element_text(size=10))
dev.off() 