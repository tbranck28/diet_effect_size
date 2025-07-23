Scripts in this directory include all processing and analysis of host taxonomic profiles and associated metadata. Input files README located in data directory. 

ADULT HUMAN SCRIPTS (HMP2, MLVS, HADZA):
hmp2_processing_file.ipynb: takes in clinical and dietary metadata, recodes to be consistent with Long's group (categorical dietary responses converted to ordinal data); creates dietary + other metadata files that are specific to disease phenotypes, as well as one dietary + other metadata file for all phenotypes. Taxonomic profiles are read in and filtered using a standard prevalence and abundance threshold.
	inputs: HMP2_MGX_metadata_09282018.csv; hmp2_metadata.csv; hmp2_species.txt
	outputs: df_CD.csv; df_UC.csv; df_nonIBD.csv; hmp2_meta_allPhenotypes.csv; hmp_species_df.csv

mlvsbb3v_preprocessing.R: takes in metaphlan taxonomic profiles which are filtered using a standard prevalence and abundance threshold; takes in 7ddr energy adjusted-diet variables (averaged over the 2 weeks 6 months apart).
	inputs: metaphlan_taxonomic_profiles.tsv; 7ddrEA_wkAvg_CovsInv.txt
	outputs: hmp_species_df.csv; mlvs_diet.csv

mantel_test_HMP2.R: bins samples into week bins (averages samples) and runs mantel test; mantel tests are stratified by disease phenotypes
	input: hmp_species_df.csv; df_CD.csv; df_UC.csv; df_nonIBD.csv
	output: HMP2_mantel_results.csv

mantel_test_mlvs.R: bins samples into week bins (averages samples) and runs mantel test
	input: mlvs_species.csv; mlvs_diet.csv
	output: MLVS_mantel_results.csv

hmp2_adonis_PCA_disease.R: calculates PCs across all samples (all phenotypes) using psych::principal, uses diet patterns (PCs) and runs adonis for disease-stratified samples
	input: hmp_species_df.csv; df_CD.csv; df_UC.csv; df_nonIBD.csv
	output: HMP2_adonisPCs_byDisease_results.csv

mlvs_adonis_PCs.R: calculates PCs across all samples and runs adonis
	input: mlvs_diet.csv, mlvs_species.csv
	output: mlvs_adonisPCs_results.csv

maaslin2_hmp2_PCs_diseaseStrat.R: calculates and uses dietary patterns (PCs - from psych::principal R package) and runs through maaslin2 for disease-stratified samples. Also runs raw variables (to be used in supplemental figure).
	input: hmp_species_df.csv; hmp2_meta_allPhenotypes.csv
	output: maaslin2 output directories (hmp2_maaslin2_PCs_diseaseStrat_CD, hmp2_maaslin2_PCs_diseaseStrat_UC, hmp2_maaslin2_PCs_diseaseStrat_a_nonIBD, hmp2_maaslin2_raw)

maaslin2_mlvs_PCs.R: calculates and uses dietary patterns (PCs - from psych::principal R package) and runs through maaslin2. Also runs raw variables (to be used in supplemental figure).
	input: mlvs_species.csv, mlvs_diet.csv, mlvs_diet_vars_list.csv (list of mlvs diet variables)
	output: maaslin2 output directories (mlvs_maaslin2_PCs,mlvs_maaslin2_raw)

Hadza_seasonal_preprocessing.R: processing and analysis (adonis, maaslin2) of Hadza seasonal samples.
	input: all_samples_taxonomy_closed_reference.tsv; hadza_seasonal_metadata.txt; 11358_prep_3753_20210706-110244.txt
	output: hadza_seasonal_species.csv; hadza_seasonal_metadata_processed.csv; maaslin2 directory (hadza_seasonal_maaslin2_output_binary)

Figure2_PCs.R: combines all components from mantel, adonis, and maaslin2 scripts (outputs from HMP2, MLVS and Hadza scripts) into main Figure 2, in addition to supplemental figures (XYZ)
	Fig. 2A MANTEL input: HMP2_mantel_results.csv, MLVS_mantel_results.csv
	MANTEL output: Combined_plot_HMP2_MLVS_stratified.png
	SUPPLEMENTAL adonis plot (individual components) input: hmp2CD_adonis_results.csv, hmp2UC_adonis_results.csv, hmp2nonIBD_adonis_results.csv, mlvs_adonis_results.csv
	SUPPLEMENTAL adonis plot (individual components) output: adonis_plot_IndividualComponents.png
	Fig. 2B ADONIS input: mlvs_adonisPCs_results.csv, HMP2_adonisPCs_byDisease_results.csv
	Fig. 2B ADONIS output: Combined_plot_adonis.png
	SUPPLEMENTAL maaslin2 heatmap input: /hmp2_maaslin2_PCs_diseaseStrat_CD/all_results.tsv, /hmp2_maaslin2_PCs_diseaseStrat_UC/all_results.tsv, /hmp2_maaslin2_PCs_diseaseStrat_a_nonIBD/all_results.tsv, /mlvs_maaslin2_PCs/all_results.tsv 
	SUPPLEMENTAL maaslin2 heatmap output: Combined_maaslin2heatmap.pdf
	Fig. 2D specific species interaction with diet (A. hadrus) input: HMP2_PCA_df.csv
	Fig. 2D output: need to save the individual scatters XXXXXX
	Fig. 2F Distribution of beta coefficients input: same inputs as for supplemental heatmap plus /hadza_seasonal_maaslin2_output_binary/all_results.tsv, /nonhumanPrimates_maaslin2_output_filt_by_pt60prev/all_results.tsv, /mouse_complexdiet_maaslin2_binned/all_results.tsv
	Fig. 2F Distribution of beta coefficients output: XXXXXX
	SUPPLEMENTAL Distribution of betas of PCs vs. raw diet components input: same as SUPPLEMENTAL maaslin2 heatmap input plus /hmp2_maaslin2_raw/all_results.tsv and /mlvs_maaslin2_raw/all_results.tsv
	SUPPLEMENTAL Distribution of betas of PCs vs. raw diet components ouput: WesternAdult_PCVsRawBox.png
	MAIN FIGURE OUT: figure2_PCs.pdf
Western_adult_maasline_discretePCs.R: generates SUPPLEMENTAL figure 4 (range of effects when diet is measure as discretized dietary patterns)

INFANT SCRIPTS (BACKHED et al. and MURPHY et al.):
backhed_analysis.R: processing and analysis (adonis and maaslin2 plus supplemental heatmap and line plots) of backhed cohort samples
	inputs: /infant_backhed/metaphlan_taxonomic_profiles.tsv, /infant_backhed/mmc2.csv, /infant_backhed/studyID_person_association.csv, filereport_read_run_PRJEB6456_tsv.txt
	output: backhed_species_for_joint_ord.csv, backhed_meta_joint_ord.csv, backhed_adonis_results.csv, maaslin2 directories (backhed_maaslin2_baseline_8.31, backhed_maaslin2_four_8.31, backhed_maaslin2_twelve_8.31), BackhedBaseline_bugs.csv, BackhedFour_bugs.csv, BackhedTwelve_bugs.csv, BackhedBaseline_meta.csv, BackhedFour_meta.csv, BackhedTwelve_meta.csv
	SUPPLEMENTAL figure heatmap for taxonomic profiles: backhed_heatmap_vertical.pdf
	SUPPLEMENTAL figure bug abundance overtime lineplots: backhed_abund_over_time.pdf

murphy_analysis.R: processing and analysis (adonis and maaslin2 plus supplemental heatmap) of murphy cohort samples.
	inputs: /Infants_murphy/ids.txt, /Infants_murphy/metadata.csv, /Infants_murphy/metaphlan_taxonomic_profiles.tsv, /Infants_murphy/SraRunTable.txt
	output: murphy_species.csv, murphy_species_for_joint_ord.csv, murphy_meta_joint_ord.csv, maaslin2 directories (murphy_maaslin2_baseline_8.31, murphy_maaslin2_3_month_8.31, murphy_maaslin2_twelve_months_8.31), MurphyBaseline_bugs.csv, MurphyThreeMonth_bugs.csv, MurphyTwelveMonth_bugs.csv, MurphyBaseline_meta.csv, MurphyThreeMonth_meta.csv, MurphyTwelveMonth_meta.csv
	SUPPLEMENTAL figure heatmap for taxonomic profiles: murphy_heatmap.pdf

Figure3_infants.R: combines components from Backhed and Murphy individual scripts into Figure 3.
	input: backhed_adonis_results.csv, infants_murphy_adonis_results.csv, all "all_results" from maaslin2 output directories for each age (three time points) for both Backhed and Murphy (6 "all_results" files total) 
	adonis + feeding description table: fig4_adonis_table.csv 
	Fig. 3A output: infants_betas.pdf
	Fig. 3B maaslin2 heatmap (truncated so that only rows/bugs with at least one significant hit are shown): infants_maaslinheatmap_subset.pdf
	Fig. 3C output: forest plot
	SUPPLEMENTAL figure maaslin2 heatmap all species: infants_maaslinheatmap_supplement.pdf
	SUPPLEMENTAL figure spearman for effect sizes measured in both datasets: corr_plots.pdf
	SUPPLEMENRAL figure line plots
	***components put together in adobe illustrator

ANIMAL SCRIPTS (NON-HUMAN PRIMATES, MICE):
dada2_nonhumanPrimate_test.R: processing and analysis (adonis, maaslin2) of NHP samples. Plots quality of reads which helped us filter out samples with low quality reads
	input: nonhumanPrimate_metadata.txt, all_samples_taxonomy_closed_reference_withseqs.tsv
	output: nonHumanPrim_Phenotp_.001abun.25prev_table.csv, nonHumanPrim_processedmetadata.csv, maaslin2 output directory (nonhumanPrimates_maaslin2_output)

mouse_preprocessing.R: processing and analysis of mouse samples with complex (baby food) and simple (refined) diets. Generates main Figure 4.
	input: mouse_faith/metaphlan_taxonomic_profiles.tsv, mouse_simplediet_exp1.csv, mouse_simplediet_exp2.csv, mouse_complexdiet.csv, mouse_simple_diets.txt, mouse_complex_diet.txt
	output: mouse_species.csv, mouse_complex_metadata.csv, mouse_complex_species.csv, maaslin2 directories (mouse_complexdiet_maaslin2_binned,mouse_simple_diets_combined_character_binned)
	Fig. 4A complex diet adonis
	Fig. 4B beta coefficient distributions
	SUPPLEMENTAL figure sample summary ordination: mouse_sampsummary_ord.pdf
	SUPPLEMENTAL figure complex diet heatmap
	SUPPLEMENTAL figure simple diet heatmap
	Fig. 4C E. rectale plots
	SUPPLEMENTAL figure simple combined diets ordination: simpleExperCombined.pdf
	SUPPLEMENTAL figure simple refined diets (adonis and beta distributions): mouse_refined_combined_supplement.png
	Figure 4 main: Figure_4.pdf

OVERVIEW FIGURE (FIGURE ABSTRACT, STUDY CHARACTERISTICS, JOINT AND INDIVIDUAL STUDY ORDINATIONS):
Figure_1.R:
	input: hmp_species_df.csv, mlvs_species.csv, mouse_species.csv, hadza_species.csv, mouse_complex_species.csv, mouse_complex_metadata.csv, backhed_species_for_joint_ord.csv, backhed_meta_joint_ord.csv, murphy_species_for_joint_ord.csv, murphy_meta_joint_ord.csv, nonHumanPrim_Phenotp_.001abun.25prev_table.csv, hadza_seasonal_species.csv, hadza_seasonal_metadata_processed.csv, overview6_21.png, hmp2_meta_allPhenotypes.csv, nonhumanPrimate_metadata.txt
	output: Figure_1_Jun21.pdf
	***need to save and add components to dropbox


