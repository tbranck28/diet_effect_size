##### [Replace the contents this README.MD file with the appropriate "Users manual" needed for the workflow.]

# Code and material repository for "Quantifying the effects of diet on microbial communities across different gut contexts" 


#### Analysis pipeline
1. ##### biobakery workflows v3.0.0
   - ##### KneadData v0.7.7
   - ##### MetaPhlAn 3 v3.0
   - ##### HUMAnN 3 v3.0.0
2. ##### Basic statistics using various R package's listed throughout the manuscript and methods
3. ##### Custom R scripts (src directory)


#### Directory Structure

- The `/input` directory contains the by-sample taxonomic profiles, metadata, and intermediate files generated for each study used in the analysis. This directory is organized by study and each study's directory contains `/Initial` and `/Intermediate` sub-directories.

- The src directory is the main source directory for all the custom scripts used for manuscript. The README in this directory describes the code and corresponding input and output files.
