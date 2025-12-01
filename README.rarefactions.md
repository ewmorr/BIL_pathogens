# Basic analyses of metabarcoding data

This repo contains R scripts for basic/initial analysis of microbial metabarcoding data (e.g., fungal ITS). The input data is an ASV table of samples x ASVs, such as the ASVs_counts.tsv table output [here](https://github.com/ewmorr/fungi_ITS2_dada2_protocol_07162024) (note that typical bioinformatic pipelines by convention output ASVs/OTUs in rows and samples in columns, which is transposed from the standard layout in most other statistical software). 

The first steps in this protocol are to: 
1. perform rarefaction, i.e., subsample the sequence count table multiple times (e.g., 1000x) saving the subsampled tables in a list and saving as an RDS for future use, and 
2. calculate alpha- and beta-diversity metrics across each of the subsampled tables, average these metrics across all of the subsamples, and then output the average diversity metrics and an average counts table. <br>

The goal of this subsampling procedure is to control for the effects of unequal (highly variable) sequencing depth which is known to bias diversity metrics ([Schloss 2023](https://journals.asm.org/doi/10.1128/msphere.00355-23), [Schloss 2024](https://journals.asm.org/doi/10.1128/msphere.00354-23)). These steps are accomplished by running the scripts [run_and_save_rarefaction.R](./ITS_analysis/run_and_save_rarefaction.R) and [avg_dist_and_div.R](./ITS_analysis/avg_dist_and_div.R), respectively. Note that basic versions of these scripts are in the [`library`](./library) folder, while versions of these scripts that have been modified to process an experiment describing soil fungal community composition at the HBEF site are in the folder [`ITS_analysis`](./ITS_analysis). <br>

The scripts in this repo are constructed assuming that the data is hosted in the same folder as your R project (`.Rproj` file) within a subfolder called `data`, and the HBEF scripts in the `ITS_analysis` directory further draw data from a `data/ITS` directory. See the filepaths in the scripts for guidance on how to setup your directory structure, or modify the file paths accordingly to reflect your preferred directory structure.

After calculation of average diversity metrics across the subsampled data standard ecological analyses can be performed, for example ordination to examine community composition [nmds.R](./ITS_analysis/nmds.R) and ANOVA on alpha-diversity [diversity_plots.R](./ITS_analysis/diversity_plots.R)
