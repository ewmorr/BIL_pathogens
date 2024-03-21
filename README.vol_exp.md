
Reads download on premise in EDRR_patho dir
```
cd EDRR_patho
mv cobb.sr.unh.edu/managed/240126_A01346_0126_AHJKWVDRX3_16Mer012624-EM-EDRR-patho/reads ./reads
```
initial qual filter test with bbduk
```
sbatch ~/repo/BIL_pathogens/premise/bbduk.slurm
```
read numbers look good
### Beginnning dada2 workflow by checking primer orientation and read quality
```
sbatch ~/repo/BIL_pathogens/premise/dada2_workflow.primer_and_qual_checks.slurm
```

DADA2 error model issue NovaSeq. Ben Callahan doesn't think it's a problem
https://forum.qiime2.org/t/novaseq-and-dada2-incompatibility/25865/12

## Initial processing in R. Checking seq counts and running rarefaction/calculation of diversity
Run interactively. Dir is `R_scripts/vol_test`
```
vol_test.exploratory_seq_counts.R
run_and_save_rarefactions.R
```
Then run div calcs (run on command line from top level of repo)
```
Rscript R_scripts/vol_test/avg_dist_and_div.R data/vol_test/rarefactions.vol.rds data/vol_test/vol. &
Rscript R_scripts/vol_test/avg_dist_and_div.R data/vol_test/rarefactions.vol.no_min.rds data/vol_test/vol.no_min. &
Rscript R_scripts/vol_test/avg_dist_and_div.R data/vol_test/rarefactions.wash.rds data/vol_test/wash. &
```
Exploratory analysis of richness and beta-div
```
richness_stats_plots.R
vol_nmds.R
wash_nmds.R
```

## Taxonomic assignment and analysis
### We ran the usual RDP-NBC classifier implemented in dada2, but we would like to be able to report a %identity as a measure of confidence of the tax assignment (the NBC bootstrap values are difficult to understand and may not be super meaningful e.g., https://www.drive5.com/usearch/manual/cvi.html)
### BLAST would be fine, but we may implement vsearch --usearch_global functions as a stand-on for blast by performing global seq similarity assignment. For global (and not local search as implemented in BLAST) we will ideally trim the reference database to include the same sites as the query (i.e., positional homology). We therefore run UNITE through ITSx to produce ITS2 seqs.

