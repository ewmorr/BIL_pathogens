
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
