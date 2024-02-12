

### Initial data testing effects of sample volume across 8 samples (5, 10, 20,40, 80 ml) was downloaded and saved to `EDRR/vol_test` on premise
### First run dada2 workflow for primer checks and primer removal
```
dada2_workflow.primer_and_qual_checks.slurm
```
## ITSxpress v>2 **still** does not trim 3' ends of sequences
### We set up to conda envs `itsxpress_2x` and `itsxpress_2x_mod` to test
### The modification is at line 89 of Dedup.py
```
cd repo
git clone https://github.com/USDA-ARS-GBRU/itsxpress.git itsxpress_2x_mod
conda create --name itsxpress_2x_mod --clone template
conda activate itsxpress_2x_mod
conda install conda-build
conda develop itsxpress_2x_mod
