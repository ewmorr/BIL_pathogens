

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
```
The above is not working well. Will have to retry. Old version of itsxpress mod is installed at .local/bin
```
conda create --name trim_3p --clone template #on premise clone the template
#conda create -n trim_3p pip hmmer bbmap vsearch biopython
conda activate trim_3p
conda install pip hmmer bbmap vsearch=2.19 biopython=1.60 python=3.6

git clone https://github.com/ewmorr/itsxpress.git
cd ~/repo/itsxpress
git branch
git branch 3p_trim origin/3p_trim
git checkout 3p_trim
pip install --user -e .

