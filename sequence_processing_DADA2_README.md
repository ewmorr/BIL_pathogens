

### Initial data testing effects of sample volume across 8 samples (5, 10, 20,40, 80 ml) was downloaded and saved to `EDRR/vol_test` on premise
### First run dada2 workflow for primer checks and primer removal
```
dada2_workflow.primer_and_qual_checks.slurm
```
## ITSxpress v>2 **still** does not trim 3' ends of sequences
### We set up to conda envs `itsxpress_2x` and `itsxpress_2x_mod` to test
### The modification is at line 89 of Dedup.py
https://github.com/conda/conda-build/issues/4251
```
conda create --name itsxpress_2x_mod --clone template
conda activate itsxpress_2x_mod
mamba install -c bioconda itsxpress --only-deps

cd ~/repo
git clone https://github.com/USDA-ARS-GBRU/itsxpress.git itsxpress_2x_mod
cd itsxpress_2x_mod
git branch 3p_trim
git checkout 3p_trim
#make the desired mods to Dedup.py

pip install --no-build-isolation --no-deps -e .
#Successfully installed itsxpress-2.0.1
```
 Old version of itsxpress mod is installed at .local/bin and 3p_trim env
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
```
The local (my laptop) version (itsxpress 1.8) of the 3' trim mod is installed locally and working under itsxpress_3p conda env. The executable is within the conda env. Can do local install of the newer version and mod
``` 
cd repo
git clone https://github.com/USDA-ARS-GBRU/itsxpress.git itsxpress_2x_mod
cd itsxpress_2x_mod
git branch 3p_trim
git checkout 3p_trim
#make the desired mods to Dedup.py

conda env list
conda create -n itsxpress_2x_mod 
conda activate itsxpress_2x_mod
conda install itsxpress --only-deps

cd ~/repo/itsxpress_2x_mod/
pip install --no-build-isolation --no-deps -e .
#Successfully installed itsxpress-2.0.1
```
Now running itsxpress mod on cutadapt processed files
```
cd ~/repo/BIL_pathogens/data/vol_test
mkdir itsxpress_mod_out

#if necessary load conda env
conda activate itsxpress_2x_mod

for i in cutadapt/*R1*
do(

    dir=${i%/*}
    r1File=${i##*/}
    pre=${r1File%R1*}
    post=${r1File##*R1}
    r2File=${pre}R2${post}

    itsxpress --fastq $dir/$r1File --fastq2 $dir/$r2File \
        --reversed_primers --threads 10 --cluster_id 1 --taxa Fungi --region ITS2 \
        --outfile itsxpress_mod_out/$r1File --outfile2 itsxpress_mod_out/$r2File
)
done
```
merging and derep with vsearch. Note that without --allow-merge-stagger about 80-90% of reads fail merging due to staggered merge. This is true for both the mod and unmod version of the itsxpress pipe. The mod should probably be -1 `[start:start+tlen-1]` to avoid this. try this if can get itsxpress_2x_mod set up on server
### After adding the -1 the --alow-merge-stager flag is no longer needed. I.e., reads pass at a normal rate, there are no staggered. This in itself points to a problem in the vanilla program
```
for i in itsxpress_mod_out/*R1*
do(

    dir=${i%/*}
    r1File=${i##*/}
    pre=${r1File%R1*}
    post=${r1File##*R1}
    r2File=${pre}R2${post}

    vsearch --fastq_mergepairs $dir/$r1File --reverse $dir/$r2File \
        --fastqout itsx_test/$pre.merged.fastq.gz \
        --fastq_allowmergestagger
    vsearch --fastx_uniques itsx_test/$pre.merged.fastq.gz --fastaout itsx_test/$pre.unique.fasta
)
done

