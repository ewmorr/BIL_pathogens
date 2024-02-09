#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="dada2"
#SBATCH --output=dada2_workflow.primer_qual_check.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa
conda activate conda activate r-dada2_new

cd /mnt/home/garnas/ewj4/EDRR_patho

#Run first R scripts
Rscript ~/repo/BIL_pathogens/sequence_processing_and_dada2/1_dada2_primer_check.r

conda deactivate

module purge
module load linuxbrew/colsa

mkdir /mnt/home/garnas/ewj4/EDRR_patho/cutadapt

#fwd
fwd=AGCCTCCGCTTATTGATATGCTTAART
#fwd.rc
fwd.rc=AYTTAAGCATATCAATAAGCGGAGGCT
#rev
rev=AACTTTYRRCAAYGGATCWCT
#rev.rc
rev.rc=AGWGATCCRTTGYYRAAAGTT

for i in filtN/*
do(

    dir=${i%/*}
    r1File=${i##*/}
    pre=${r1File%R1*}
    post=${r1File##*R1}
    r2File=${pre}R2${post}

    cutadapt -g $fwd -a $rev.rc \
        -G $rev -A $fwd.rc \
        -m 50 -n 2 \
        -o cutadapt/$r1File -p cutadapt/$r2File \
        $r1File $r2File
)
done

module purge
module load anaconda/colsa
conda activate r-dada2_new

Rscript ~/repo/BIL_pathogens/sequence_processing_and_dada2/2_post_cutadapt_primer_check.r
#this is a good place to stop to check primer orientation and quality profiles
