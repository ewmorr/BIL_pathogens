#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="dada2"
#SBATCH --output=dada2_workflow.itsxpress.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24

module purge
module load anaconda/colsa
conda activate itsxpress_2x

cd /mnt/home/garnas/ewj4/EDRR_patho

mkdir itsx_test

for i in itsxpress_out/*R1*
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
