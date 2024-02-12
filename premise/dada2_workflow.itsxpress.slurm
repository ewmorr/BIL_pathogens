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

mkdir itsxpress_out

for i in cutadapt/*R1*
do(

    dir=${i%/*}
    r1File=${i##*/}
    pre=${r1File%R1*}
    post=${r1File##*R1}
    r2File=${pre}R2${post}

    itsxpress --fastq $dir/$r1File --fastq2 $dir/$r2File \
        --reversed_primers --threads 24 --cluster_id 1 --taxa Fungi \
        --outfile itsxpress_out/$r1File --outfile2 itsxpress_out/$r2File
)
done

