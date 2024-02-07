#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="bbduk"
#SBATCH --output=bbduk.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa
conda activate qiime2-2019.4


adapters=/mnt/lustre/software/anaconda/colsa/envs/qiime2-2019.4/opt/bbmap-38.57-0/resources/adapters.fa

cd ~/EDRR_patho

for i in reads/*R1*
do(
    dir=${i%/*}
    r1File=${i##*/}
    pre=${r1File%R1*}
    post=${r1File##*R1}
    r2File=${pre}R2{post}
    
    mkdir bbtrimmed_$dir
    
    srun bbduk.sh in=$dir/$r1File in2=$dir/$r2File \
        out=bbtrimmed_$dir/$r1File out2=bbtrimmed_$dir/$r2File \
        -ref=$adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo \
        qtrim=r trimq=10 -minlength=36
)
done

