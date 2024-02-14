#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="hmmsearch"
#SBATCH --output=dada2_workflow.itsxpress.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24

module purge
module load anaconda/colsa

conda activate itsxpress_2x

cd /mnt/home/garnas/ewj4/EDRR_patho/itsx_test

for i in EPV_2_*unique.fasta
do(

    sample=${i%.unique.fasta}
        #--F1 1e-6 --F2 1e-6 --F3 1e-6 \ #this is the itsxpress way
        #--max is the4 ITSx way

    hmmsearch -T 10 --domtblout $sample.domtbl.txt \
        --cpu 24 --tformat fasta \
        --max \
        ~/repo/itsxpress_2x_mod/itsxpress/ITSx_db/HMMs/F.hmm \
        $i
)
done

    
for i in *unique.fasta
do(

    sample=${i%.unique.fasta}
        #--F1 1e-6 --F2 1e-6 --F3 1e-6 \ #this is the itsxpress way
        #--max is the4 ITSx way

    perl ~/perl_scripts/get_seq_lens.pl $i > $sample.lens.txt
)
done
