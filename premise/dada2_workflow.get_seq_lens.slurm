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
    
for i in *unique.fasta
do(

    sample=${i%.unique.fasta}
    perl ~/perl_scripts/get_seq_lens.pl $i > $sample.lens.txt
)
done

cat *lens.txt > all.lens.txt


cd /mnt/home/garnas/ewj4/EDRR_patho/itsx_mod_test
    
for i in *unique.fasta
do(

    sample=${i%.unique.fasta}
    perl ~/perl_scripts/get_seq_lens.pl $i > $sample.lens.txt
)
done

cat *lens.txt > all_mod.lens.txt
