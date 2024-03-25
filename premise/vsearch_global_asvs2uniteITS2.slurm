#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vsearch"
#SBATCH --output=vsearch_global.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24

module purge
module load linuxbrew/colsa

cd /mnt/home/garnas/ewj4/blast_dbs/unite_07252023

fastaIn=/mnt/home/garnas/ewj4/EDRR_patho/itsxpress_mod_dada2_core/ASVs.fa
db=/mnt/home/garnas/ewj4/blast_dbs/unite_07252023/sh_general_allEuk_dynamic_singletons_25072023.ITS2.udb
outFile=/mnt/home/garnas/ewj4/EDRR_patho/itsxpress_mod_dada2_core/taxonomy_vsearch_global.tsv

#we are making user dfined output fields instead of using the --blast6out format.
# the b6 format has six fields that are meaningless with global nt alignments
#b6 fields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits
# 'qlo+qhi+tlo+thi are always 1 (lo) or aln length (hi) for both query and target
#evalue+bits are undefined for global pairwise nt alignments https://vcru.wisc.edu/simonlab/bioinformatics/programs/vsearch/vsearch_manual.pdf

# --userfields query+target+id+alnlen+mism+opens
# 1 = query ID, 2 = target ID, 3 = % identity (see --iddef), 4 = aln len, 5 = number mismatches, 6 = number gaps
# based on def %id should equal (alnlen - mism - opens)/alnlen

#also set --iddef to 1 or 4 to match BLAST output i.e., matching columns/alignment length
#the manual states that iddef1 == 4 (blast) in case of global pairwise alns
#note the default is matching columns/(alignment length - terminal gaps), which sure sounds like the BLAST def, but in any case we are trimming the seqs to positional homology so we will want to count terminal gaps as mismatches

# --maxaccepts determines the number of good matches (i.e., above the %id threshold) found before terminating the search (and the number of hits reported if maxaccepts is not set); --maxrejects determines the numbner of failed matches before termination. The search order is determined by first 'U sorting', i.e., sorting based on the number of shared 8-mers in attempt to look for the best hits first, and then results are reported in order of %ID. Because UNITE is already clustered the top hit *should* be better than all others (by the OTU clustering threshold of the DB) in cases where there is a *good* match (like within the clustering threshold or greater) but this fails when there are no good hits. In this case we may like to use LCA taxonomic assignment so we will want to report more than just the top hit. Increase --maxaccepts to 25 from the default 1, and maxrejects to 100 to increase sensitivity in the case of poor matches.

vsearch --usearch_global $fastaIn --userout $outFile \
    --userfields query+target+id+alnlen+mism+opens --db $db --threads 24 \
    --maxaccepts 25 --maxrejects 100 --id 0.75 --iddef 1 \
    --notrunclabels --output_no_hits --minseqlength 16
