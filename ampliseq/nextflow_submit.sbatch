#!/bin/bash
##SBATCH --ntasks=1
#SBATCH --job-name="ampliseq"
#SBATCH --output=ampliseq.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24
##excluding a bunch bc they are often used by priority users and we don't
## need the special resources
##https://premise.sr.unh.edu/slurm.html#srun-excludenode117node118-hostname
##  the node range syntax isn't working so we will list it out
##SBATCH --exclude=node\[101-112,117-125\]
#SBATCH --exclude=node105,node106,node107,node108,node117,node118,node119,node167

## This method reserves 24 cpus via SLURM and then nextflow submits jobs to the
## reserved CPUs using local job submission approach (i.e., ignoring SLURM).
## Therefore nf param `process.executor` should be set to `local` (typically the default).
## We use this approach because using nextflow to manage SLURM submissions resulted
## in jobs hitting priority nodes and getting killed by priority user jobs in previous trials.

cd ~/test_dir/

module purge
module load anaconda/colsa
conda activate nf_ampliseq

export NXF_SINGULARITY_CACHEDIR=$HOME/nf_singularity_images
WORKFLOW_DIR=$HOME/nf_ampliseq_workflow/2_16_1

# Create the config file
echo "
process.executor = 'local'
" > nextflow.config

#ITS4-Fun
FWD_primer=AGCCTCCGCTTATTGATATGCTTAART
#5.8S-Fun
REV_primer=AACTTTYRRCAAYGGATCWCT

nextflow run "$WORKFLOW_DIR" \
    -profile singularity \
    --input "/mnt/gpfs01/home/garnas/ewj4/FEDRR_project/FEDRR_2025_all_samples_ampliseq.tsv" \
    --FW_primer "$FWD_primer" \
    --RV_primer "$REV_primer" \
    --outdir "/mnt/gpfs01/home/garnas/ewj4/FEDRR_project/FEDRR_2025_ampliseq_out" \
    --multiple_sequencing_runs \
    --ignore_empty_input_files \
    --illumina_pe_its \
    --dada_ref_tax_custom "/mnt/gpfs01/home/garnas/ewj4/blast_dbs/unite_07252023/sh_general_allEuk_dynamic_singletons_25072023.ITS2.SPP_OF_CONCERN_CORRECTIONS.fasta" \
    --skip_dada_addspecies \
    --cut_its 'its2' \
    --illumina_novaseq \
    --skip_qiime \
    --skip_barplot \
    --skip_abundance_tables \
    --skip_alpha_rarefaction \
    --skip_diversity_indices \
    --skip_ancom
    #--multiple_sequencing_runs

# note that running the test profile (-profile test,singularity) without additional flags resulted in an error
# at "NFCORE_AMPLISEQ:AMPLISEQ:QIIME2_BARPLOT":
# We don't care about the qiime parts so moving on, but if they are desired for some reason this will likely
# require troubleshooting to fix the qiime install
#
