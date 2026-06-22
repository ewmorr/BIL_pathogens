#!/bin/bash
#SBATCH --account=PAS0247
#SBATCH --mail-type=END,FAIL

# First load OSC's (mini)Conda module
module load miniconda3
# Then activate the Nextflow conda environment 
source activate /users/PAS0247/torres704/.conda/envs/nextflow
# Create an environment variable for the container directory
export NXF_SINGULARITY_CACHEDIR=/fs/ess/PAS0247/containers

#Download the nf-core ampliseq pipeline 
nf-core pipelines download ampliseq \
    --revision 2.12.0 \
    --outdir software_2/nfc-ampliseq\
    --compress none \
    --container-system singularity \
    --singularity-cache-only
    
    --container-cache-utilisation amend \

# Define the input and output dir 
outdir=results/nfc_ampliseq
input="sample_sheet.csv"
#work-dir=results/raw
# Define other paramethers
FW_primer="CTTGGTCATTTAGAGGAAGTAA"
RV_primer="GCTGCGTTCTTCATCGATGC"
#dada_ref_taxonomy='unite-fungi'
#extension="/*.raw_{1,2}.fastq.gz"
