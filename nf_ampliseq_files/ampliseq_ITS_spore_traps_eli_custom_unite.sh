#!/bin/bash
#SBATCH --account=PAS0247
#SBATCH --time=70:00:00
#SBATCH --mem=96G
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-nfc_ampliseq-%j.out

# Settings and constants
WORKFLOW_DIR=software/nfc-ampliseq/2_12_0

# Load the Nextflow Conda environment
module load miniconda3
source activate /users/PAS0247/torres704/.conda/envs/nextflow
export NXF_SINGULARITY_CACHEDIR=/users/PAS0247/torres704/containers

# Strict Bash settings
set -euo pipefail

# Process command-line arguments
input=$1
FW_primer=$2
RV_primer=$3
outdir=$4
dada_ref_tax_custom=$5

# Report
echo "Starting script ampliseq_ITS.sh"
date
echo "Samplesheet:                   $input"
echo "Forward primer:                $FW_primer"
echo "Reverse primer:                $RV_primer"
echo "Output dir:                    $outdir"
echo "Taxonomy db:                   $dada_ref_tax_custom"
echo

# Create the config file
echo "
process.executor = 'slurm'
process.clusterOptions='--account=PAS0247'
" > nextflow.config

# Run the workflow
nextflow run "$WORKFLOW_DIR" \
    --input "sample_sheet.csv" \
    --FW_primer "AACTTTYRRCAAYGGATCWCT" \
    --RV_primer "AGCCTCCGCTTATTGATATGCTTAART" \
    --outdir "./results_ITS" \
    --dada_ref_tax_custom "sh_general_allEuk_dynamic_singletons_25072023.ITS2.SPP_OF_CONCERN_CORRECTIONS_4.fasta" \
    --dada_ref_tax_custom_sp "sh_general_allEuk_dynamic_singletons_25072023.ITS2.SPP_OF_CONCERN_CORRECTIONS_4.fasta"\
    --illumina_pe_its \
    --illumina_novaseq \
    --skip_barrnap\
    -profile singularity \
    -ansi-log false \
    -resume

# Report
echo "Done with script ampliseq_ITS.sh"
date