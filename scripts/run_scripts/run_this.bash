#!/bin/bash
#SBATCH --job-name=scrna_analysis
#SBATCH --output=/gpfs/scratch/ajn8926/spring_semester/scOmics/scripts/outputs/scrna_%j.out
#SBATCH --error=/gpfs/scratch/ajn8926/spring_semester/scOmics/scripts/outputs/scrna_%j.err
#SBATCH --time=0:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ajn8926@nyu.edu

# Check if script name was provided
if [ -z "$1" ]; then
    echo "Error: No script name provided"
    echo "Usage: sbatch run_analysis.bash <script_name>"
    exit 1
fi

SCRIPT_NAME="$1"
BASE_DIR="/gpfs/scratch/ajn8926/spring_semester/scOmics/scripts"

# Check if R script exists
if [ ! -f "${BASE_DIR}/R_scripts/${SCRIPT_NAME}.R" ]; then
    echo "Error: ${BASE_DIR}/R_scripts/${SCRIPT_NAME}.R not found"
    exit 1
fi

module load r/4.5.0
Rscript ${BASE_DIR}/R_scripts/${SCRIPT_NAME}.R

