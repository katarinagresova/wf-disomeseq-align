#!/bin/bash
# Driver job: runs snakemake with the SLURM profile, which submits one SLURM
# job per rule instance. Submit from the repo root (logs/ must exist) with:
#   sbatch slurm_job.sh [snakemake args], e.g. --configfile path/to/config.yaml

#SBATCH --job-name=wf-disomeseq-align
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate snake

snakemake --workflow-profile workflow/profiles/slurm -s workflow/Snakefile "$@"
