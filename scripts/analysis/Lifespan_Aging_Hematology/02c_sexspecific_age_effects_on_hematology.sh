#!/bin/bash
#SBATCH -p compute -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 64 # number of cores
#SBATCH --mem=128GB # total memory
#SBATCH -t 72:00:00 # time (D-HH:MM)
#SBATCH -o slurm.%j.out # STDOUT
#SBATCH -e slurm.%j.err # STDERR

module load singularity
singularity exec do_hema_ls.sif Rscript 02c_sexspecific_age_effects_on_hematology.R
