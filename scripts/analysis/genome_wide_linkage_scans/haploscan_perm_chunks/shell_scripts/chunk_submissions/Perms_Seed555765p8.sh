#!/bin/bash
#SBATCH -p compute -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 14 # number of cores
#SBATCH --mem=64GB # total memory
#SBATCH -t 72:00:00 # time (D-HH:MM)
#SBATCH -o slurm.%j.out # STDOUT
#SBATCH -e slurm.%j.err # STDERR

module load singularity
singularity exec do_hema_ls.sif Rscript haploscan_perm_chunks/LifespanHema_Perms_Seed555765p8.R


