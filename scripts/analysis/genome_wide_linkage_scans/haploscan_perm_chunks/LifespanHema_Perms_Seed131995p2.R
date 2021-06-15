# 2021-06-11



################################################################################
## load libraries ###########################################

options(stringsAsFactors = FALSE)
library(qtl2)
library(tidyverse)

#####



################################################################################
## Load data etc #############


# load('data/processed/qtl_analysis_data/environments/LifespanHema_Environment.RData')

load('data/LifespanHema_Environment.RData')


#####


################################################################################
## Set up #############

SEED <- 131995 + 2
PERMS <- 100
CORES <- 13
# OUTFILE_UNCOND_NAME <- paste0(
#   'results/genome_wide_linkage_scans/haploscan_perm_chunks/LifespanHema_',
#   PERMS, 'Perms_Seed',
#   SEED,
#   '.rds'
# )
OUTFILE_UNCOND_NAME <- paste0(
  'out/haploscan_perm_chunks/LifespanHema_',
  PERMS, 'Perms_Seed',
  SEED,
  '.rds'
)


#####


################################################################################
## univariate genome scans #############

cat('\n Starting unconditioned scan.\n')

set.seed(seed = SEED)

UNCOND_SCAN <- scan1perm(
  genoprobs = genoprobs, 
  pheno = data.matrix(
    dataset.LifespanHemaCovSexGen.phenotype$pheno[,-1,drop = FALSE]
  ),
  kinship = K, 
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  n_perm = PERMS,
  cores = CORES
)

cat(paste0(
  '\nSaving unconditioned scan to file:\n',
  OUTFILE_UNCOND_NAME,
  '\n'
))

saveRDS(
  UNCOND_SCAN,
  file = OUTFILE_UNCOND_NAME
)

#####



################################################################################
## clear workspace ###########################################

rm(list = ls())
graphics.off()

#####

