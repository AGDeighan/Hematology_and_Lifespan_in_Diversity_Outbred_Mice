# 2021-06-11
################################################################################
#
#   This script performs haplotype (linkage) scans for lifespan and each CBC
#   trait across the genome. Sex and generation are included as covariates.
#   For lifespan, an additional three scans are performed. In each of which one
#   of the aging biomarker CBC traits (RDW, HDW, and NLR) are also included as
#   a covariate.
# 
#   This script is designed to by run on a high performance computing cluster 
#   with at least 13 cores.
#
#   Author:  Andrew Deighan
#   E-mail: andrew.deighan@jax.org
#
################################################################################



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

CORES <- 13

#####


################################################################################
## univariate genome scans #############

# Lifespan and hematology biomarkers
UNCONDITIONED_SCAN <- scan1(
  genoprobs = genoprobs, 
  pheno = data.matrix(
    dataset.LifespanHemaCovSexGen.phenotype$pheno[,-1, drop = FALSE]
  ),
  kinship = K, 
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  cores = CORES
)

# Lifespan conditioned on RDW
RDW_CONDITIONED_SCAN <- scan1(
  genoprobs = genoprobs, 
  pheno = data.matrix(
    dataset.LifespanCovSexGenRDW.phenotype$pheno[,-1, drop = FALSE]
  ),
  kinship = K, 
  addcovar = dataset.LifespanCovSexGenRDW.phenotype$covar,
  cores = CORES
)
colnames(RDW_CONDITIONED_SCAN) <- c('LifespanCondRDW')

# Lifespan conditioned on HDW
HDW_CONDITIONED_SCAN <- scan1(
  genoprobs = genoprobs, 
  pheno = data.matrix(
    dataset.LifespanCovSexGenHDW.phenotype$pheno[,-1, drop = FALSE]
  ),
  kinship = K, 
  addcovar = dataset.LifespanCovSexGenHDW.phenotype$covar,
  cores = CORES
)
colnames(HDW_CONDITIONED_SCAN) <- c('LifespanCondHDW')

# Lifespan conditioned on NLR
NLR_CONDITIONED_SCAN <- scan1(
  genoprobs = genoprobs, 
  pheno = data.matrix(
    dataset.LifespanCovSexGenNLR.phenotype$pheno[,-1, drop = FALSE]
  ),
  kinship = K, 
  addcovar = dataset.LifespanCovSexGenNLR.phenotype$covar,
  cores = CORES
)
colnames(NLR_CONDITIONED_SCAN) <- c('LifespanCondNLR')

# Combine scans
UNIVARIATE_SCAN <- cbind(
  UNCONDITIONED_SCAN, 
  RDW_CONDITIONED_SCAN,
  HDW_CONDITIONED_SCAN,
  NLR_CONDITIONED_SCAN
)

#####



################################################################################
## save ###########################################


# saveRDS(
#   UNIVARIATE_SCAN,
#   'results/genome_wide_linkage_scans/LifespanHema_HaploScans.rds'
# )

saveRDS(
  UNIVARIATE_SCAN,
  'out/LifespanHema_HaploScans.rds'
)

#####



################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####

