# 2021-06-23



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

# CORES <- 6
CORES <- 13

# QUERY_VAR <- create_variant_query_func(dbfile = 'data/raw/genetic/cc_variants.sqlite')
QUERY_VAR <- create_variant_query_func(dbfile = 'data/cc_variants.sqlite')

#####


################################################################################
## univariate genome scans #############


cat('\nStarting unconditioned lifespan and hematology scan\n')
# Lifespan and hematology biomarkers
UNCONDITIONED_SCAN <- scan1snps(
  genoprobs = genoprobs, 
  map = map,
  pheno = data.matrix(
    dataset.LifespanHemaCovSexGen.phenotype$pheno[,c('Lifespan', 'rdw_07', 'hdw_07'), drop = FALSE]
  ),
  kinship = K, 
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  query_func = QUERY_VAR,
  chr = names(map),
  cores = CORES
)


cat('\nStarting RDW-conditioned lifespan scan\n')
# Lifespan conditioned on RDW
RDW_CONDITIONED_SCAN <- scan1snps(
  genoprobs = genoprobs, 
  map = map,
  pheno = data.matrix(
    dataset.LifespanCovSexGenRDW.phenotype$pheno[,-1, drop = FALSE]
  ),
  kinship = K, 
  addcovar = dataset.LifespanCovSexGenRDW.phenotype$covar,
  query_func = QUERY_VAR,
  chr = names(map),
  cores = CORES
)
colnames(RDW_CONDITIONED_SCAN$lod) <- c('LifespanCondRDW')


cat('\nStarting HDW-conditioned lifespan scan\n')
# Lifespan conditioned on HDW
HDW_CONDITIONED_SCAN <- scan1snps(
  genoprobs = genoprobs, 
  map = map,
  pheno = data.matrix(
    dataset.LifespanCovSexGenHDW.phenotype$pheno[,-1, drop = FALSE]
  ),
  kinship = K, 
  addcovar = dataset.LifespanCovSexGenHDW.phenotype$covar,
  query_func = QUERY_VAR,
  chr = names(map),
  cores = CORES
)
colnames(HDW_CONDITIONED_SCAN$lod) <- c('LifespanCondHDW')

cat('\nConcatenating scans\n')
# Combine scans
UNIVARIATE_SCAN <- list(
  lod = cbind(
    UNCONDITIONED_SCAN$lod, 
    RDW_CONDITIONED_SCAN$lod,
    HDW_CONDITIONED_SCAN$lod
  ),
  snpinfo = UNCONDITIONED_SCAN$snpinfo
)

#####



################################################################################
## save ###########################################

# cat('\nSaving scan results to:\n   results/VariantScans/LifespanHema_GenomeWideVariantScans.rds\n')
# saveRDS(
#   UNIVARIATE_SCAN,
#   'results/VariantScans/LifespanHema_GenomeWideVariantScans.rds'
# )
cat('\nSaving scan results to:\n   out/LifespanHema_GenomeWideVariantScans.rds\n')
saveRDS(
  UNIVARIATE_SCAN,
  'out/LifespanHema_GenomeWideVariantScans.rds'
)


cat('\n\nDone\n#---------------------------------------------#\n\n\n\n\n')
#####



################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####

