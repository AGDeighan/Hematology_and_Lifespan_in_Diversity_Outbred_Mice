# 2021-06-23



################################################################################
## load libraries ###########################################

options(stringsAsFactors = FALSE)
library(qtl2)
library(tidyverse)

#####


################################################################################
## helper functions ###########################################

scan1snps_perm <- function(
  genoprobs, map, pheno, kinship, addcovar, query_func, batch_length = 20, n_perm, cores
){
  SNP_INFO <- index_snps(
    map = map,
    snpinfo = query_func(
      chr = 1,
      start = 0,
      end = 0 + batch_length
    )
  )
  
  SNP_PROBS <- genoprob_to_snpprob(
    genoprobs = genoprobs, 
    snpinfo = SNP_INFO
  )
  
  PERM_RESULTS <- scan1perm(
    genoprobs = SNP_PROBS,
    pheno = pheno,
    kinship = kinship['1'],
    addcovar = addcovar,
    n_perm = n_perm,
    cores = cores
  )
  
  for(CHR in names(map)){
    for(POS in seq(0, max(map[[CHR]]), by = batch_length)){
      SNP_INFO <- index_snps(
        map = map,
        snpinfo = query_func(
          chr = CHR,
          start = POS,
          end = POS + batch_length
        )
      )
      
      SNP_PROBS <- genoprob_to_snpprob(
        genoprobs = genoprobs, 
        snpinfo = SNP_INFO
      )
      
      PERM_RESULTS <- pmax(
        PERM_RESULTS,
        scan1perm(
          genoprobs = SNP_PROBS,
          pheno = pheno,
          kinship = kinship[CHR],
          addcovar = addcovar,
          n_perm = n_perm,
          cores = cores
        )
      )
    }
  }
  
  return(PERM_RESULTS)
}


#####


################################################################################
## Load data etc #############

# load('data/processed/qtl_analysis_data/environments/LifespanHema_Environment.RData')
load('data/LifespanHema_Environment.RData')

#####


################################################################################
## Set up #############

SEED <- 501162 + 5
# QUERY_VAR <- create_variant_query_func(dbfile = 'data/original/genetic/cc_variants.sqlite')
QUERY_VAR <- create_variant_query_func(dbfile = 'data/cc_variants.sqlite')
PERMS <- 100
CORES <- 13
# OUTFILE_NAME <- paste0(
#   'results/VariantScans/varscan_perm_chunks/LifespanHema_',
#   PERMS, 'Perms_Seed',
#   SEED,
#   '.rds'
# )
OUTFILE_NAME <- paste0(
  'out/varscan_perm_chunks/LifespanHema_',
  PERMS, 'Perms_Seed',
  SEED,
  '.rds'
)

#####



################################################################################
## univariate genome scans #############

set.seed(seed = SEED)

cat('\nStarting unconditioned lifespan and hematology permutatation scans\n')
# Lifespan and hematology biomarkers
UNCONDITIONED_SCAN <- scan1snps_perm(
  genoprobs = genoprobs,
  map = map,
  pheno = data.matrix(
    dataset.LifespanHemaCovSexGen.phenotype$pheno[,c('Lifespan', 'rdw_07', 'hdw_07'), drop = FALSE]
  ),
  kinship = K,
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  n_perm = PERMS,
  batch_length = 20,
  query_func = QUERY_VAR,
  cores = CORES
)



#####



################################################################################
## save ###########################################

cat(paste0(
  '\nSaving results to file:\n',
  OUTFILE_NAME,
  '\n'
))

saveRDS(
  UNCONDITIONED_SCAN,
  OUTFILE_NAME
)


cat('\n\nDone\n#---------------------------------------------#\n\n\n\n\n')
#####



################################################################################
## clear workspace ###########################################

rm(list = ls())
graphics.off()

#####
