# 20201-01-28
################################################################################
#
#
#   Author:  Andrew Deighan
#   E-mail: andrew.deighan@jax.org
#
################################################################################


################################################################################
## libraries etc #########################################################
set.seed(19940418)
options(stringsAsFactors = FALSE)
options(na.action = 'na.exclude')
library(qtl2)
library(qtl2convert)
library(tidyverse)

#####


################################################################################
## helper functions #########################################################

#####


################################################################################
## load phenotype data ##############################################

PHENODATA <- read_csv('data/processed/phenotypic/ClumpDateAdj_HemaLS.csv') %>% 
  mutate(
    Sex = factor(Sex),
    Generation = factor(
      Generation,
      levels = c('G7', 'G8', 'G9', 'G10', 'G11')
    ),
    Timepoint = factor(
      Timepoint, 
      levels = c('7 Months', '13 Months', '19 Months')
    )
  )

#####


################################################################################
## subset and clean the 8-state alleleprobs #############################

# The phenotype data has already been filtered to only include mice that
# lived more than 8 months, received the first round of CBC, and have genotype
# data. Thus, if we simple need to subset the genoprobs to include only mice
# present in the phenotype data.

# We will also exclude the Y chromosome, mitochondrial DNA, and the 
# pseudoautosomal region from analysis

ALLELE_PROBS <- readRDS('data/raw/genetic/Long_AlleleProbs.rds')

ALLELE_PROBS <- subset(ALLELE_PROBS, ind = unique(PHENODATA$MouseID))

ALLELE_PROBS <- subset(ALLELE_PROBS, chr = names(ALLELE_PROBS)[1:20])

ALLELE_PROBS <- clean_genoprob(
  ALLELE_PROBS,
  value_threshold = 0.001,
  column_threshold = 0.1
)

saveRDS(
  ALLELE_PROBS,
  'data/processed/genetic/AlleleProbs_8State.rds'
)

rm(ALLELE_PROBS)

#####


################################################################################
## subset and clean the 36-state genoprobs #############################

# The phenotype data has already been filtered to only include mice that
# lived more than 8 months, received the first round of CBC, and have genotype
# data. Thus, if we simple need to subset the genoprobs to include only mice
# present in the phenotype data.

# We will also exclude the Y chromosome, mitochondrial DNA, and the 
# pseudoautosomal region from analysis

GENO_PROBS <- readRDS('data/raw/genetic/Long_GenoProbs.rds')

GENO_PROBS <- subset(GENO_PROBS, ind = unique(PHENODATA$MouseID))

GENO_PROBS <- subset(GENO_PROBS, chr = names(GENO_PROBS)[1:20])

GENO_PROBS <- clean_genoprob(
  GENO_PROBS,
  value_threshold = 0.0005,
  column_threshold = 0.05
)

saveRDS(
  GENO_PROBS,
  'data/processed/genetic/GenoProbs_36State.rds'
)

#####


################################################################################
## load maps from cross2 object ##############################################

# Subset to only the autosomes and the X chromosome

CROSS2_OBJ <- readRDS('data/raw/genetic/QTL2_cross2_LongAndCS.rds')
CROSS2_OBJ <- subset(CROSS2_OBJ, chr = chr_names(CROSS2_OBJ)[1:20])
GMAP_LIST <- CROSS2_OBJ$gmap
PMAP_LIST <- CROSS2_OBJ$pmap

GMAP_DF <- map_list_to_df(GMAP_LIST)
PMAP_DF <- map_list_to_df(PMAP_LIST)

saveRDS(GMAP_LIST, 'data/processed/genetic/GeneticMap_List.rds')
saveRDS(GMAP_DF, 'data/processed/genetic/GeneticMap_Dataframe.rds')
saveRDS(PMAP_LIST, 'data/processed/genetic/PhysicalMap_List.rds')
saveRDS(PMAP_DF, 'data/processed/genetic/PhysicalMap_Dataframe.rds')


rm(CROSS2_OBJ, GMAP_DF, GMAP_LIST, PMAP_DF, PMAP_LIST)

#####


################################################################################
## Cleaning workspace #########################################################


rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####