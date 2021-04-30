# 2021-01-28
################################################################################
#
#   This script creates the RData environment that will be used for analysis of 
#   the hematology data.
#
#   Author:  Andrew Deighan
#   E-mail: andrew.deighan@jax.org
#
################################################################################

options(stringsAsFactors = FALSE)

library(tidyverse)

RBCphenos <- c(
  'rbc', 'ch', 'chcm', 'hdw',
  'mcv', 'rdw', 
  'hct', 'hgb'
)
WBCphenos <- c(
  'wbc', 'nlr',
  'n.lymph', 
  'n.neut',
  'n.mono',
  'n.eos'
)
Platelets <- c('plt', 'mpv', 'mpm')
CBCphenos <- c(RBCphenos, WBCphenos, Platelets)

Clump_Affected_Traits <- c('n.eos', Platelets)

Info_Table <- readRDS('data/processed/misc/ComponentInfoTable.rds')

PhenoData_ClumpDateAdj <- read_csv('data/processed/phenotypic/ClumpDateAdj_HemaLS.csv') %>% 
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

PhenoData_Raw <- read_csv('data/processed/phenotypic/Raw_HemaLS.csv') %>% 
  mutate(
    Sex = factor(Sex),
    Generation = factor(
      Generation,
      levels = c('G7', 'G8', 'G9', 'G10', 'G11')
    ),
    Timepoint = factor(
      Timepoint, 
      levels = c('7 Months', '13 Months', '19 Months')
    ),
    DOT = as.character(DOT)
  )

save(list = ls(), file = 'data/processed/phenotypic/CBC_AnalysisEnvironment.RData')

#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####