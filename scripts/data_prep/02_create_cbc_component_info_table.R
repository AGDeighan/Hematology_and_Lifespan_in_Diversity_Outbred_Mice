# 2021-01-28
################################################################################
#
#   This script creates a table with information on all the CBC components
#
#   Author:  Andrew Deighan
#   E-mail: andrew.deighan@jax.org
#
################################################################################


################################################################################
## libraries etc #########################################################
options(stringsAsFactors = FALSE)
library(tidyverse)

#####


################################################################################
## helper functions #########################################################

#####


################################################################################
## Setup #########################################################

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

#####


################################################################################
## create TF data #########################################################

INFO_TABLE <- tibble(
  Phenotype = c(RBCphenos, WBCphenos, Platelets)
)

CBC_NAME <- c(
  'Red Blood Cell Count', 'Mean Hemoglobin Content', 
  'Corpuscular Hemoglobin Concentration Mean', 'Hemoglobin Distribution Width', 
  'Mean Corpuscular Volume', 'Red Cell Distribution Width', 
  'Hematocrit', 'Hemoglobin',
  'White Blood Cell Count', 'Neutrophil to Lymphocyte Ratio', 
  'Lymphocyte Count', 
  'Neutrophil Count', 
  'Monocyte Count', 
  'Eosinophil Count', 
  'Platelet Count', 'Mean Platelet Volume', 'Mean Platelet Mass'
)

ABBREVIATION <- c(
  'RBC count', 'CH', 
  'CHCM', 
  'HDW', 
  'MCV', 'RDW', 
  'Hct', 'Hgb',
  'WBC count', 'NLR',
  'Lymph. count', 
  'Neut. count', 
  'Mono. coumt', 
  'Eos. coumt', 
  'Platelet count', 'MPV', 'MPM'
)

FULL_NAME <- c(
  'RBC Count', 
  'Mean RBC Hgb Content (CH)',
  'Mean RBC Hgb Concentration (CHCM)',
  'Hemoglobin Distribution Width (HDW)',
  'Mean RBC Volume (MCV)', 'Red Cell Distribution Width (RDW)',
  'Hematocrit', 'Blood Hemoglobin',
  'WBC Count', 'Neutrophil to Lymphocyte Ratio (NLR)', 
  'Lymphocyte Count', 
  'Neutrophil Count', 
  'Monocyte Count', 
  'Eosinophil Count',
  'Platelet Count', 'Mean Platelet Volume (MPV)', 'Mean Platelet Mass (MPM)'
)

DICTIONARY <- read_csv('data/raw/phenotypic/JAC_DO_HemaLS_Phenotypes_descr.csv') %>% 
  filter(category == 'Hematology') %>% 
  mutate(
    Phenotype = gsub("\\.6", '', data_name),
    Phenotype = gsub("\\.12", '', Phenotype),
    Phenotype = gsub("\\.18", '', Phenotype)
  ) %>% 
  group_by(Phenotype) %>% 
  summarise(Units = unique(units)) %>% 
  ungroup()

INFO_TABLE <- INFO_TABLE %>% 
  left_join(DICTIONARY, by = 'Phenotype') %>% 
  mutate(
    TraditionalName = CBC_NAME,
    Abbreviation = ABBREVIATION,
    Name = FULL_NAME,
    Group = ifelse(
      Phenotype %in% RBCphenos, 'RBC',
      ifelse(
        Phenotype %in% WBCphenos, 'WBC',
        'Platelet'
      )
    ),
    Group = factor(Group, levels = c('RBC', 'WBC', 'Platelet')),
    Phenotype = factor(Phenotype, levels = c(RBCphenos, WBCphenos, Platelets))
  ) %>%
  arrange(Phenotype) %>% 
  select(Phenotype, Group, TraditionalName, Name, Abbreviation, Units) 

INFO_TABLE$Units[INFO_TABLE$Phenotype == 'hdw'] <- '%'
INFO_TABLE$Units[INFO_TABLE$Phenotype == 'nlr'] <- 'ratio'

INFO_TABLE %>% 
  saveRDS('data/processed/misc/ComponentInfoTable.rds')
INFO_TABLE %>% 
  write_csv('tables/misc/ComponentInfoTable.csv')

#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####
