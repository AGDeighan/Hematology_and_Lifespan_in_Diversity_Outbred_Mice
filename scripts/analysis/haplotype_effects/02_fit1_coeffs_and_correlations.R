# 2021-06-22
################################################################################
#
#   This script 
# 
#
#   Author:  Andrew Deighan
#   E-mail: andrew.deighan@jax.org
#
################################################################################


################################################################################
## libraries etc ###########################################

library(qtl2)
library(tidyverse)

#####


################################################################################
## load data  ###########################################

load('data/processed/qtl_analysis_data/environments/LifespanHema_Environment.RData')
rm(
  dataset.LifespanCovSexGenHDW.phenotype, 
  dataset.LifespanCovSexGenNLR.phenotype, 
  dataset.LifespanCovSexGenRDW.phenotype
)
load('data/processed/phenotypic/CBC_AnalysisEnvironment.RData')


LIFESPAN_PEAK_TABLE <- read_csv(
  'tables/genome_wide_linkage_scans/Lifespan_Peak_Table.csv'
)
AGING_BIOM_PEAK_TABLE <- read_csv(
  'tables/genome_wide_linkage_scans/Aging_Biomarker_Hema_Peak_Table.csv'
)

#####


################################################################################
## setup  ###########################################

data.frame(LIFESPAN_PEAK_TABLE)
#   Chromosome  Position    Trait      LOD LODDropRDW LODDropHDW  LODDropNLR PValue  PeakMarker       lbSI      ubSI
# 1          2 128.76415 Lifespan 7.220318 -0.1324881  0.5770305 -0.04400962 0.0648  UNC3946248 128.372768 129.36445
# 2          6  25.22833 Lifespan 7.018060  0.3318909  0.3219315  0.54619685 0.0904 UNC10749252  23.581197  27.65261
# 3         18  10.77744 Lifespan 6.591569  1.3006534  1.9110734 -0.25081632 0.1853 UNC28718524   8.523683  14.39734

data.frame(AGING_BIOM_PEAK_TABLE)
#    Chromosome  Position Trait Timepoint       LOD  PValue         PeakMarker      lbSI      ubSI
# 1           7 105.60220   hdw        07 20.409209 < 1e-04 backupUNC070393214 101.66161 105.71390
# 2           7 102.38993   hdw        13 10.313110   2e-04        UNC13501160  99.62355 109.07832
# 3           7 102.76635   hdw        19  8.213306  0.0115        UNC13509119  98.35179 107.37202
# 4           9  73.87704   hdw        07  8.689333  0.0043        UNC16632607  69.49226  78.33424
# 5           9 108.09224   hdw        07 13.433244 < 1e-04        UNC17091640 107.58254 108.34486
# 6           9  74.90381   hdw        13  8.100343  0.0153        UNC16645564  71.48548  78.24372
# 7           9 108.06443   hdw        13 15.147539 < 1e-04        UNC17091435 107.58254 108.32187
# 8           9 108.14416   hdw        19 11.725464   1e-04        UNC17092085 107.60656 108.32187
# 9          13  34.30827   nlr        07  6.596819  0.1879        UNC22392281  20.13939  35.08892
# 10         18  17.06342   hdw        07  6.994611   0.096        UNC28802837  13.46265  24.62421
# 11         18  16.12121   hdw        13  7.638294   0.032        UNC28790404  13.40594  23.77372
# 12         18  15.83136   rdw        07  7.486349  0.0436       UNC180052940  14.05343  18.91441
# 13         18  17.64595   rdw        19  7.809898  0.0248 backupUNC180056465  13.46265  20.60528

#####

COEF_TABLE <- tibble(
  Peak = vector(),
  Marker = vector(),
  Phenotype = vector(),
  Founder = vector(),
  Coef = vector(),
  SE = vector()
)

################################################################################
## Get haplotype effects at lifespan peak #################################

MARKER <- 'UNC28718524'
PEAK <- 'Lifespan'



PHENO <- 'Lifespan'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = dataset.LifespanHemaCovSexGen.phenotype$pheno[,PHENO, drop = FALSE],
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = PHENO,
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '7 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '13 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '19 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '7 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '13 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '19 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)



rm(
  MARKER, PEAK, PHENO, TP,
  FITRESULTS
)

#####


################################################################################
## Get haplotype effects at 7-month RDW peak #################################

MARKER <- 'UNC180052940'
PEAK <- 'RDW at 7 Months'



PHENO <- 'Lifespan'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = dataset.LifespanHemaCovSexGen.phenotype$pheno[,PHENO, drop = FALSE],
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = PHENO,
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '7 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '13 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '19 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '7 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '13 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '19 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


rm(
  MARKER, PEAK, PHENO, TP,
  FITRESULTS
)

#####


################################################################################
## Get haplotype effects at 19-month RDW peak #################################

MARKER <- 'backupUNC180056465'
PEAK <- 'RDW at 19 Months'



PHENO <- 'Lifespan'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = dataset.LifespanHemaCovSexGen.phenotype$pheno[,PHENO, drop = FALSE],
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = PHENO,
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '7 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '13 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '19 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '7 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '13 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '19 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


rm(
  MARKER, PEAK, PHENO, TP,
  FITRESULTS
)

#####


################################################################################
## Get haplotype effects at 7-month HDW peak #################################

MARKER <- 'UNC28802837'
PEAK <- 'HDW at 7 Months'



PHENO <- 'Lifespan'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = dataset.LifespanHemaCovSexGen.phenotype$pheno[,PHENO, drop = FALSE],
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = PHENO,
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '7 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '13 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '19 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '7 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '13 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '19 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


rm(
  MARKER, PEAK, PHENO, TP,
  FITRESULTS
)

#####


################################################################################
## Get haplotype effects at 13-month HDW peak #################################

MARKER <- 'UNC28790404'
PEAK <- 'HDW at 13 Months'



PHENO <- 'Lifespan'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = dataset.LifespanHemaCovSexGen.phenotype$pheno[,PHENO, drop = FALSE],
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = PHENO,
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '7 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '13 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'RDW'
TP <- '19 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '7 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '13 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


PHENO <- 'HDW'
TP <- '19 Months'
FITRESULTS <- fit1(
  genoprobs = genoprobs[['18']][,,MARKER],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == TP) %>% 
    select(MouseID, all_of(tolower(PHENO))) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  se = TRUE,
  blup = TRUE
)
COEF_TABLE <- rbind(
  COEF_TABLE,
  tibble(
    Peak = PEAK,
    Marker = MARKER,
    Phenotype = paste(PHENO, TP, sep = ' at '),
    Founder = LETTERS[1:8],
    Coef = as.numeric(FITRESULTS$coef[LETTERS[1:8]]),
    SE = as.numeric(FITRESULTS$SE[LETTERS[1:8]])
  )
)


rm(
  MARKER, PEAK, PHENO, TP,
  FITRESULTS
)

#####


################################################################################
## Correlations of founder effects ###########################################

COEF_COR_TABLE <- COEF_TABLE %>% 
  select(
    -SE
  ) %>% 
  pivot_wider(
    values_from = 'Coef',
    names_from = 'Phenotype'
  ) %>% 
  group_by(
    Peak, Marker
  ) %>% 
  summarise(
    Lifespan.RDWat7mo_Correlation = cor.test(Lifespan, `RDW at 7 Months`)$estimate,
    Lifespan.RDWat7mo_PValue = cor.test(Lifespan, `RDW at 7 Months`)$p.value,
    Lifespan.RDWat13mo_Correlation = cor.test(Lifespan, `RDW at 13 Months`)$estimate,
    Lifespan.RDWat13mo_PValue = cor.test(Lifespan, `RDW at 13 Months`)$p.value,
    Lifespan.RDWat19mo_Correlation = cor.test(Lifespan, `RDW at 19 Months`)$estimate,
    Lifespan.RDWat19mo_PValue = cor.test(Lifespan, `RDW at 19 Months`)$p.value,
    
    Lifespan.HDWat7mo_Correlation = cor.test(Lifespan, `HDW at 7 Months`)$estimate,
    Lifespan.HDWat7mo_PValue = cor.test(Lifespan, `HDW at 7 Months`)$p.value,
    Lifespan.HDWat13mo_Correlation = cor.test(Lifespan, `HDW at 13 Months`)$estimate,
    Lifespan.HDWat13mo_PValue = cor.test(Lifespan, `HDW at 13 Months`)$p.value,
    Lifespan.HDWat19mo_Correlation = cor.test(Lifespan, `RDW at 19 Months`)$estimate,
    Lifespan.HDWat19mo_PValue = cor.test(Lifespan, `HDW at 19 Months`)$p.value,
  ) %>% 
  ungroup() %>% 
  pivot_longer(
    cols = Lifespan.RDWat7mo_Correlation:Lifespan.HDWat19mo_PValue,
    names_to = 'K',
    values_to = 'V'
  ) %>% 
  separate(
    K,
    into = c('Pair', 'K'),
    sep = '_'
  ) %>% 
  pivot_wider(
    names_from = K,
    values_from = V
  ) %>% 
  separate(
    Pair,
    into = c('PhenoA', 'PhenoB'),
    sep = '\\.'
  ) %>% 
  mutate(
    PhenoB = gsub('at', ' at ', PhenoB),
    PhenoB = gsub('mo', ' Months', PhenoB)
  )

#####


################################################################################
## Plot for peak of lifespan locus ###########################################

PLOT_DATA <- COEF_TABLE %>% 
  filter(
    Peak == 'Lifespan'
  ) %>% 
  filter(
    Phenotype != 'Lifespan'
  ) %>% 
  separate(
    Phenotype,
    into = c('CBCPheno', 'Timepoint'),
    sep = ' at '
  ) %>% 
  rename(
    CBC_Coef = Coef,
    CBC_SE = SE
  ) %>% 
  mutate(
    Timepoint = gsub('Mon', 'mon', Timepoint),
    Timepoint = factor(Timepoint, levels = unique(Timepoint)),
    CBCPheno = factor(CBCPheno, levels = unique(CBCPheno))
  ) %>% 
  left_join(
    COEF_TABLE %>% 
      filter(
        Peak == 'Lifespan'
      ) %>% 
      filter(
        Phenotype == 'Lifespan'
      ) %>% 
      rename(
        LS_Coef = Coef,
        LS_SE = SE
      ),
    by = c('Peak', 'Marker', 'Founder')
  )

COR_DATA <- COEF_COR_TABLE %>% 
  separate(
    PhenoB,
    into = c('CBCPheno', 'Timepoint'),
    sep = ' at '
  ) %>% 
  filter(
    Peak == 'Lifespan',
    PhenoA == 'Lifespan'
  ) %>% 
  mutate(
    Timepoint = gsub('Mon', 'mon', Timepoint),
    Timepoint = factor(Timepoint, levels = unique(Timepoint)),
    CBCPheno = factor(CBCPheno, levels = unique(CBCPheno)),
    YPos = ifelse(
      CBCPheno == 'RDW', 1.05, 0.55
    ),
    XPos = 1,
    Text = paste0(
      'Cor ~ ', signif(Correlation, digits = 3),
      '\n',
      '    p ~ ', signif(PValue, digits = 3)
    )
  )

PLOT <- PLOT_DATA %>% 
  ggplot() +
  theme_minimal(
    base_size = 9
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 3
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 3
  ) +
  geom_errorbarh(
    aes(
      xmin = LS_Coef - LS_SE, xmax = LS_Coef + LS_SE,
      y = CBC_Coef,
      color = Founder
    )
  ) +
  geom_errorbar(
    aes(
      x = LS_Coef,
      ymin = CBC_Coef - CBC_SE, ymax = CBC_Coef + CBC_SE,
      color = Founder
    )
  ) +
  geom_point(
    aes(
      x = LS_Coef, 
      y = CBC_Coef,
      color = Founder
    )
  ) +
  geom_text(
    data = COR_DATA,
    aes(
      x = XPos, y = YPos,
      label = Text
    ),
    hjust = 'left',
    vjust = 'top',
    size = 9 * 5 / 14
  ) +
  scale_color_manual(
    values = as.character(CCcolors),
    labels = names(CCcolors)
  ) +
  scale_x_continuous(
    limits = c(-6, 6),
    breaks = seq(-20, 20, 2),
    exp = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(-2, 2, 0.25)
  ) +
  facet_grid(
    CBCPheno ~ Timepoint,
    scales = 'free_y',
    switch = 'y',
    labeller = as_labeller(c(
      RDW = 'RDW (%)', HDW = 'HDW (%)',
      `7 months` = '7 months',
      `13 months` = '13 months',
      `19 months` = '19 months'
    ) )
  ) +
  theme(
    legend.position = 'top',
    strip.placement = 'outside',
    axis.title.x = element_text(size = 9),
    strip.text.y = element_text(size = 9)
  ) +
  labs(
    title = 'Effect of founder strain alleles on lifespan, RDW, and HDW',
    subtitle = 'Effects estimated at the peak marker of the locus for lifespan (chr 18, 10.78 Mb)',
    x = 'Lifespan (months)', 
    y = NULL
  )


pdf(
  'figures/haplotype_effects/allele_coefs_at_lifespan_locus_peak.pdf',
  width = 7, height = 5
)
plot(PLOT)
dev.off()

rm(PLOT_DATA, COR_DATA, PLOT)


#####


################################################################################
## Plot for peak of RDW at 7 months locus ###############################

PLOT_DATA <- COEF_TABLE %>% 
  filter(
    Peak == 'RDW at 7 Months'
  ) %>% 
  filter(
    Phenotype != 'Lifespan'
  ) %>% 
  separate(
    Phenotype,
    into = c('CBCPheno', 'Timepoint'),
    sep = ' at '
  ) %>% 
  rename(
    CBC_Coef = Coef,
    CBC_SE = SE
  ) %>% 
  mutate(
    Timepoint = gsub('Mon', 'mon', Timepoint),
    Timepoint = factor(Timepoint, levels = unique(Timepoint)),
    CBCPheno = factor(CBCPheno, levels = unique(CBCPheno))
  ) %>% 
  left_join(
    COEF_TABLE %>% 
      filter(
        Peak == 'RDW at 7 Months'
      ) %>% 
      filter(
        Phenotype == 'Lifespan'
      ) %>% 
      rename(
        LS_Coef = Coef,
        LS_SE = SE
      ),
    by = c('Peak', 'Marker', 'Founder')
  )

COR_DATA <- COEF_COR_TABLE %>% 
  separate(
    PhenoB,
    into = c('CBCPheno', 'Timepoint'),
    sep = ' at '
  ) %>% 
  filter(
    Peak == 'RDW at 7 Months',
    PhenoA == 'Lifespan'
  ) %>% 
  mutate(
    Timepoint = gsub('Mon', 'mon', Timepoint),
    Timepoint = factor(Timepoint, levels = unique(Timepoint)),
    CBCPheno = factor(CBCPheno, levels = unique(CBCPheno)),
    YPos = ifelse(
      CBCPheno == 'RDW', 1.05, 0.55
    ),
    XPos = 1,
    Text = paste0(
      'Cor ~ ', signif(Correlation, digits = 3),
      '\n',
      '    p ~ ', signif(PValue, digits = 3)
    )
  )

PLOT <- PLOT_DATA %>% 
  ggplot() +
  theme_minimal(
    base_size = 9
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 3
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 3
  ) +
  geom_errorbarh(
    aes(
      xmin = LS_Coef - LS_SE, xmax = LS_Coef + LS_SE,
      y = CBC_Coef,
      color = Founder
    )
  ) +
  geom_errorbar(
    aes(
      x = LS_Coef,
      ymin = CBC_Coef - CBC_SE, ymax = CBC_Coef + CBC_SE,
      color = Founder
    )
  ) +
  geom_point(
    aes(
      x = LS_Coef, 
      y = CBC_Coef,
      color = Founder
    )
  ) +
  geom_text(
    data = COR_DATA,
    aes(
      x = XPos, y = YPos,
      label = Text
    ),
    hjust = 'left',
    vjust = 'top',
    size = 9 * 5 / 14
  ) +
  scale_color_manual(
    values = as.character(CCcolors),
    labels = names(CCcolors)
  ) +
  scale_x_continuous(
    limits = c(-6, 6),
    breaks = seq(-20, 20, 2),
    exp = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(-2, 2, 0.25)
  ) +
  facet_grid(
    CBCPheno ~ Timepoint,
    scales = 'free_y',
    switch = 'y',
    labeller = as_labeller(c(
      RDW = 'RDW (%)', HDW = 'HDW (%)',
      `7 months` = '7 months',
      `13 months` = '13 months',
      `19 months` = '19 months'
    ) )
  ) +
  theme(
    legend.position = 'top',
    strip.placement = 'outside',
    axis.title.x = element_text(size = 9),
    strip.text.y = element_text(size = 9)
  ) +
  labs(
    title = 'Effect of founder strain alleles on lifespan, RDW, and HDW',
    subtitle = 'Effects estimated at the peak marker of the locus for RDW at 7 months (chr 18, 15.83 Mb)',
    x = 'Lifespan (months)', 
    y = NULL
  )


pdf(
  'figures/haplotype_effects/allele_coefs_at_RDW07_locus_peak.pdf',
  width = 7, height = 5
)
plot(PLOT)
dev.off()

rm(PLOT_DATA, COR_DATA, PLOT)

#####


################################################################################
## Plot for peak of RDW at 19 months locus ###############################

PLOT_DATA <- COEF_TABLE %>% 
  filter(
    Peak == 'RDW at 19 Months'
  ) %>% 
  filter(
    Phenotype != 'Lifespan'
  ) %>% 
  separate(
    Phenotype,
    into = c('CBCPheno', 'Timepoint'),
    sep = ' at '
  ) %>% 
  rename(
    CBC_Coef = Coef,
    CBC_SE = SE
  ) %>% 
  mutate(
    Timepoint = gsub('Mon', 'mon', Timepoint),
    Timepoint = factor(Timepoint, levels = unique(Timepoint)),
    CBCPheno = factor(CBCPheno, levels = unique(CBCPheno))
  ) %>% 
  left_join(
    COEF_TABLE %>% 
      filter(
        Peak == 'RDW at 19 Months'
      ) %>% 
      filter(
        Phenotype == 'Lifespan'
      ) %>% 
      rename(
        LS_Coef = Coef,
        LS_SE = SE
      ),
    by = c('Peak', 'Marker', 'Founder')
  )

COR_DATA <- COEF_COR_TABLE %>% 
  separate(
    PhenoB,
    into = c('CBCPheno', 'Timepoint'),
    sep = ' at '
  ) %>% 
  filter(
    Peak == 'RDW at 19 Months',
    PhenoA == 'Lifespan'
  ) %>% 
  mutate(
    Timepoint = gsub('Mon', 'mon', Timepoint),
    Timepoint = factor(Timepoint, levels = unique(Timepoint)),
    CBCPheno = factor(CBCPheno, levels = unique(CBCPheno)),
    YPos = ifelse(
      CBCPheno == 'RDW', 1.05, 0.55
    ),
    XPos = 1,
    Text = paste0(
      'Cor ~ ', signif(Correlation, digits = 3),
      '\n',
      '    p ~ ', signif(PValue, digits = 3)
    )
  )

PLOT <- PLOT_DATA %>% 
  ggplot() +
  theme_minimal(
    base_size = 9
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 3
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 3
  ) +
  geom_errorbarh(
    aes(
      xmin = LS_Coef - LS_SE, xmax = LS_Coef + LS_SE,
      y = CBC_Coef,
      color = Founder
    )
  ) +
  geom_errorbar(
    aes(
      x = LS_Coef,
      ymin = CBC_Coef - CBC_SE, ymax = CBC_Coef + CBC_SE,
      color = Founder
    )
  ) +
  geom_point(
    aes(
      x = LS_Coef, 
      y = CBC_Coef,
      color = Founder
    )
  ) +
  geom_text(
    data = COR_DATA,
    aes(
      x = XPos, y = YPos,
      label = Text
    ),
    hjust = 'left',
    vjust = 'top',
    size = 9 * 5 / 14
  ) +
  scale_color_manual(
    values = as.character(CCcolors),
    labels = names(CCcolors)
  ) +
  scale_x_continuous(
    limits = c(-6, 6),
    breaks = seq(-20, 20, 2),
    exp = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(-2, 2, 0.25)
  ) +
  facet_grid(
    CBCPheno ~ Timepoint,
    scales = 'free_y',
    switch = 'y',
    labeller = as_labeller(c(
      RDW = 'RDW (%)', HDW = 'HDW (%)',
      `7 months` = '7 months',
      `13 months` = '13 months',
      `19 months` = '19 months'
    ) )
  ) +
  theme(
    legend.position = 'top',
    strip.placement = 'outside',
    axis.title.x = element_text(size = 9),
    strip.text.y = element_text(size = 9)
  ) +
  labs(
    title = 'Effect of founder strain alleles on lifespan, RDW, and HDW',
    subtitle = 'Effects estimated at the peak marker of the locus for RDW at 19 months (chr 18, 17.65 Mb)',
    x = 'Lifespan (months)', 
    y = NULL
  )


pdf(
  'figures/haplotype_effects/allele_coefs_at_RDW19_locus_peak.pdf',
  width = 7, height = 5
)
plot(PLOT)
dev.off()

rm(PLOT_DATA, COR_DATA, PLOT)

#####


################################################################################
## Plot for peak of HDW at 7 months locus ###############################

PLOT_DATA <- COEF_TABLE %>% 
  filter(
    Peak == 'HDW at 7 Months'
  ) %>% 
  filter(
    Phenotype != 'Lifespan'
  ) %>% 
  separate(
    Phenotype,
    into = c('CBCPheno', 'Timepoint'),
    sep = ' at '
  ) %>% 
  rename(
    CBC_Coef = Coef,
    CBC_SE = SE
  ) %>% 
  mutate(
    Timepoint = gsub('Mon', 'mon', Timepoint),
    Timepoint = factor(Timepoint, levels = unique(Timepoint)),
    CBCPheno = factor(CBCPheno, levels = unique(CBCPheno))
  ) %>% 
  left_join(
    COEF_TABLE %>% 
      filter(
        Peak == 'HDW at 7 Months'
      ) %>% 
      filter(
        Phenotype == 'Lifespan'
      ) %>% 
      rename(
        LS_Coef = Coef,
        LS_SE = SE
      ),
    by = c('Peak', 'Marker', 'Founder')
  )

COR_DATA <- COEF_COR_TABLE %>% 
  separate(
    PhenoB,
    into = c('CBCPheno', 'Timepoint'),
    sep = ' at '
  ) %>% 
  filter(
    Peak == 'HDW at 7 Months',
    PhenoA == 'Lifespan'
  ) %>% 
  mutate(
    Timepoint = gsub('Mon', 'mon', Timepoint),
    Timepoint = factor(Timepoint, levels = unique(Timepoint)),
    CBCPheno = factor(CBCPheno, levels = unique(CBCPheno)),
    YPos = ifelse(
      CBCPheno == 'RDW', 1.05, 0.55
    ),
    XPos = 1,
    Text = paste0(
      'Cor ~ ', signif(Correlation, digits = 3),
      '\n',
      '    p ~ ', signif(PValue, digits = 3)
    )
  )

PLOT <- PLOT_DATA %>% 
  ggplot() +
  theme_minimal(
    base_size = 9
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 3
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 3
  ) +
  geom_errorbarh(
    aes(
      xmin = LS_Coef - LS_SE, xmax = LS_Coef + LS_SE,
      y = CBC_Coef,
      color = Founder
    )
  ) +
  geom_errorbar(
    aes(
      x = LS_Coef,
      ymin = CBC_Coef - CBC_SE, ymax = CBC_Coef + CBC_SE,
      color = Founder
    )
  ) +
  geom_point(
    aes(
      x = LS_Coef, 
      y = CBC_Coef,
      color = Founder
    )
  ) +
  geom_text(
    data = COR_DATA,
    aes(
      x = XPos, y = YPos,
      label = Text
    ),
    hjust = 'left',
    vjust = 'top',
    size = 9 * 5 / 14
  ) +
  scale_color_manual(
    values = as.character(CCcolors),
    labels = names(CCcolors)
  ) +
  scale_x_continuous(
    limits = c(-6, 6),
    breaks = seq(-20, 20, 2),
    exp = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(-2, 2, 0.25)
  ) +
  facet_grid(
    CBCPheno ~ Timepoint,
    scales = 'free_y',
    switch = 'y',
    labeller = as_labeller(c(
      RDW = 'RDW (%)', HDW = 'HDW (%)',
      `7 months` = '7 months',
      `13 months` = '13 months',
      `19 months` = '19 months'
    ) )
  ) +
  theme(
    legend.position = 'top',
    strip.placement = 'outside',
    axis.title.x = element_text(size = 9),
    strip.text.y = element_text(size = 9)
  ) +
  labs(
    title = 'Effect of founder strain alleles on lifespan, RDW, and HDW',
    subtitle = 'Effects estimated at the peak marker of the locus for HDW at 7 months (chr 18, 17.06 Mb)',
    x = 'Lifespan (months)', 
    y = NULL
  )


pdf(
  'figures/haplotype_effects/allele_coefs_at_HDW07_locus_peak.pdf',
  width = 7, height = 5
)
plot(PLOT)
dev.off()

rm(PLOT_DATA, COR_DATA, PLOT)

#####


################################################################################
## Plot for peak of HDW at 13 months locus ###############################

PLOT_DATA <- COEF_TABLE %>% 
  filter(
    Peak == 'HDW at 13 Months'
  ) %>% 
  filter(
    Phenotype != 'Lifespan'
  ) %>% 
  separate(
    Phenotype,
    into = c('CBCPheno', 'Timepoint'),
    sep = ' at '
  ) %>% 
  rename(
    CBC_Coef = Coef,
    CBC_SE = SE
  ) %>% 
  mutate(
    Timepoint = gsub('Mon', 'mon', Timepoint),
    Timepoint = factor(Timepoint, levels = unique(Timepoint)),
    CBCPheno = factor(CBCPheno, levels = unique(CBCPheno))
  ) %>% 
  left_join(
    COEF_TABLE %>% 
      filter(
        Peak == 'HDW at 13 Months'
      ) %>% 
      filter(
        Phenotype == 'Lifespan'
      ) %>% 
      rename(
        LS_Coef = Coef,
        LS_SE = SE
      ),
    by = c('Peak', 'Marker', 'Founder')
  )

COR_DATA <- COEF_COR_TABLE %>% 
  separate(
    PhenoB,
    into = c('CBCPheno', 'Timepoint'),
    sep = ' at '
  ) %>% 
  filter(
    Peak == 'HDW at 13 Months',
    PhenoA == 'Lifespan'
  ) %>% 
  mutate(
    Timepoint = gsub('Mon', 'mon', Timepoint),
    Timepoint = factor(Timepoint, levels = unique(Timepoint)),
    CBCPheno = factor(CBCPheno, levels = unique(CBCPheno)),
    YPos = ifelse(
      CBCPheno == 'RDW', 1.05, 0.55
    ),
    XPos = 1,
    Text = paste0(
      'Cor ~ ', signif(Correlation, digits = 3),
      '\n',
      '    p ~ ', signif(PValue, digits = 3)
    )
  )

PLOT <- PLOT_DATA %>% 
  ggplot() +
  theme_minimal(
    base_size = 9
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 3
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 3
  ) +
  geom_errorbarh(
    aes(
      xmin = LS_Coef - LS_SE, xmax = LS_Coef + LS_SE,
      y = CBC_Coef,
      color = Founder
    )
  ) +
  geom_errorbar(
    aes(
      x = LS_Coef,
      ymin = CBC_Coef - CBC_SE, ymax = CBC_Coef + CBC_SE,
      color = Founder
    )
  ) +
  geom_point(
    aes(
      x = LS_Coef, 
      y = CBC_Coef,
      color = Founder
    )
  ) +
  geom_text(
    data = COR_DATA,
    aes(
      x = XPos, y = YPos,
      label = Text
    ),
    hjust = 'left',
    vjust = 'top',
    size = 9 * 5 / 14
  ) +
  scale_color_manual(
    values = as.character(CCcolors),
    labels = names(CCcolors)
  ) +
  scale_x_continuous(
    limits = c(-6, 6),
    breaks = seq(-20, 20, 2),
    exp = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(-2, 2, 0.25)
  ) +
  facet_grid(
    CBCPheno ~ Timepoint,
    scales = 'free_y',
    switch = 'y',
    labeller = as_labeller(c(
      RDW = 'RDW (%)', HDW = 'HDW (%)',
      `7 months` = '7 months',
      `13 months` = '13 months',
      `19 months` = '19 months'
    ) )
  ) +
  theme(
    legend.position = 'top',
    strip.placement = 'outside',
    axis.title.x = element_text(size = 9),
    strip.text.y = element_text(size = 9)
  ) +
  labs(
    title = 'Effect of founder strain alleles on lifespan, RDW, and HDW',
    subtitle = 'Effects estimated at the peak marker of the locus for HDW at 13 months (chr 18, 16.12 Mb)',
    x = 'Lifespan (months)', 
    y = NULL
  )


pdf(
  'figures/haplotype_effects/allele_coefs_at_HDW13_locus_peak.pdf',
  width = 7, height = 5
)
plot(PLOT)
dev.off()

rm(PLOT_DATA, COR_DATA, PLOT)

#####


################################################################################
##  ###########################################
13

#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####
