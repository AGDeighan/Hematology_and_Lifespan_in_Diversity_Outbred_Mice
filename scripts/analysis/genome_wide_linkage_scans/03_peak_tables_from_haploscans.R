# 2021-01-07



################################################################################
## load libraries ###########################################

options(stringsAsFactors = FALSE)
library(qtl2)
library(tidyverse)

#####


################################################################################
## accessory functions #############



#####



################################################################################
## Load data etc #############


MAP <- readRDS('data/working/genetic/PhysicalMap_List.rds')
ADD_SCANS <- readRDS('results/HaplotypeGenomeScans/LifespanHema_HaploScans.rds')
PERM_SCANS <- readRDS('results/HaplotypeGenomeScans/LifespanHema_10000Perms.rds')



#####


################################################################################
## Calculate permutation thresholds #############

ALPHA <- c(0.05, 0.2)
PERM_THRESH <- summary_scan1perm(PERM_SCANS, alpha = ALPHA)
#      Lifespan  rdw  hdw  nlr
# 0.05     7.38 7.37 7.36 7.30
# 0.2      6.53 6.53 6.50 6.49

#####


################################################################################
## Create peak table  #############

PEAK_DROP <- 5
LOD_DROP <- 1.5

(
  LIFESPAN_PEAK_TABLE <- find_peaks(
    subset_scan1(ADD_SCANS, lodcolumn = 'Lifespan'),
    map = MAP,
    threshold = PERM_THRESH['0.2', 'Lifespan'],
    peakdrop = PEAK_DROP,
    drop = LOD_DROP
  ) %>% 
    rename(
      Trait = lodcolumn,
      Chromosome = chr,
      Position = pos,
      LOD = lod,
      lbSI = ci_lo,
      ubSI = ci_hi
    ) %>% 
    rowwise() %>% 
    mutate(
      Threshold80 = PERM_THRESH[2,Trait],
      Threshold95 = PERM_THRESH[1,Trait],
      PValue = mean(PERM_SCANS[,Trait] > LOD, na.rm = TRUE),
      PValue = ifelse(
        PValue == 0, 
        paste0('< ', 1/nrow(PERM_SCANS)), 
        as.character(PValue)
      ),
      PeakMarker = find_marker(map = MAP, chr = Chromosome, pos = Position),
      Chromosome = as.character(Chromosome),
      LODDropRDW = LOD - ADD_SCANS[PeakMarker, 'LifespanCondRDW'],
      LODDropHDW = LOD - ADD_SCANS[PeakMarker, 'LifespanCondHDW'],
      LODDropNLR = LOD - ADD_SCANS[PeakMarker, 'LifespanCondNLR']
    ) %>% 
    select(
      Chromosome, Position, Trait, 
      LOD, LODDropRDW, LODDropHDW, LODDropNLR, PValue, 
      PeakMarker, lbSI, ubSI
    ) %>%
    arrange(as.numeric(Chromosome), Trait, Position) %>% 
    data.frame()
)
#   Chromosome  Position    Trait      LOD   LODDropRDW LODDropHDW  LODDropNLR PValue  PeakMarker       lbSI      ubSI
# 1          2 128.76415 Lifespan 7.572850 -0.228302174 0.59122269 -0.01117861 0.0364  UNC3946248 128.372768 129.29544
# 2          6  25.22833 Lifespan 7.095377  0.334361187 0.31112368  0.62115483 0.0809 UNC10749252  23.581197  27.29246
# 3         11  65.09165 Lifespan 6.793891 -0.005805706 0.06052521  0.25390365 0.1351 UNC19805699  63.074182  67.33882
# 4         18  10.72873 Lifespan 6.540183  1.240995612 1.86692308 -0.06789056 0.1957 JAX00080317   8.589508  14.39734

(
  HEMA_PEAK_TABLE <- find_peaks(
    subset_scan1(ADD_SCANS, lodcolumn = c('rdw', 'hdw', 'nlr')),
    map = MAP,
    threshold = PERM_THRESH['0.2', c('rdw', 'hdw', 'nlr')],
    peakdrop = PEAK_DROP,
    drop = LOD_DROP
  ) %>% 
    rename(
      Trait = lodcolumn,
      Chromosome = chr,
      Position = pos,
      LOD = lod,
      lbSI = ci_lo,
      ubSI = ci_hi
    ) %>% 
    rowwise() %>% 
    mutate(
      PValue = mean(PERM_SCANS[,Trait] > LOD, na.rm = TRUE),
      PValue = ifelse(
        PValue == 0, 
        paste0('< ', 1/nrow(PERM_SCANS)), 
        as.character(PValue)
      ),
      PeakMarker = find_marker(map = MAP, chr = Chromosome, pos = Position),
      Chromosome = as.character(Chromosome)
    ) %>% 
    select(
      Chromosome, Position, Trait, 
      LOD, PValue, 
      PeakMarker, lbSI, ubSI
    ) %>%
    arrange(as.numeric(Chromosome), Trait, Position) %>% 
    data.frame()
)
#   Chromosome  Position Trait       LOD  PValue  PeakMarker      lbSI      ubSI
# 1          7 102.08464   hdw 20.021751 < 1e-04 UNC13497365 101.66161 105.70587
# 2          9  75.30248   hdw  9.035527   0.002 UNC16650140  69.49213  78.33424
# 3          9 108.09224   hdw 13.572662 < 1e-04 UNC17091640 107.58254 108.34486
# 4         13  31.66024   nlr  6.678560  0.1487 UNC22367602  20.14612  34.67895
# 5         18  16.58863   hdw  6.961403  0.0984 UNC28796722  13.46265  24.83662
# 6         18  16.13494   rdw  7.344930  0.0515 JAX00452857  14.04108  18.99703

#####


################################################################################
## save ###########################################

write_csv(
  LIFESPAN_PEAK_TABLE,
  'tables/HaplotypeGenomeScans/Lifespan_Peak_Table.csv'
)

write_csv(
  HEMA_PEAK_TABLE,
  'tables/HaplotypeGenomeScans/Hema_Peak_Table.csv'
)

#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####

