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

MAP <- readRDS('data/processed/genetic/PhysicalMap_List.rds')
ADD_SCANS <- readRDS('results/genome_wide_linkage_scans/LifespanHema_HaploScans.rds')
PERM_SCANS <- readRDS('results/genome_wide_linkage_scans/LifespanHema_10000Perms.rds')

#####


################################################################################
## Calculate permutation thresholds #############

ALPHA <- c(0.05, 0.2)
PERM_THRESH <- summary_scan1perm(PERM_SCANS, alpha = ALPHA)
#      Lifespan  rdw  hdw  nlr
# 0.05     7.38 7.37 7.36 7.30
# 0.2      6.53 6.53 6.50 6.49


# LOD thresholds (10000 permutations)
#       Lifespan    rbc_07     rbc_13     rbc_19      ch_07     ch_13      ch_19    chcm_07    chcm_13   chcm_19    hdw_07    hdw_13   hdw_19
# 0.05      7.36      7.37       7.36       7.43       7.35      7.39       7.46       7.35       7.38      7.39      7.36      7.39     7.39
# 0.2       6.55      6.54       6.53       6.58       6.50      6.56       6.56       6.54       6.55      6.57      6.53      6.55     6.59
#          mcv_07   mcv_13     mcv_19     rdw_07     rdw_13    rdw_19     hct_07     hct_13     hct_19    hgb_07    hgb_13    hgb_19
# 0.05       7.36     7.40       7.42       7.39       7.43      7.41       7.41       7.42       7.39      7.36      7.40      7.46
# 0.2        6.51     6.55       6.57       6.54       6.57      6.56       6.55       6.54       6.56      6.54      6.55      6.56
#          wbc_07   wbc_13     wbc_19     nlr_07     nlr_13    nlr_19 n.lymph_07 n.lymph_13 n.lymph_19 n.neut_07 n.neut_13 n.neut_19 n.mono_07 n.mono_13 n.mono_19 
# 0.05       7.40     7.35       7.38       7.38       7.38      7.45       7.35       7.37       7.40      7.41      7.38      7.39      7.35      7.40      7.36 
# 0.2        6.58     6.54       6.55       6.56       6.51      6.60       6.55       6.54       6.57      6.55      6.54      6.54      6.53      6.55      6.56 
#        n.eos_07 n.eos_13   n.eos_19     plt_07    plt_13     plt_19     mpv_07     mpv_13     mpv_19    mpm_07    mpm_13    mpm_19
# 0.05       7.37     7.39       7.45       7.41      7.40       7.44       7.40       7.39       7.43      7.31      7.38      7.40
# 0.2        6.55     6.54       6.58       6.54      6.55       6.58       6.54       6.53       6.57      6.51      6.52      6.55

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
#   Chromosome  Position    Trait      LOD LODDropRDW LODDropHDW  LODDropNLR PValue  PeakMarker       lbSI      ubSI
# 1          2 128.76415 Lifespan 7.220318 -0.1324881  0.5770305 -0.04400962 0.0648  UNC3946248 128.372768 129.36445
# 2          6  25.22833 Lifespan 7.018060  0.3318909  0.3219315  0.54619685 0.0904 UNC10749252  23.581197  27.65261
# 3         18  10.77744 Lifespan 6.591569  1.3006534  1.9110734 -0.25081632 0.1853 UNC28718524   8.523683  14.39734

(
  AGING_BIOM_PEAK_TABLE <- find_peaks(
    subset_scan1(
      ADD_SCANS, 
      lodcolumn = paste(
        rep(c('rdw', 'hdw', 'nlr'), each = 3), 
        c('07', '13', '19'),
        sep = '_'
      )
    ),
    map = MAP,
    threshold = PERM_THRESH[
      '0.2', 
      paste(
        rep(c('rdw', 'hdw', 'nlr'), each = 3), 
        c('07', '13', '19'),
        sep = '_'
      )
    ],
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
    separate(
      Trait,
      into = c('Trait', 'Timepoint'),
      sep = '_'
    ) %>% 
    select(
      Chromosome, Position, Trait, Timepoint,
      LOD, PValue, 
      PeakMarker, lbSI, ubSI
    ) %>%
    arrange(as.numeric(Chromosome), Trait, Timepoint, Position) %>% 
    data.frame()
)
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

(
  OTHER_HEMA_PEAK_TABLE <- find_peaks(
    subset_scan1(
      ADD_SCANS, 
      lodcolumn = colnames(PERM_THRESH)[-1]
    ),
    map = MAP,
    threshold = PERM_THRESH[
      '0.2', 
      colnames(PERM_THRESH)[-1]
    ],
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
    separate(
      Trait,
      into = c('Trait', 'Timepoint'),
      sep = '_'
    ) %>% 
    select(
      Chromosome, Position, Trait, Timepoint,
      LOD, PValue, 
      PeakMarker, lbSI, ubSI
    ) %>%
    arrange(as.numeric(Chromosome), Trait, Timepoint, Position) %>% 
    data.frame()
)
#    Chromosome  Position   Trait Timepoint       LOD  PValue         PeakMarker       lbSI      ubSI
# 1           1 158.23139   n.eos        13  7.138857   0.076         UNC2031646 156.497611 166.88328
# 2           1  78.64981  n.neut        19  6.932242  0.1103          UNC995947  76.523368  84.23175
# 3           1 143.73991     wbc        07  7.153857  0.0775         UNC1841019 135.816613 147.08021
# 4           2 165.49863    chcm        13  7.133569  0.0772         UNC4418426 165.289259 165.88089
# 5           2 110.26378     hgb        07  6.584675  0.1865         UNC3744236 108.384206 168.92841
# 6           3  57.15931     mpv        07 10.995629   1e-04        JAX00107762  53.040565  63.03999
# 7           3  61.39382     mpv        13  6.670517  0.1606         UNC5375208  53.040565  63.03999
# 8           3  61.39382     mpv        19  9.433236  0.0013         UNC5375208  56.339548  63.03999
# 9           4  58.60305     hct        07  6.814430  0.1339         UNC7334884  46.537864  63.47401
# 10          4  30.78690  n.neut        07  6.967874  0.1033         UNC6994954  24.596410  33.87110
# 11          4 138.46331  n.neut        07  6.584368  0.1918         UNC8328214 133.650139 140.29139
# 12          4  28.23435     plt        19  6.818886   0.138         UNC6958946  17.671082  29.21059
# 13          4 132.95504     wbc        13  7.001567  0.0925         UNC8250192 132.879727 134.54316
# 14          5 137.58353     mcv        13  7.155106  0.0794        UNC10269463 103.337750 139.30484
# 15          7 104.23904    chcm        07 37.856081 < 1e-04        UNC13537149 102.533667 105.60220
# 16          7 104.14897    chcm        13 31.710737 < 1e-04        UNC13534893 102.389926 105.60220
# 17          7 104.23904    chcm        19 22.838725 < 1e-04        UNC13537149 101.775441 105.70587
# 18          7 105.60220     hdw        07 20.409209 < 1e-04 backupUNC070393214 101.661608 105.71390
# 19          7 102.38993     hdw        13 10.313110   2e-04        UNC13501160  99.623551 109.07832
# 20          7 102.76635     hdw        19  8.213306  0.0115        UNC13509119  98.351788 107.37202
# 21          8  12.56355      ch        13  6.933325  0.1075        UNC14177910  10.738292  13.47682
# 22          9 107.57624    chcm        07 14.326926 < 1e-04        JAX00175387 106.563624 108.51902
# 23          9 111.17766    chcm        13  9.722578  0.0011        UNC17122103 107.278677 111.48467
# 24          9 111.43643    chcm        19 11.598809 < 1e-04        UNC17125737 108.973633 111.56281
# 25          9  73.87704     hdw        07  8.689333  0.0043        UNC16632607  69.492256  78.33424
# 26          9 108.09224     hdw        07 13.433244 < 1e-04        UNC17091640 107.582543 108.34486
# 27          9  74.90381     hdw        13  8.100343  0.0153        UNC16645564  71.485485  78.24372
# 28          9 108.06443     hdw        13 15.147539 < 1e-04        UNC17091435 107.582543 108.32187
# 29          9 108.14416     hdw        19 11.725464   1e-04        UNC17092085 107.606558 108.32187
# 30          9 108.94490     mcv        07  9.446472   7e-04        JAX00705053 108.015342 110.91745
# 31          9   6.25458     rbc        07  7.641779  0.0314       UNC090007514   3.369866  17.89850
# 32          9 111.88493     rbc        07  7.966726  0.0184        UNC17134162 110.505368 113.01624
# 33          9  21.49732     rbc        13  7.445345  0.0422        UNC15968208  15.172830  23.20059
# 34         10  65.76825      ch        13  7.214548  0.0664        UNC18123321  58.908569  68.86711
# 35         10  61.98939     mcv        13  6.709482  0.1574        UNC18075499  52.886204  69.25914
# 36         11  88.93354      ch        07  7.466815  0.0413        UNC20082690  82.146142  97.06966
# 37         11 111.29691      ch        07  6.580797  0.1749        UNC20388850 109.296395 112.92104
# 38         11  91.74932      ch        13  7.674637  0.0285        UNC20128412  83.656832  97.14728
# 39         11 111.29691      ch        13  6.895180  0.1145        UNC20388850 108.172117 112.46240
# 40         11  92.20418      ch        19  8.253742  0.0124       UNC111437875  88.921751  97.02702
# 41         11  88.93354     mcv        07  7.391423  0.0466        UNC20082690  87.857447  92.74347
# 42         11 111.29691     mcv        07  8.703195  0.0034        UNC20388850 111.241337 112.42449
# 43         11 111.29691     mcv        13  6.760784  0.1448        UNC20388850  80.937399 111.61020
# 44         11  21.66367     mpm        07  6.900792  0.1067  backupJAX00025163  18.533106  40.83727
# 45         12  77.76394     wbc        13  6.815529  0.1291        UNC21446355  77.023767  83.04458
# 46         13  34.30827     nlr        07  6.596819  0.1879        UNC22392281  20.139385  35.08892
# 47         16  14.78359     plt        07  7.622534  0.0336        UNC26379353  13.479130  16.41683
# 48         16  96.12278     rbc        19  6.932936  0.1151        JAX00072810  93.074097  96.61885
# 49         17  66.17727    chcm        13  7.221415  0.0662        UNC28207384  65.396598  67.40928
# 50         17  32.29756     wbc        07  8.285864    0.01        UNC27781681  31.329496  43.50301
# 51         18  17.06342     hdw        07  6.994611   0.096        UNC28802837  13.462651  24.62421
# 52         18  16.12121     hdw        13  7.638294   0.032        UNC28790404  13.405937  23.77372
# 53         18  10.84216   n.eos        19  6.596502  0.1951        UNC28719154  10.631570  12.40209
# 54         18  15.83413 n.lymph        13  7.057449  0.0887        UNC28786580  13.462651  17.06342
# 55         18  15.83136     rdw        07  7.486349  0.0436       UNC180052940  14.053427  18.91441
# 56         18  17.64595     rdw        19  7.809898  0.0248 backupUNC180056465  13.462651  20.60528
# 57         19  10.16726  n.mono        07  8.377103  0.0083        JAX00470180   7.347888  10.32098
# 58          X 142.88053      ch        07  7.147433  0.0712        UNC31408194 141.241721 153.84331
# 59          X  95.90337  n.neut        13  6.550093  0.1978        UNC31113437  93.901906  99.15429
# 60          X 140.10545     wbc        07  6.618580  0.1893        UNC31386546 136.599742 140.93831

#####


################################################################################
## save ###########################################

write_csv(
  LIFESPAN_PEAK_TABLE,
  'tables/genome_wide_linkage_scans/Lifespan_Peak_Table.csv'
)

write_csv(
  AGING_BIOM_PEAK_TABLE,
  'tables/genome_wide_linkage_scans/Aging_Biomarker_Hema_Peak_Table.csv'
)

write_csv(
  OTHER_HEMA_PEAK_TABLE,
  'tables/genome_wide_linkage_scans/Other_Hema_Peak_Table.csv'
)

#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####

