# 2021-06-22
################################################################################
#
#   This script performs the BLUP scans across chromosome 18 for lifespan, RDW,
#   and HDW to estimate the haplotype effects on these traits
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


################################################################################
## Run scans ###########################################

LS_BLUP_SCAN_CHR18 <- scan1blup(
  genoprobs = genoprobs[,'18'],
  pheno = dataset.LifespanHemaCovSexGen.phenotype$pheno[,'Lifespan', drop = FALSE],
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  cores = 6,
  quiet = FALSE
)

RDW_BLUP_SCAN_CHR18 <- scan1blup(
  genoprobs = genoprobs[,'18'],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == '7 Months') %>% 
    select(MouseID, rdw) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  cores = 6,
  quiet = FALSE
)

HDW_BLUP_SCAN_CHR18 <- scan1blup(
  genoprobs = genoprobs[,'18'],
  pheno = PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == '7 Months') %>% 
    select(MouseID, hdw) %>% 
    column_to_rownames('MouseID') %>%
    data.matrix(),
  kinship = K[['18']],
  addcovar = dataset.LifespanHemaCovSexGen.phenotype$covar,
  cores = 6,
  quiet = FALSE
)

#####


################################################################################
## Save results of scan ###########################################

saveRDS(
  LS_BLUP_SCAN_CHR18,
  'results/haplotype_effects/Chr18_BLUPs_Lifespan.rds'
)

saveRDS(
  RDW_BLUP_SCAN_CHR18,
  'results/haplotype_effects/Chr18_BLUPs_RDW07.rds'
)

saveRDS(
  HDW_BLUP_SCAN_CHR18,
  'results/haplotype_effects/Chr18_BLUPs_HDW07.rds'
)

#####


################################################################################
## Plot ###########################################

pdf(
  'figures/haplotype_effects/full_chr18_BLUP_scans_LS_RDW_HDW.pdf',
  width = 7, height = 9
)
cowplot::plot_grid(
  qtl2ggplot::ggplot_coef(
    LS_BLUP_SCAN_CHR18,
    map = map,
    ylim = c(-4, 4)
  ) +
    theme_minimal(
      base_size = 9
    ) +
    geom_vline(
      xintercept = c(8.523683, 24.62421),
      linetype = 3
    ) +
    geom_vline(
      xintercept = c(10.77744),
      linetype = 2
    ) +
    scale_color_manual(
      values = CCcolors
    ) +
    scale_x_continuous(
      expand = c(0, 0)
    ) +
    theme(
      legend.position = 'top',
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")
    ) +
    labs(
      title = 'Estimated allele effects (BLUP) on lifespan, RDW, and HDW',
      color = 'Allele',
      x = NULL,
      y = 'Effect on lifespan (months)'
    ),
  
  qtl2ggplot::ggplot_coef(
    RDW_BLUP_SCAN_CHR18,
    map = map,
    ylim = c(-1, 1)
  ) +
    theme_minimal(
      base_size = 9
    ) +
    geom_vline(
      xintercept = c(8.523683, 24.62421),
      linetype = 3
    ) +
    geom_vline(
      xintercept = c(10.77744),
      linetype = 2
    ) +
    scale_color_manual(
      values = CCcolors
    ) +
    scale_x_continuous(
      expand = c(0, 0)
    ) +
    theme(
      legend.position = 'none',
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = unit(c(2.1,0.1,2.1,0.1), "lines")
    ) +
    labs(
      x = NULL,
      y = 'Effect on 7-month RDW (%)'
    ),
  
  qtl2ggplot::ggplot_coef(
    HDW_BLUP_SCAN_CHR18,
    map = map,
    ylim = c(-0.4, 0.4)
  ) +
    theme_minimal(
      base_size = 9
    ) +
    geom_vline(
      xintercept = c(8.523683, 24.62421),
      linetype = 3
    ) +
    geom_vline(
      xintercept = c(10.77744),
      linetype = 2
    ) +
    scale_color_manual(
      values = CCcolors
    ) +
    scale_x_continuous(
      expand = c(0, 0)
    ) +
    theme(
      legend.position = 'none',
      plot.caption = element_text(hjust = 0),
      plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")
    ) +
    labs(
      x = 'Position on chromosome 18 (Mb)',
      y = 'Effect on 7-month HDW (%)',
      caption = 'The vertical dashed line indicates the location of the lifespan peak, and the dotted lines indicate the upper and lower bounds of the merged LOD support\nintervals for lifespan, RDW, and HDW'
    ),
  align = 'v',
  axis = 'lr',
  ncol = 1,
  rel_heights = c(10, 9.6, 8.75)
)
dev.off()

#####


################################################################################
##  ###########################################



#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####
