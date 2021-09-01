# 2021-08-30



################################################################################
## load libraries ###########################################

options(stringsAsFactors = FALSE)
library(qtl2)
library(ppcor)
library(tidyverse)

#####


################################################################################
#### Helper functions ###################

sdp_names <- function(sdps, sep = ':'){
  
  SPD_PATTERN <- qtl2::invert_sdp(sdps, 8) == 3
  
  get_names <- function(r){
    names <- paste(names(qtl2::CCcolors)[SPD_PATTERN[r,]], collapse = sep)
    return(names)
  }
  
  NAMES <- unlist(lapply(1:nrow(SPD_PATTERN), get_names))
  
  return(NAMES)
}

get_minorallele_probs <- function(id, ...){
  x <- genoprob_to_snpprob(...)
  if(length(id) == 1){
    x <- x[[1]][,,id]
    x <- x[,'B', drop = FALSE]
    colnames(x) <- id
  } else{
    id <- id[id %in% dimnames(x[[1]])[[3]]]
    x <- x[[1]][,,id]
    x <- x[,'B',]
  }
  return(x)
}

#####


################################################################################
## Load data etc #############

load(
  'data/processed/qtl_analysis_data/environments/LifespanHema_Environment.RData'
)
rm(
  dataset.LifespanCovSexGenHDW.phenotype, 
  dataset.LifespanCovSexGenNLR.phenotype,
  dataset.LifespanCovSexGenRDW.phenotype,
  markers
)

QUERY_VAR <- create_variant_query_func(dbfile = 'data/raw/genetic/cc_variants.sqlite')
QUERY_GENES <- create_gene_query_func(dbfile = 'data/raw/genetic/mouse_genes_mgi.sqlite')

ADD_SCANS <- readRDS('results/variant_associations/LifespanHema_GenomeWideVariantScans.rds')
PERM_SCANS <- readRDS('results/variant_associations/LifespanHema_10000Perms.rds')


read.csv('tables/genome_wide_linkage_scans/Lifespan_Peak_Table.csv')
#   Chromosome  Position    Trait      LOD LODDropRDW LODDropHDW  LODDropNLR PValue  PeakMarker       lbSI      ubSI
# 1          2 128.76415 Lifespan 7.220318 -0.1324881  0.5770305 -0.04400962 0.0648  UNC3946248 128.372768 129.36445
# 2          6  25.22833 Lifespan 7.018060  0.3318909  0.3219315  0.54619685 0.0904 UNC10749252  23.581197  27.65261
# 3         18  10.77744 Lifespan 6.591569  1.3006534  1.9110734 -0.25081632 0.1853 UNC28718524   8.523683  14.39734

read.csv('tables/genome_wide_linkage_scans/Aging_Biomarker_Hema_Peak_Table.csv')
#    Chromosome  Position Trait Timepoint       LOD  PValue         PeakMarker      lbSI      ubSI
# 1           7 105.60220   hdw         7 20.409209 < 1e-04 backupUNC070393214 101.66161 105.71390
# 2           7 102.38993   hdw        13 10.313110   2e-04        UNC13501160  99.62355 109.07832
# 3           7 102.76635   hdw        19  8.213306  0.0115        UNC13509119  98.35179 107.37202
# 4           9  73.87704   hdw         7  8.689333  0.0043        UNC16632607  69.49226  78.33424
# 5           9 108.09224   hdw         7 13.433244 < 1e-04        UNC17091640 107.58254 108.34486
# 6           9  74.90381   hdw        13  8.100343  0.0153        UNC16645564  71.48548  78.24372
# 7           9 108.06443   hdw        13 15.147539 < 1e-04        UNC17091435 107.58254 108.32187
# 8           9 108.14416   hdw        19 11.725464   1e-04        UNC17092085 107.60656 108.32187
# 9          13  34.30827   nlr         7  6.596819  0.1879        UNC22392281  20.13939  35.08892
# 10         18  17.06342   hdw         7  6.994611   0.096        UNC28802837  13.46265  24.62421
# 11         18  16.12121   hdw        13  7.638294   0.032        UNC28790404  13.40594  23.77372
# 12         18  15.83136   rdw         7  7.486349  0.0436       UNC180052940  14.05343  18.91441
# 13         18  17.64595   rdw        19  7.809898  0.0248 backupUNC180056465  13.46265  20.60528


#####


################################################################################
## subset variant scans to chromosome 18 #############

SNP_INFO_TABLE <- ADD_SCANS$snpinfo %>% 
  filter(
    chr == '18'
  ) %>% 
  index_snps(
    map = map
  )

ADD_SCANS <- ADD_SCANS$lod[SNP_INFO_TABLE$snp_id, , drop = FALSE]


#####


################################################################################
## variant scan peak table #############

PEAK_TABLE <- find_peaks(
  scan1_output = ADD_SCANS,
  map = SNP_INFO_TABLE,
  threshold = 2,
  peakdrop = Inf,
  drop = 1.5
) 
#   lodindex       lodcolumn chr      pos      lod     ci_lo    ci_hi
# 1        1        Lifespan  18 12.57542 4.324162  3.007194 23.37704
# 2        2          rdw_07  18 16.14743 6.021129 10.553241 19.95190
# 3        3          hdw_07  18 17.06348 6.579012 14.120072 20.31640
# 4        4 LifespanCondRDW  18 12.57542 3.747732  3.007194 18.79850
# 5        5 LifespanCondHDW  18 12.57542 3.504044  3.007194 18.79850

#####


################################################################################
## Subset to region of interest #############


SNP_INFO_TABLE <- SNP_INFO_TABLE %>% 
  filter(
    pos > 0,
    pos <= 30
  )

ADD_SCANS <- ADD_SCANS[SNP_INFO_TABLE$snp_id, , drop = FALSE]

#####


################################################################################
## Create table with LOD scores for each trait at each variant #############

LOD_TABLE <- ADD_SCANS %>% 
  data.frame() %>% 
  rownames_to_column('snp_id')%>% 
  pivot_longer(
    cols = Lifespan:LifespanCondHDW,
    names_to = 'Trait',
    values_to = 'LOD'
  ) %>% 
  left_join(
    SNP_INFO_TABLE %>% 
      select(snp_id, pos, sdp),
    by = 'snp_id'
  ) %>% 
  mutate(
    FAP = sdp_names(sdp),
    Trait = factor(
      Trait,
      levels = c('Lifespan', 'LifespanCondRDW', 'LifespanCondHDW', 'rdw_07', 'hdw_07')
    )
  )


#####


################################################################################
## Get lead variants for each FAP #############

LEAD_FAP <- LOD_TABLE %>% 
  group_by(
    Trait, sdp, FAP
  ) %>% 
  summarise(
    MaxLOD = max(LOD, na.rm = TRUE),
    LeadVar = snp_id[which.max(LOD)]
  ) %>% 
  ungroup() %>% 
  arrange(
    Trait, desc(MaxLOD)
  )

LEAD_FAP %>% 
  filter(
    Trait == 'Lifespan'
  ) %>% 
  data.frame() %>% 
  head(6)
#      Trait sdp                FAP   MaxLOD     LeadVar
# 1 Lifespan 132            129:WSB 4.324162  rs50243818
# 2 Lifespan  32               CAST 4.035748 rs251895426
# 3 Lifespan 148        129:NZO:WSB 3.957871  rs30257971
# 4 Lifespan  42        B6:NOD:CAST 3.807898  rs29814803
# 5 Lifespan 213 AJ:129:NZO:PWK:WSB 3.799089 rs252217772
# 6 Lifespan 212    129:NZO:PWK:WSB 3.658353   rs3704039

LEAD_FAP %>% 
  filter(
    Trait == 'rdw_07'
  ) %>% 
  data.frame() %>% 
  head(6)
#    Trait sdp          FAP   MaxLOD     LeadVar
# 1 rdw_07  32         CAST 6.021129  rs49905069
# 2 rdw_07  34      B6:CAST 4.967919 rs256828598
# 3 rdw_07  44 129:NOD:CAST 3.901116 rs387022687
# 4 rdw_07  50  B6:NZO:CAST 3.881201 rs253707890
# 5 rdw_07  36     129:CAST 3.628944 rs220453347
# 6 rdw_07  33      AJ:CAST 3.622094  rs51901914

LEAD_FAP %>% 
  filter(
    Trait == 'hdw_07'
  ) %>% 
  data.frame() %>% 
  head(6)
#    Trait sdp             FAP   MaxLOD     LeadVar
# 1 hdw_07  32            CAST 6.579012 rs229816196
# 2 hdw_07  34         B6:CAST 4.729331 rs256828598
# 3 hdw_07  96        CAST:PWK 4.566047 rs221451784
# 4 hdw_07  50     B6:NZO:CAST 4.172508 rs253707890
# 5 hdw_07  36        129:CAST 3.799767 rs236803032
# 6 hdw_07 114 B6:NZO:CAST:PWK 3.572731  rs31124698

#####


################################################################################
## Get SNP-probs #############

SNP_PROBS <- get_minorallele_probs(
  id = unique(LOD_TABLE$snp_id), 
  # id = LEAD_FAP$LeadVar[LEAD_FAP$Trait == 'Lifespan'], 
  genoprobs = genoprobs,
  snpinfo = SNP_INFO_TABLE
)

#####


################################################################################
## Correlation with lead 129:WSB SNP  #############

COR_DISTANCE_129WSB <- SNP_PROBS %>% 
  as_tibble() %>% 
  mutate(
    MouseID = rownames(SNP_PROBS)
  ) %>% 
  select(
    MouseID, rs50243818,
    everything()
  ) %>% 
  pivot_longer(
    cols = -c(MouseID, rs50243818),
    names_to = 'snp_id',
    values_to = 'minor_prob'
  ) %>% 
  group_by(
    snp_id
  ) %>% 
  summarise(
    Correlation = abs(cor(rs50243818, minor_prob))
  ) %>% 
  left_join(
    LOD_TABLE %>% 
      filter(
        Trait == 'Lifespan'
      ) %>% 
      select(
        snp_id, pos, FAP
      ),
    by = 'snp_id'
  ) %>% 
  mutate(
    Distance = 12.57542 - pos
  ) %>% 
  select(
    -pos
  )

CORPLOT_129WSB <- COR_DISTANCE_129WSB %>% 
  filter(
    !FAP %in% c('129:CAST:WSB', '129:WSB', 'CAST', '129', 'WSB', 'NOD')
  ) %>% 
  ggplot() +
  theme_minimal() +
  geom_vline(xintercept = 0) +
  geom_point(
    aes(x = Distance, y = Correlation),
    color = 'gray60',
    size = 0.5
  ) +
  geom_point(
    data = COR_DISTANCE_129WSB %>% 
      filter(
        FAP == '129:CAST:WSB'
      ),
    aes(x = Distance, y = Correlation),
    color = 'red',
    size = 0.5
  ) +
  geom_point(
    data = COR_DISTANCE_129WSB %>% 
      filter(
        FAP == '129:WSB'
      ),
    aes(x = Distance, y = Correlation),
    color = 'black',
    size = 0.5
  ) +
  geom_point(
    data = COR_DISTANCE_129WSB %>% 
      filter(
        FAP == 'CAST'
      ),
    aes(x = Distance, y = Correlation),
    color = '#2ECC40',
    size = 0.5
  ) +
  geom_point(
    data = COR_DISTANCE_129WSB %>% 
      filter(
        FAP == '129'
      ),
    aes(x = Distance, y = Correlation),
    color = '#F08080',
    size = 0.5
  ) +
  geom_point(
    data = COR_DISTANCE_129WSB %>% 
      filter(
        FAP == 'WSB'
      ),
    aes(x = Distance, y = Correlation),
    color = '#B10DC9',
    size = 0.5
  ) +
  geom_point(
    data = COR_DISTANCE_129WSB %>% 
      filter(
        FAP == 'NOD'
      ),
    aes(x = Distance, y = Correlation),
    color = '#0064C9',
    size = 0.5
  ) +
  labs(
    title = 'Correlation with lead 129/WSB Variant',
    subtitle = 'rs50243818',
    x = 'Distance (Mb) from lead variant',
    caption = 'Black = 129/WSB, Red = 120/WSB/CAST, Green = CAST, Purple = WSB, Peach = 129, Blue = NOD'
  ); CORPLOT_129WSB

pdf(
  'figures/variant_associations/Variant_correlation_and_distance_plot_129-WSB.pdf',
  width = 7, height = 5
)
plot(CORPLOT_129WSB)
dev.off()

#####


################################################################################
## Correlation with lead CAST SNP  #############

COR_DISTANCE_CAST <- SNP_PROBS %>% 
  as_tibble() %>% 
  mutate(
    MouseID = rownames(SNP_PROBS)
  ) %>% 
  select(
    MouseID, rs251895426,
    everything()
  ) %>% 
  pivot_longer(
    cols = -c(MouseID, rs251895426),
    names_to = 'snp_id',
    values_to = 'minor_prob'
  ) %>% 
  group_by(
    snp_id
  ) %>% 
  summarise(
    Correlation = abs(cor(rs251895426, minor_prob))
  ) %>% 
  left_join(
    LOD_TABLE %>% 
      filter(
        Trait == 'Lifespan'
      ) %>% 
      select(
        snp_id, pos, FAP
      ),
    by = 'snp_id'
  ) %>% 
  mutate(
    Distance = 17.71143 - pos
  ) %>% 
  select(
    -pos
  )

CORPLOT_CAST <- COR_DISTANCE_CAST %>% 
  filter(
    !FAP %in% c('129:CAST:WSB', '129:WSB', 'CAST', '129', 'WSB', 'NOD')
  ) %>% 
  ggplot() +
  theme_minimal() +
  geom_vline(xintercept = 0) +
  geom_point(
    aes(x = Distance, y = Correlation),
    color = 'gray60',
    size = 0.5
  ) +
  geom_point(
    data = COR_DISTANCE_CAST %>% 
      filter(
        FAP == '129:CAST:WSB'
      ),
    aes(x = Distance, y = Correlation),
    color = 'red',
    size = 0.5
  ) +
  geom_point(
    data = COR_DISTANCE_CAST %>% 
      filter(
        FAP == '129:WSB'
      ),
    aes(x = Distance, y = Correlation),
    color = 'black',
    size = 0.5
  ) +
  geom_point(
    data = COR_DISTANCE_CAST %>% 
      filter(
        FAP == 'CAST'
      ),
    aes(x = Distance, y = Correlation),
    color = '#2ECC40',
    size = 0.5
  ) +
  geom_point(
    data = COR_DISTANCE_CAST %>% 
      filter(
        FAP == '129'
      ),
    aes(x = Distance, y = Correlation),
    color = '#F08080',
    size = 0.5
  ) +
  geom_point(
    data = COR_DISTANCE_CAST %>% 
      filter(
        FAP == 'WSB'
      ),
    aes(x = Distance, y = Correlation),
    color = '#B10DC9',
    size = 0.5
  ) +
  geom_point(
    data = COR_DISTANCE_CAST %>% 
      filter(
        FAP == 'NOD'
      ),
    aes(x = Distance, y = Correlation),
    color = '#0064C9',
    size = 0.5
  ) +
  labs(
    title = 'Correlation with lead CAST Variant',
    subtitle = 'rs251895426',
    x = 'Distance (Mb) from lead variant',
    caption = 'Black = 129/WSB, Red = 120/WSB/CAST, Green = CAST, Purple = WSB, Peach = 129, Blue = NOD'
  ); CORPLOT_CAST

pdf(
  'figures/variant_associations/Variant_correlation_and_distance_plot_CAST.pdf',
  width = 7, height = 5
)
plot(CORPLOT_CAST)
dev.off()

#####


################################################################################
## save ###########################################


#####



################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####

