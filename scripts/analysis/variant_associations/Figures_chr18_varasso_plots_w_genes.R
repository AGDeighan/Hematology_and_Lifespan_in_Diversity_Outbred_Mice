# 2021-01-15



################################################################################
## load libraries ###########################################

options(stringsAsFactors = FALSE)
library(qtl2)
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

protein_coding_genes <- function(
  chr, lb, ub,
  query_genes = QUERY_GENES
){
  GENES <- query_genes(
    chr = chr,
    start = lb,
    end = ub
  )
  
  PROTEIN_CODING_GENES <- GENES %>% 
    filter(mgi_type == 'protein coding gene') %>% 
    filter(substr(tolower(description), 1, 14) != 'predicted gene') %>% 
    filter(substr(tolower(description), 1, 5) != 'riken') %>% 
    filter(substr(tolower(description), 1, 8) != 'microrna') %>% 
    filter(substr(tolower(description), 1, 9) != 'predicted') %>% 
    filter(substr(tolower(description), 1, 13) != 'expressed seq') %>% 
    select(start, stop, strand, Name)
  
  return(PROTEIN_CODING_GENES)
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

save_pdf <- function(plot, ...){
  pdf(...)
  plot(plot)
  dev.off()
}

#####


################################################################################
## Load data etc #############

load('data/processed/qtl_analysis_data/environments/LifespanHema_Environment.RData')
rm(
  K, 
  dataset.LifespanCovSexGenHDW.phenotype, 
  dataset.LifespanCovSexGenNLR.phenotype,
  dataset.LifespanCovSexGenRDW.phenotype,
  dataset.LifespanHemaCovSexGen.phenotype,
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
## Set up #############

# Subset scan
SNP_INFO_TABLE <- ADD_SCANS$snpinfo %>% 
  filter(
    chr == '18',
    pos >= 0,
    pos <= 30
  ) %>% 
  index_snps(
    map = map
  )
ADD_SCANS <- ADD_SCANS$lod[SNP_INFO_TABLE$snp_id, , drop = FALSE]


# Get peak table
PEAK_TABLE <- find_peaks(
  ADD_SCANS,
  SNP_INFO_TABLE,
  drop = 1.5,
  peakdrop = Inf,
  threshold = 2
)
#   lodindex       lodcolumn chr      pos      lod     ci_lo    ci_hi
# 1        1        Lifespan  18 12.57542 4.324162  3.007194 23.37704
# 2        2          rdw_07  18 16.14743 6.021129 10.553241 19.95190
# 3        3          hdw_07  18 17.06348 6.579012 14.120072 20.31640
# 4        4 LifespanCondRDW  18 12.57542 3.747732  3.007194 18.79850
# 5        5 LifespanCondHDW  18 12.57542 3.504044  3.007194 18.79850


# Create a LOD-score table
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
    FAP = sdp_names(sdp)
  )

# Pull genes in region
GENES <- protein_coding_genes(
  chr = '18',
  lb = 0, 
  ub = 28
) %>% 
  mutate(
    PLOT_POS = rep(1:7, length.out = n())
  )

# Find leading founder allele patterns for each trait
LEAD_FAP <- LOD_TABLE %>% 
  filter(
    Trait %in% c('Lifespan', 'rdw_07', 'hdw_07')
  ) %>% 
  group_by(Trait) %>% 
  filter(
    LOD > quantile(LOD, 0.995, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  group_by(
    Trait, FAP
  ) %>% 
  summarise(
    MaxLOD = max(LOD, na.rm = TRUE),
    snp_id = snp_id[which.max(LOD)]
  ) %>% 
  ungroup() %>% 
  arrange(
    factor(
      Trait, levels = c('Lifespan', 'rdw_07', 'hdw_07')
    ),
    desc(MaxLOD)
  )
# # A tibble: 7 x 4
#    Trait    FAP                MaxLOD snp_id     
#    <chr>    <chr>               <dbl> <chr>      
#  1 Lifespan 129:WSB              4.32 rs50243818 
#  2 Lifespan CAST                 4.04 rs251895426
#  3 Lifespan 129:NZO:WSB          3.96 rs30257971 
#  4 Lifespan B6:NOD:CAST          3.81 rs29814803 
#  5 Lifespan AJ:129:NZO:PWK:WSB   3.80 rs252217772
#  6 rdw_07   CAST                 6.02 rs49905069 
#  7 hdw_07   CAST                 6.58 rs229816196

#####


# Creating plots:

LEADING_SDP <- c('129:WSB', 'CAST',     '129:NZO:WSB', 'B6:NOD:CAST')
# SDP_COLORS  <- c('#D046A4', '#2ECC40',  '#B17DC6',     '#179884')
# SDP_COLORS  <- c('magenta', '#2ECC40',  'plum4',     'paleturquoise4')
SDP_COLORS  <- c('#D046A4', '#2ECC40',  'plum4',     'paleturquoise4')
PLOT_LOD_LIMITS <- c(0, 6.7)
LOD_MAJ_BREAKS <- seq(0, 100, 1)
LOD_MIN_BREAKS <- seq(0, 100, 0.2)
POS_LB <- 0
POS_UB <- 30

################################################################################
## Genes ###########################################

PLOT_MB_LIMITS <- c(
  min(c(POS_LB, GENES$start, GENES$stop)),
  max(c(POS_UB, GENES$start, GENES$stop)) + 0.5
)
MB_MAJ_BREAKS <- seq(0, 500, 1)
MB_MIN_BREAKS <- seq(0, 500, 0.25)

GENE_PLOT <- GENES %>% 
  ggplot() +
  theme_minimal() +
  geom_rect(
    aes(
      ymin = start, ymax = stop, 
      xmin = PLOT_POS - 0.49, xmax = PLOT_POS + 0.49
    ),
    alpha = 0.75,
    fill = 'grey50'
  ) +
  geom_text(
    aes(
      x = PLOT_POS,
      y = stop + .02,
      label = Name
    ),
    size = 2.5,
    vjust = 0
  ) +
  scale_x_continuous(
    limits = c(0.3, 7.7),
    breaks = NULL,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    limits = PLOT_MB_LIMITS,
    breaks = MB_MAJ_BREAKS,
    minor_breaks = MB_MIN_BREAKS,
    exp = c(0,0)
  ) +
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    plot.margin = unit(c(0, 0, .1, 0), "cm")
  ) +
  labs(
    x = "Protein Coding Genes"
  ); GENE_PLOT

#####


################################################################################
## RDW plots ###########################################

HAPLO_SI <- c(13.46265, 20.60528)

# Unconditioned
PLOT_TABLE <- LOD_TABLE %>% 
  filter(
    Trait == 'rdw_07'
  ) %>% 
  mutate(
    FAP = factor(FAP, levels = LEADING_SDP),
    LeadSNP = LOD >= max(LOD)
  ) %>% 
  select(
    snp_id, LOD, pos, FAP, LeadSNP
  )
VARIANT_SI <- c(10.553241, 19.95190)


RDW_UNCOND_PLOT <- PLOT_TABLE %>% 
  filter(!is.na(FAP), !LeadSNP) %>% 
  ggplot() +
  theme_minimal() +
  geom_hline(
    yintercept = c(HAPLO_SI[1], HAPLO_SI[2]),
    linetype = 1
  ) +
  geom_hline(
    yintercept = c(VARIANT_SI[1], VARIANT_SI[2]),
    linetype = 3
  ) +
  geom_point(
    data = PLOT_TABLE %>% 
      filter(is.na(FAP)),
    aes(x = LOD, y = pos),
    color = 'grey90', alpha = 1, shape = 20
  ) +
  geom_point(
    aes(x = LOD, y = pos, col = FAP), 
    alpha = .85, shape = 20
  ) +
  geom_point(
    data = PLOT_TABLE %>% 
      filter(LeadSNP),
    aes(x = LOD, y = pos, fill = FAP), 
    alpha = 1, shape = 23, color = 'black', size = 2
  ) +
  scale_color_manual(
    values = SDP_COLORS,
    drop = FALSE
  ) + 
  scale_fill_manual(
    values = SDP_COLORS,
    na.value = 'grey90',
    drop = FALSE,
    guide = FALSE
  ) + 
  scale_y_continuous(
    limits = PLOT_MB_LIMITS,
    breaks = MB_MAJ_BREAKS,
    minor_breaks = MB_MIN_BREAKS,
    exp = c(0,0)
  ) +
  scale_x_continuous(
    limits = PLOT_LOD_LIMITS,
    breaks = LOD_MAJ_BREAKS,
    minor_breaks = LOD_MIN_BREAKS,
    exp = c(0,0)
  ) +
  theme(
    legend.position = 'top',
    plot.margin = unit(c(0, 0, .1, .11), "cm")
  ) +
  labs(
    x = 'Association with RDW (LOD)',
    y = 'Position on Chromosome 18 (Mbp)',
    color = 'Minor Allele:  '
  ); RDW_UNCOND_PLOT

#####


################################################################################
## HDW plots ###########################################

HAPLO_SI <- c(13.40594, 24.62421)

# Unconditioned
PLOT_TABLE <- LOD_TABLE %>% 
  filter(
    Trait == 'hdw_07'
  ) %>% 
  mutate(
    FAP = factor(FAP, levels = LEADING_SDP),
    LeadSNP = LOD >= max(LOD)
  ) %>% 
  select(
    snp_id, LOD, pos, FAP, LeadSNP
  )
VARIANT_SI <- c(14.120072, 20.31640)


HDW_UNCOND_PLOT <- PLOT_TABLE %>% 
  filter(!is.na(FAP)) %>% 
  ggplot() +
  theme_minimal() +
  geom_hline(
    yintercept = c(HAPLO_SI[1], HAPLO_SI[2]),
    linetype = 1
  ) +
  geom_hline(
    yintercept = c(VARIANT_SI[1], VARIANT_SI[2]),
    linetype = 3
  ) +
  geom_point(
    data = PLOT_TABLE %>% 
      filter(is.na(FAP)),
    aes(x = LOD, y = pos),
    color = 'grey90', alpha = 1, shape = 20
  ) +
  geom_point(
    aes(x = LOD, y = pos, col = FAP), 
    alpha = .85, shape = 20
  ) +
  geom_point(
    data = PLOT_TABLE %>% 
      filter(LeadSNP),
    aes(x = LOD, y = pos, fill = FAP), 
    alpha = 1, shape = 23, color = 'black', size = 2
  ) +
  scale_color_manual(
    values = SDP_COLORS,
    drop = FALSE
  ) + 
  scale_fill_manual(
    values = SDP_COLORS,
    na.value = 'grey90',
    drop = FALSE,
    guide = FALSE
  ) + 
  scale_y_continuous(
    limits = PLOT_MB_LIMITS,
    breaks = MB_MAJ_BREAKS,
    minor_breaks = MB_MIN_BREAKS,
    exp = c(0,0)
  ) +
  scale_x_continuous(
    limits = PLOT_LOD_LIMITS,
    breaks = LOD_MAJ_BREAKS,
    minor_breaks = LOD_MIN_BREAKS,
    exp = c(0,0)
  ) +
  theme(
    legend.position = 'top',
    plot.margin = unit(c(0, 0, .1, .11), "cm")
  ) +
  labs(
    x = 'Association with HDW (LOD)',
    y = 'Position on Chromosome 18 (Mbp)',
    color = 'Minor Allele:  '
  ); HDW_UNCOND_PLOT


#####


################################################################################
## Lifespan plots ###########################################

HAPLO_SI <- c(8.523683,  14.39734)

# Unconditioned
PLOT_TABLE <- LOD_TABLE %>% 
  filter(
    Trait == 'Lifespan'
  ) %>% 
  mutate(
    FAP = factor(FAP, levels = LEADING_SDP),
    LeadSNP = LOD >= max(LOD)
  ) %>% 
  select(
    snp_id, LOD, pos, FAP, LeadSNP
  )
VARIANT_SI <- c(3.007194, 23.37704)


LS_UNCOND_PLOT <- PLOT_TABLE %>% 
  filter(!is.na(FAP), !LeadSNP) %>% 
  ggplot() +
  theme_minimal() +
  geom_hline(
    yintercept = c(HAPLO_SI[1], HAPLO_SI[2]),
    linetype = 1
  ) +
  geom_hline(
    yintercept = c(VARIANT_SI[1], VARIANT_SI[2]),
    linetype = 3
  ) +
  geom_point(
    data = PLOT_TABLE %>% 
      filter(is.na(FAP)),
    aes(x = LOD, y = pos),
    color = 'grey90', alpha = 1, shape = 20
  ) +
  geom_point(
    aes(x = LOD, y = pos, col = FAP), 
    alpha = .85, shape = 20
  ) +
  geom_point(
    data = PLOT_TABLE %>% 
      filter(LeadSNP),
    aes(x = LOD, y = pos, fill = FAP), 
    alpha = 1, shape = 23, color = 'black', size = 2
  ) +
  scale_color_manual(
    values = SDP_COLORS,
    drop = FALSE
  ) + 
  scale_fill_manual(
    values = SDP_COLORS,
    na.value = 'grey90',
    drop = FALSE,
    guide = FALSE
  ) + 
  scale_y_continuous(
    limits = PLOT_MB_LIMITS,
    breaks = MB_MAJ_BREAKS,
    minor_breaks = MB_MIN_BREAKS,
    exp = c(0,0)
  ) +
  scale_x_continuous(
    limits = PLOT_LOD_LIMITS,
    breaks = LOD_MAJ_BREAKS,
    minor_breaks = LOD_MIN_BREAKS,
    exp = c(0,0)
  ) +
  theme(
    legend.position = 'top',
    plot.margin = unit(c(0, 0, .1, .11), "cm")
  ) +
  labs(
    x = 'Association with Lifespan (LOD)',
    y = 'Position on Chromosome 18 (Mbp)',
    color = 'Minor Allele:  '
  ); LS_UNCOND_PLOT

#####

LEGEND <- cowplot::plot_grid(
  cowplot::get_legend(RDW_UNCOND_PLOT), NULL, nrow = 1, 
  rel_widths = c(1,2/5)
)

################################################################################
## Lifespan-RDW plots ###########################################

# Unconditioned
cowplot::plot_grid(
  LEGEND, 
  cowplot::plot_grid(
    LS_UNCOND_PLOT + 
      theme(legend.position = 'none'),
    GENE_PLOT,
    RDW_UNCOND_PLOT + 
      scale_x_reverse(
        limits = rev(PLOT_LOD_LIMITS),
        breaks = LOD_MAJ_BREAKS,
        minor_breaks = LOD_MIN_BREAKS,
        exp = c(0,0)
      ) +
      theme(
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = 'none',
        plot.margin = unit(c(0, 0.1, .1, .11), "cm")
      ), 
    nrow = 1, align = 'h', rel_widths = c(100, 85, 95)
  ),
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
) %>% save_pdf(
  'figures/variant_associations/LifespanRDW_unconditioned_variant_scan.pdf',
  width = 9, height = 6
)


#####



################################################################################
## Lifespan-HDW plots ###########################################

# Unconditioned
cowplot::plot_grid(
  LEGEND, 
  cowplot::plot_grid(
    LS_UNCOND_PLOT + 
      theme(legend.position = 'none'),
    GENE_PLOT,
    HDW_UNCOND_PLOT + 
      scale_x_reverse(
        limits = rev(PLOT_LOD_LIMITS),
        breaks = LOD_MAJ_BREAKS,
        minor_breaks = LOD_MIN_BREAKS,
        exp = c(0,0)
      ) +
      theme(
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = 'none',
        plot.margin = unit(c(0, 0.1, .1, .11), "cm")
      ), 
    nrow = 1, align = 'h', rel_widths = c(100, 85, 95)
  ),
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
) %>% save_pdf(
  'figures/variant_associations/LifespanHDW_unconditioned_variant_scan.pdf',
  width = 9, height = 6
)


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

