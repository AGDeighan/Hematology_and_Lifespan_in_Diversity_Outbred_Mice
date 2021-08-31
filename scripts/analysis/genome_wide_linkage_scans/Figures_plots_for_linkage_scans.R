# 2021-06-16
################################################################################
#
#   This script creates plots summarizing the results of the genome-wide
#   linkage (haplotype) scans for lifespan and the blood traits
# 
#   Author:  Andrew Deighan
#   E-mail: andrew.deighan@jax.org
#
################################################################################

options(na.action = 'na.exclude')
options(stringsAsFactors = FALSE)
DEFAULT_PLOT_PAR <- par()

#####


################################################################################
## libraries etc #########################################################

library(qtl2)
library(tidyverse)
library(qtl2ggplot)
library(cowplot)

#####


################################################################################
## helper functions #########################################################

scale_0_1 <- function(x){
  MIN <- min(x, na.rm = TRUE)
  MAX <- max(x, na.rm = TRUE)
  return((x - MIN) / (MAX - MIN))
}

create_pheno_lod_table <- function(
  add_scan, map_df, 
  id_string = NULL, 
  pos_precision = .005,
  lod_precision = 0.1,
  lod_limits = c(5,10),
  pos_offset = 0.5
){
  
  map_df <- map_df %>% 
    mutate(chr = factor(chr, levels = unique(chr))) %>% 
    group_by(chr) %>% 
    mutate(
      SIZE = diff(range(pos, na.rm = TRUE)),
      pos = scale_0_1(pos)
    ) %>% 
    ungroup() %>% 
    mutate(
      SIZE = scale_0_1(SIZE) + 0.5,
      SIZE = SIZE / 5 + 0.8,
      pos = round(pos/(pos_precision/SIZE))*(pos_precision/SIZE),
      pos = pos + as.numeric(chr)
    )
  
  PHENO_LOD <- add_scan %>% 
    data.frame() %>% 
    rownames_to_column('marker') %>% 
    select(marker, everything()) 
  
  if(!is.null(id_string)){
    PHENO_LOD <- PHENO_LOD %>% 
      select(marker, contains(id_string)) 
    names(PHENO_LOD) <- gsub(id_string, '', names(PHENO_LOD))
  }
  
  FEATURES <- names(PHENO_LOD)[-1]
  
  PHENO_LOD <- PHENO_LOD %>% 
    select(
      marker, 
      all_of(FEATURES)
    ) %>% 
    pivot_longer(
      cols = all_of(FEATURES),
      values_to = 'LOD',
      names_to = 'Phenotype'
    ) %>% 
    mutate(
      Phenotype = factor(Phenotype, levels = FEATURES)
    ) %>% 
    left_join(
      map_df, 
      by = 'marker'
    ) %>% 
    as_tibble()
  
  PHENO_LOD <- PHENO_LOD %>% 
    mutate(pos = as.character(pos)) %>% 
    group_by(Phenotype, chr, pos) %>% 
    summarise(LOD = mean(LOD)) %>% 
    ungroup() %>% 
    mutate(pos = as.numeric(pos)) %>% 
    arrange(chr, pos) %>% 
    group_by(Phenotype, chr) %>% 
    mutate(
      LeadLOD1 = lead(LOD, 1), 
      LeadLOD2 = lead(LOD, 2), 
      LeadLOD3 = lead(LOD, 3), 
      LagLOD1 = lag(LOD, 1), 
      LagLOD2 = lag(LOD, 2), 
      LagLOD3 = lag(LOD, 3)
    ) %>% 
    ungroup() %>% 
    rowwise() %>% 
    mutate(
      LOD = weighted.mean(
        c(LeadLOD3,LeadLOD2,LeadLOD1, LOD, LagLOD1,LagLOD2,LagLOD3), 
        c(0.1, 0.15, 0.25, 1, 0.25, 0.15, 0.1), 
        na.rm = TRUE
      )
    ) %>% 
    ungroup() %>% 
    mutate(LOD = round(LOD/lod_precision)*lod_precision) %>% 
    mutate(LOD = ifelse(LOD > lod_limits[2], lod_limits[2], LOD)) %>% 
    mutate(LOD = ifelse(LOD < lod_limits[1], NA, LOD)) %>% 
    mutate(pos = pos - pos_offset) %>% 
    select(-starts_with('Lead'), -starts_with('Lag'))
  
  return(PHENO_LOD)
}

lod_heat_plot <- function(
  pheno_lod_table,
  xlim = c(0.5, 20.5), x_breaks = c(1:20), x_labels = c(1:19,'X'),
  ylim, y_breaks, y_labels,
  colors2 = c('darkorange', 'firebrick4'),
  nacolor = '#ffffe5',
  color.midpoint = 7.5,
  divider_color = 'gray50'
){
  pheno_lod_table <- pheno_lod_table %>% 
    mutate(MIN_X = lag(pos))
  
  PLOT <- pheno_lod_table %>% 
    filter(!is.na(MIN_X)) %>% 
    ggplot(aes(color = LOD, fill = LOD)) +
    geom_rect(
      aes(xmin = MIN_X, xmax = pos,
          ymin = as.numeric(Phenotype) - 0.5, ymax = as.numeric(Phenotype) + 0.5
      )
    ) +
    scale_x_continuous(
      limits = xlim,
      expand = c(0,0),
      breaks = x_breaks,
      labels = x_labels
    ) +
    scale_y_continuous(
      limits = ylim,
      expand = c(0,0),
      breaks = y_breaks,
      labels = y_labels
    ) +
    scale_colour_gradient2(
      low = colors2[1], mid = colorRampPalette(colors2)(3)[2],  high = colors2[2], 
      midpoint = color.midpoint,
      breaks = c(5, 6, 7, 8, 9, 10),
      limits = c(5,10),
      labels = c(5, 6, 7, 8, 9, '> 10'),
      na.value = nacolor
    ) +
    scale_fill_gradient2(
      low = colors2[1], mid = colorRampPalette(colors2)(3)[2], high = colors2[2], 
      midpoint = color.midpoint,
      breaks = c(5, 6, 7, 8, 9, 10),
      limits = c(5,10),
      labels = c(5, 6, 7, 8, 9, '> 10'),
      na.value = nacolor, 
      guide = FALSE
    ) +
    theme(
      panel.background = element_rect(
        fill = nacolor,
        colour = nacolor,
        size = 0.5, linetype = "solid"
      ),
      axis.ticks = element_blank()
    ) + 
    labs(
      x = 'Chromosome'
    )
  
  for(i in x_breaks[-1]-1){
    PLOT <- PLOT +
      geom_vline(xintercept = i + 0.5, color = divider_color, size = .5)
  }
  return(PLOT)
}

#####


################################################################################
## load data etc #########################################################

load('data/processed/phenotypic/CBC_AnalysisEnvironment.RData')

MAP <- readRDS(
  'data/processed/genetic/PhysicalMap_List.rds'
)
MAP_DF <- readRDS(
  'data/processed/genetic/PhysicalMap_Dataframe.rds'
)
ADD_SCANS <- readRDS(
  'results/genome_wide_linkage_scans/LifespanHema_HaploScans.rds'
)
PERM_SCANS <- readRDS(
  'results/genome_wide_linkage_scans/LifespanHema_10000Perms.rds'
)

LIFESPAN_PEAK_TABLE <- read_csv(
  'tables/genome_wide_linkage_scans/Lifespan_Peak_Table.csv'
)
AGING_BIOM_PEAK_TABLE <- read_csv(
  'tables/genome_wide_linkage_scans/Aging_Biomarker_Hema_Peak_Table.csv'
)
OTHER_HEMA_PEAK_TABLE <- read_csv(
  'tables/genome_wide_linkage_scans/Other_Hema_Peak_Table.csv'
)

#####


################################################################################
##  #################

ALPHA <- c(0.05, 0.2)
PERM_THRESH <- summary_scan1perm(PERM_SCANS, alpha = ALPHA)

#####


################################################################################
## Plot of peaks #################

PEAK_PLOT <- qtl2ggplot::ggplot_peaks(
  peaks = find_peaks(
    subset_scan1(ADD_SCANS, lodcolumn = colnames(PERM_SCANS)),
    map = MAP,
    threshold = PERM_THRESH['0.2',],
    peakdrop = 5,
    drop = 1.5
  ) %>% 
    mutate(
      lodcolumn = factor(lodcolumn, levels = unique(lodcolumn))
    ),
  map = MAP
)
pdf(
  'figures/genome_wide_linkage_scans/Peak_Plot.pdf',
  width = 7, height = 7
)
plot(PEAK_PLOT)
dev.off()

rm(PEAK_PLOT)


#####


################################################################################
## CBC trait linkage scan plots #################

for(PHENO in CBCphenos){
  # Create initial addscan plots
  TRAITS <- paste(PHENO, c('07', '13', '19'), sep = '_')
  PT <- PERM_THRESH[, TRAITS]
  SCAN1 <- ADD_SCANS[, TRAITS]
  MAXLOD <- maxlod(SCAN1)
  MAXLOD <- ceiling(MAXLOD/2.5)*2.5
  
  PLOT_07 <- ggplot_scan1(
    x = ADD_SCANS,
    map = MAP,
    lodcolumn = TRAITS[1],
    gap = 0,
    col = 'cadetblue',
    lwd = 0.5
  ) +
    theme_minimal(
      base_size = 9
    ) +
    geom_hline(
      yintercept = PT[as.character(ALPHA[1]), TRAITS[1]],
      color = 'firebrick',
      linetype = c(2)
    ) +
    geom_hline(
      yintercept = PT[as.character(ALPHA[2]), TRAITS[1]],
      color = 'firebrick',
      linetype = c(3)
    ) +
    scale_y_continuous(
      limits = c(0, MAXLOD)
    ) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = 'none'
    ) +
    labs(
      y = NULL
    )
  PLOT_13 <- ggplot_scan1(
    x = ADD_SCANS,
    map = MAP,
    lodcolumn = TRAITS[2],
    gap = 0,
    col = 'cadetblue',
    lwd = 0.5
  ) +
    theme_minimal(
      base_size = 9
    ) +
    geom_hline(
      yintercept = PT[as.character(ALPHA[1]), TRAITS[2]],
      color = 'firebrick',
      linetype = c(2)
    ) +
    geom_hline(
      yintercept = PT[as.character(ALPHA[2]), TRAITS[2]],
      color = 'firebrick',
      linetype = c(3)
    ) +
    scale_y_continuous(
      limits = c(0, MAXLOD)
    ) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = 'none'
    ) +
    labs(
      y = NULL
    )
  PLOT_19 <- ggplot_scan1(
    x = ADD_SCANS,
    map = MAP,
    lodcolumn = TRAITS[3],
    gap = 0,
    col = 'cadetblue',
    lwd = 0.5
  ) +
    theme_minimal(
      base_size = 9
    ) +
    geom_hline(
      yintercept = PT[as.character(ALPHA[1]), TRAITS[3]],
      color = 'firebrick',
      linetype = c(2)
    ) +
    geom_hline(
      yintercept = PT[as.character(ALPHA[2]), TRAITS[3]],
      color = 'firebrick',
      linetype = c(3)
    ) +
    scale_y_continuous(
      limits = c(0, MAXLOD)
    ) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = 'none'
    ) +
    labs(
      y = NULL
    )
  
  rm(TRAITS, PT, SCAN1, MAXLOD)
  
  
  # Add timepoint labels
  PLOT_07 <- plot_grid(
    PLOT_07 + 
      theme(
        plot.margin = unit(c(0,1,0,0), "lines")
      ),
    ggdraw() + 
      draw_label(
        "7 months",
        size = 9,
        x = 0.5,
        hjust = 0.5,
        angle = 270
      ),
    ncol = 2,
    rel_widths = c(1, 0.01)
  )
  PLOT_13 <- plot_grid(
    PLOT_13 + 
      theme(
        plot.margin = unit(c(0,1,0,0), "lines"), 
        strip.text = element_blank()
      ),
    ggdraw() + 
      draw_label(
        "13 months",
        size = 9,
        x = 0.5,
        hjust = 0.5,
        angle = 270
      ),
    ncol = 2,
    rel_widths = c(1, 0.01)
  )
  PLOT_19 <- plot_grid(
    PLOT_19 + 
      theme(
        plot.margin = unit(c(0,1,0,0), "lines"), 
        strip.text = element_blank()
      ),
    ggdraw() + 
      draw_label(
        "19 months",
        size = 9,
        x = 0.5,
        hjust = 0.5,
        angle = 270
      ),
    ncol = 2,
    rel_widths = c(1, 0.01)
  )
  
  # Combine add scan plots
  SCAN_PLOTS <- plot_grid(
    PLOT_07 + 
      theme(
        plot.margin = unit(c(0,0.1,0,0), "lines")
      ), 
    PLOT_13 + 
      theme(
        plot.margin = unit(c(2,0.1,0,0), "lines")
      ), 
    PLOT_19 + 
      theme(
        plot.margin = unit(c(2,0.1,0,0), "lines")
      ),
    ncol = 1
  )
  
  rm(PLOT_07, PLOT_13, PLOT_19)
  
  YLAB <- ggdraw() + 
    draw_label(
      "LOD score",
      size = 9,
      x = 0.5, 
      y = 0.5,
      vjust = 0.5,
      hjust = 0.5,
      angle = 90
    )
  
  SCAN_PLOTS <- plot_grid(
    YLAB + 
      theme(
        plot.margin = unit(c(2,0.25,0,0), "lines")
      ), 
    SCAN_PLOTS,
    ncol = 2,
    rel_widths = c(0.025,1)
  )
  
  rm(YLAB)
  
  TITLE <- ggdraw() + 
    draw_label(
      Info_Table$Name[Info_Table$Phenotype == PHENO],
      size = 10,
      x = 0,
      hjust = 0,
      angle = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  
  CAPTION_TEXT <- paste0(
    'The plots above show the results of the genome-wide linkage scans for ',
    ifelse(
      PHENO %in% c('hct', 'hgb', 'n.lymph', 'n.neut', 'n.mono', 'n.eos', 'plt'),
      tolower(Info_Table$TraditionalName[Info_Table$Phenotype == PHENO]),
      Info_Table$Abbreviation[Info_Table$Phenotype == PHENO]
    ),
    ' conditioned on sex and generation at each timepoint.', 
    ' The y-axis shows the LOD score ', 
    # '(base-10 logarithm of the ratio of the likelihood of a QTL at the given position to the likelihood of no QTL) ', 
    'and the x-axis indicates the position on each chromosome.',
    ' The dashed red lines are the trait and timepoint specific permutation thresholds for significance at an alpha level of 0.05 and the dotted red lines are the permutation thresholds for significance at an alpha level of 0.2.'
  )
  i <- 0
  for(INDEX in seq(150, nchar(CAPTION_TEXT), by = 150)){
    POSITION <- INDEX + i
    if(substr(CAPTION_TEXT, POSITION, POSITION) != ' '){
      for(j in 1:20){
        if(substr(CAPTION_TEXT, POSITION + j, POSITION + j) == ' '){
          POSITION <- POSITION + j
          break
        }
        if(substr(CAPTION_TEXT, POSITION - j, POSITION - j) == ' '){
          POSITION <- POSITION - j
          break
        }
      }
    }
    CAPTION_TEXT <- paste0(
      substr(CAPTION_TEXT, 1, POSITION),
      '\n',
      substr(CAPTION_TEXT, POSITION + 1, 999999)
    )
    i <- i + 2
  }
  rm(i, INDEX, POSITION, j)
  
  CAPTION <- ggdraw() + 
    draw_label(
      CAPTION_TEXT,
      size = 7,
      x = 0,
      y = 0.8,
      hjust = 0,
      vjust = 1,
      angle = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 3)
    )
  
  rm(CAPTION_TEXT)
  
  FULL_PLOT <- plot_grid(
    TITLE,
    SCAN_PLOTS,
    CAPTION,
    nrow = 3,
    rel_heights = c(0.1, 1, 0.15)
  )
  
  
  rm(TITLE, SCAN_PLOTS, CAPTION)
  
  pdf(
    paste0(
      'figures/genome_wide_linkage_scans/',
      PHENO,
      '.pdf'
    ),
    width = 7, height = 7
  )
  plot(FULL_PLOT)
  dev.off()
  
  rm(FULL_PLOT)
}

rm(PHENO)

#####


################################################################################
## Lifespan linkage scan plot #################

PHENO <- 'Lifespan'


# Create initial addscan plots
PT <- PERM_THRESH[, PHENO, drop = FALSE]
SCAN1 <- ADD_SCANS[, PHENO, drop = FALSE]
MAXLOD <- maxlod(SCAN1)
MAXLOD <- ceiling(MAXLOD/2.5)*2.5

SCAN_PLOT <- ggplot_scan1(
  x = ADD_SCANS,
  map = MAP,
  lodcolumn = PHENO,
  gap = 0,
  col = 'cadetblue',
  lwd = 0.5
) +
  theme_minimal(
    base_size = 9
  ) +
  geom_hline(
    yintercept = PT[as.character(ALPHA[1]), PHENO],
    color = 'firebrick',
    linetype = c(2)
  ) +
  geom_hline(
    yintercept = PT[as.character(ALPHA[2]), PHENO],
    color = 'firebrick',
    linetype = c(3)
  ) +
  scale_y_continuous(
    limits = c(0, MAXLOD)
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )

rm(PT, SCAN1, MAXLOD)


TITLE <- ggdraw() + 
  draw_label(
    PHENO,
    size = 10,
    x = 0,
    hjust = 0,
    angle = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

CAPTION_TEXT <- paste0(
  'The plot above shows the results of the genome-wide linkage scan for ',
  'lifespan',
  ' conditioned on sex and generation.', 
  ' The y-axis shows the LOD score ', 
  # '(base-10 logarithm of the ratio of the likelihood of a QTL at the given position to the likelihood of no QTL) ', 
  'and the x-axis indicates the position on each chromosome.',
  ' The dashed red line is the permutation threshold for significance at an alpha level of 0.05 and the dotted red line is the permutation threshold for significance at an alpha level of 0.2.'
)
i <- 0
for(INDEX in seq(155, nchar(CAPTION_TEXT), by = 155)){
  POSITION <- INDEX + i
  if(substr(CAPTION_TEXT, POSITION, POSITION) != ' '){
    for(j in 1:20){
      if(substr(CAPTION_TEXT, POSITION + j, POSITION + j) == ' '){
        POSITION <- POSITION + j
        break
      }
      if(substr(CAPTION_TEXT, POSITION - j, POSITION - j) == ' '){
        POSITION <- POSITION - j
        break
      }
    }
  }
  CAPTION_TEXT <- paste0(
    substr(CAPTION_TEXT, 1, POSITION),
    '\n',
    substr(CAPTION_TEXT, POSITION + 1, 999999)
  )
  i <- i + 2
}
rm(i, INDEX, POSITION, j)

CAPTION <- ggdraw() + 
  draw_label(
    CAPTION_TEXT,
    size = 7,
    x = 0,
    y = 0.8,
    hjust = 0,
    vjust = 1,
    angle = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 3)
  )

rm(CAPTION_TEXT)

FULL_PLOT <- plot_grid(
  TITLE,
  SCAN_PLOT,
  CAPTION,
  nrow = 3,
  rel_heights = c(0.1, 1.2/3, 0.15)
)


rm(TITLE, SCAN_PLOT, CAPTION)

pdf(
  paste0(
    'figures/genome_wide_linkage_scans/',
    PHENO,
    '.pdf'
  ),
  width = 7, height = 3.5
)
plot(FULL_PLOT)
dev.off()

rm(FULL_PLOT)

rm(PHENO)

#####


################################################################################
## Lifespan conditioned on blood traits linkage scan plot #################

PHENO <- 'Lifespan'


TRAITS <- c('LifespanCondRDW', 'LifespanCondHDW', 'LifespanCondNLR')
PT <- PERM_THRESH[, 'Lifespan']
SCAN1 <- ADD_SCANS[, c('Lifespan', TRAITS)]
MAXLOD <- maxlod(SCAN1)
MAXLOD <- ceiling(MAXLOD/2.5)*2.5


# Create initial addscan plots
PLOT_RDW <- ggplot_scan1(
  x = ADD_SCANS,
  map = MAP,
  lodcolumn = c('Lifespan', TRAITS[1]),
  gap = 0,
  col = c('gray80', 'cadetblue'),
  lwd = 0.5
) +
  theme_minimal(
    base_size = 9
  ) +
  geom_hline(
    yintercept = PT[as.character(ALPHA[1])],
    color = 'firebrick',
    linetype = c(2)
  ) +
  geom_hline(
    yintercept = PT[as.character(ALPHA[2])],
    color = 'firebrick',
    linetype = c(3)
  ) +
  scale_y_continuous(
    limits = c(0, MAXLOD)
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  ) +
  labs(
    y = NULL
  )
PLOT_HDW <- ggplot_scan1(
  x = ADD_SCANS,
  map = MAP,
  lodcolumn = c('Lifespan', TRAITS[2]),
  gap = 0,
  col = c('gray80', 'cadetblue'),
  lwd = 0.5
) +
  theme_minimal(
    base_size = 9
  ) +
  geom_hline(
    yintercept = PT[as.character(ALPHA[1])],
    color = 'firebrick',
    linetype = c(2)
  ) +
  geom_hline(
    yintercept = PT[as.character(ALPHA[2])],
    color = 'firebrick',
    linetype = c(3)
  ) +
  scale_y_continuous(
    limits = c(0, MAXLOD)
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  ) +
  labs(
    y = NULL
  )
PLOT_NLR <- ggplot_scan1(
  x = ADD_SCANS,
  map = MAP,
  lodcolumn = c('Lifespan', TRAITS[3]),
  gap = 0,
  col = c('gray80', 'cadetblue'),
  lwd = 0.5
) +
  theme_minimal(
    base_size = 9
  ) +
  geom_hline(
    yintercept = PT[as.character(ALPHA[1])],
    color = 'firebrick',
    linetype = c(2)
  ) +
  geom_hline(
    yintercept = PT[as.character(ALPHA[2])],
    color = 'firebrick',
    linetype = c(3)
  ) +
  scale_y_continuous(
    limits = c(0, MAXLOD)
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  ) +
  labs(
    y = NULL
  )

rm(TRAITS, PT, SCAN1, MAXLOD)


# Add timepoint labels
PLOT_RDW <- plot_grid(
  PLOT_RDW + 
    theme(
      plot.margin = unit(c(0,1,0,0), "lines")
    ),
  ggdraw() + 
    draw_label(
      "RDW",
      size = 9,
      x = 0.5,
      hjust = 0.5,
      angle = 270
    ) + 
    theme(
      plot.margin = unit(c(0,0.1,0,0), "lines")
    ),
  ncol = 2,
  rel_widths = c(1, 0.01)
)
PLOT_HDW <- plot_grid(
  PLOT_HDW + 
    theme(
      plot.margin = unit(c(0,1,0,0), "lines"), 
      strip.text = element_blank()
    ),
  ggdraw() + 
    draw_label(
      "HDW",
      size = 9,
      x = 0.5,
      hjust = 0.5,
      angle = 270
    ) + 
    theme(
      plot.margin = unit(c(0,0.1,0,0), "lines")
    ),
  ncol = 2,
  rel_widths = c(1, 0.01)
)
PLOT_NLR <- plot_grid(
  PLOT_NLR + 
    theme(
      plot.margin = unit(c(0,1,0,0), "lines"), 
      strip.text = element_blank()
    ),
  ggdraw() + 
    draw_label(
      "NLR",
      size = 9,
      x = 0.5,
      hjust = 0.5,
      angle = 270
    ) + 
    theme(
      plot.margin = unit(c(0,0.1,0,0), "lines")
    ),
  ncol = 2,
  rel_widths = c(1, 0.01)
)

# Combine add scan plots
SCAN_PLOTS <- plot_grid(
  PLOT_RDW + 
    theme(
      plot.margin = unit(c(0,0.1,0,0), "lines")
    ), 
  PLOT_HDW + 
    theme(
      plot.margin = unit(c(2,0.1,0,0), "lines")
    ), 
  PLOT_NLR + 
    theme(
      plot.margin = unit(c(2,0.1,0,0), "lines")
    ),
  ncol = 1
)

rm(PLOT_RDW, PLOT_HDW, PLOT_NLR)

YLAB <- ggdraw() + 
  draw_label(
    "LOD score",
    size = 9,
    x = 0.5, 
    y = 0.5,
    vjust = 0.5,
    hjust = 0.5,
    angle = 90
  )

SCAN_PLOTS <- plot_grid(
  YLAB + 
    theme(
      plot.margin = unit(c(2,0.25,0,0), "lines")
    ), 
  SCAN_PLOTS,
  ncol = 2,
  rel_widths = c(0.025,1)
)

rm(YLAB)


TITLE <- ggdraw() + 
  draw_label(
    'Lifespan conditioned on 7-month RDW, HDW, or NLR',
    size = 10,
    x = 0,
    hjust = 0,
    angle = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

CAPTION_TEXT <- paste0(
  'The plot above shows the results of the genome-wide linkage scans for ',
  'lifespan',
  ' conditioned on sex, generation and one of the three blood traits shown to change with age and be predictive of lifespan (RDW, HDW, or NLR).', 
  ' The y-axis shows the LOD score ', 
  # '(base-10 logarithm of the ratio of the likelihood of a QTL at the given position to the likelihood of no QTL) ',
  '(the strength of the genetic association) ',
  'and the x-axis indicates the position on each chromosome.',
  ' The blue-green lines indicate the LOD score for the genetic association of lifespan conditioned on the corresponding blood trait (indicated to the right of the plot) and the light gray lines indicate the genetic association of lifespan ignoring the value of the blood trait.',
  ' The dashed red line is the permutation threshold for significance at an alpha level of 0.05 and the dotted red line is the permutation threshold for significance at an alpha level of 0.2.'
)
i <- 0
for(INDEX in seq(150, nchar(CAPTION_TEXT), by = 150)){
  POSITION <- INDEX + i
  if(substr(CAPTION_TEXT, POSITION, POSITION) != ' '){
    for(j in 1:20){
      if(substr(CAPTION_TEXT, POSITION + j, POSITION + j) == ' '){
        POSITION <- POSITION + j
        break
      }
      if(substr(CAPTION_TEXT, POSITION - j, POSITION - j) == ' '){
        POSITION <- POSITION - j
        break
      }
    }
  }
  CAPTION_TEXT <- paste0(
    substr(CAPTION_TEXT, 1, POSITION),
    '\n',
    substr(CAPTION_TEXT, POSITION + 1, 999999)
  )
  i <- i + 2
}
rm(i, INDEX, POSITION, j)

CAPTION <- ggdraw() + 
  draw_label(
    CAPTION_TEXT,
    size = 7,
    x = 0,
    y = 0.8,
    hjust = 0,
    vjust = 1,
    angle = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 3)
  )

rm(CAPTION_TEXT)

FULL_PLOT <- plot_grid(
  TITLE,
  SCAN_PLOTS,
  CAPTION,
  nrow = 3,
  rel_heights = c(0.1, 1, 0.15)
)


rm(TITLE, SCAN_PLOTS, CAPTION)

pdf(
  'figures/genome_wide_linkage_scans/Lifespan_conditioned_on_blood_traits.pdf',
  width = 7, height = 7
)
plot(FULL_PLOT)
dev.off()

rm(FULL_PLOT)

rm(PHENO)

#####

rm(ALPHA, PERM_THRESH)


################################################################################
## Heatmap ################################

PHENO_LOD_TABLE <- create_pheno_lod_table(
  add_scan = subset_scan1(ADD_SCANS, lodcolumn = colnames(PERM_SCANS)), 
  map_df = MAP_DF,
  id_string = NULL
)

LEVELS <- c(
  'Lifespan',
  paste(
    rep(CBCphenos, each = 3), 
    rep(c('07', '13', '19'), length(CBCphenos)), 
    sep = '_'
  )
)

YLIM <- c(0.5, length(LEVELS) + 0.5)
YBREAKS <- c(1, seq(3, 51, 3))
YLABELS <- c(
  'Lifespan',
  Info_Table$Abbreviation
)

LOD_PLOT <- PHENO_LOD_TABLE %>% 
  mutate(
    Phenotype = factor(Phenotype, levels = LEVELS)
  ) %>% 
  lod_heat_plot(
    ylim = YLIM, y_breaks = YBREAKS, y_labels = YLABELS
  )

rm(PHENO_LOD_TABLE, LEVELS, YLIM, YBREAKS, YLABELS)

# Add horizontal lines between phenotypes

# Fine lines between same trait at different timepoint
for(i in seq(1.5, 52.5, 1)){
  LOD_PLOT <- LOD_PLOT +
    geom_hline(yintercept = i, color = 'gray50', size = .025)
}

# Large lines between totally different phenotypes
LOD_PLOT <- LOD_PLOT +
  geom_hline(yintercept = 1.5, color = 'gray50', size = .5)

MAJOR_BREAKS <- c(seq(3, 51, 3))
for(i in MAJOR_BREAKS[1:(length(MAJOR_BREAKS) - 1)]){
  LOD_PLOT <- LOD_PLOT +
    geom_hline(yintercept = i + 1.5, color = 'gray50', size = .5)
}

rm(i, MAJOR_BREAKS)

pdf(
  'figures/genome_wide_linkage_scans/Heatmap_of_LOD_scores.pdf',
  width = 7, height = 7
)
LOD_PLOT +
  labs(
    title = 'Heatmap of LOD scores for lifespan and blood traits',
    caption = paste0(
      'LOD scores above ten are truncated and LOD scores below five are dropped (indicated by the pale tan background color).',
      '\nFor each blood trait, the three timepoints are shown by three individual rows separated by fine gray lines.', 
      ' The earliest timepoint\n(7 months) is represented by the bottommost row for the corresponding blood trait,', 
      ' and the latest timepoint (19 months) is\nindicated by the topmost row.'
    )
  ) +
  theme(
    plot.caption = element_text(size = 7, hjust = 0)
  )
dev.off()

rm(LOD_PLOT)

#####


################################################################################
##  ################################

#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####