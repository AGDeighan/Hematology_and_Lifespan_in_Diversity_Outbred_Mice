# 2021-06-29
################################################################################
#
#   This script creates genome-wide variant-association scan plots for lifespan,
#   RDW, and HDW
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
library(qtl2ggplot)
library(cowplot)

#####


################################################################################
## load data  ###########################################

ADD_SCANS <- readRDS('results/variant_associations/LifespanHema_GenomeWideVariantScans.rds')
PERM_SCANS <- readRDS('results/variant_associations/LifespanHema_10000Perms.rds')

#####


################################################################################
## setup  ###########################################

ALPHA <- c(0.05, 0.2)
PERM_THRESH <- summary_scan1perm(PERM_SCANS, alpha = ALPHA)
#      Lifespan rdw_07 hdw_07
# 0.05     5.50   5.43   5.51
# 0.2      4.82   4.80   4.82

MAP <- vector('list', length = length(unique(ADD_SCANS$snpinfo$chr)))
names(MAP) <- unique(ADD_SCANS$snpinfo$chr)
for(i in unique(ADD_SCANS$snpinfo$chr)){
  POS <- ADD_SCANS$snpinfo %>% 
    filter(chr == i) %>% 
    select(pos) %>% 
    unlist() %>% 
    as.numeric()
  
  names(POS) <-  ADD_SCANS$snpinfo %>% 
    filter(chr == i) %>% 
    select(snp_id) %>% 
    unlist() %>% 
    as.character()
  
  MAP[[i]] <- POS
}
rm(i, POS)

#####


################################################################################
## Lifespan, covar Sex and Gen  ###########################################

PHENO <- 'Lifespan'

PT <- PERM_THRESH[, PHENO]
SCAN1 <- ADD_SCANS$lod[, PHENO, drop = FALSE]
MAXLOD <- maxlod(SCAN1)
MAXLOD <- ceiling(max(c(MAXLOD, PT))/2.5)*2.5

SCAN_PLOT <- ggplot_scan1(
  x = SCAN1,
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
  'The plot above shows the results of the genome-wide variant-association scan for ',
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

tiff(
  'figures/variant_associations/GenomeWide_VarScan_Lifespan_CondSexGen.tiff',
  width = 7, height = 3.5, units = 'in',
  res = 200,
  type = 'quartz'
)
plot(FULL_PLOT)
dev.off()

rm(FULL_PLOT)

rm(PHENO)

#####


################################################################################
## Lifespan conditioned on blood traits linkage scan plot #################

PHENO <- 'Lifespan'

TRAITS <- c('LifespanCondRDW', 'LifespanCondHDW')
PT <- PERM_THRESH[, 'Lifespan']
SCAN1 <- ADD_SCANS$lod[, c('Lifespan', TRAITS), drop = FALSE]
MAXLOD <- maxlod(SCAN1)
MAXLOD <- ceiling(max(c(MAXLOD, PT))/2.5)*2.5


# Create initial addscan plots
PLOT_RDW <- ggplot_scan1(
  x = SCAN1,
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
  x = SCAN1,
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
# PLOT_NLR <- ggplot_scan1(
#   x = ADD_SCANS,
#   map = MAP,
#   lodcolumn = c('Lifespan', TRAITS[3]),
#   gap = 0,
#   col = c('gray80', 'cadetblue'),
#   lwd = 0.5
# ) +
#   theme_minimal(
#     base_size = 9
#   ) +
#   geom_hline(
#     yintercept = PT[as.character(ALPHA[1])],
#     color = 'firebrick',
#     linetype = c(2)
#   ) +
#   geom_hline(
#     yintercept = PT[as.character(ALPHA[2])],
#     color = 'firebrick',
#     linetype = c(3)
#   ) +
#   scale_y_continuous(
#     limits = c(0, MAXLOD)
#   ) +
#   theme(
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     legend.position = 'none'
#   ) +
#   labs(
#     y = NULL
#   )

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
# PLOT_NLR <- plot_grid(
#   PLOT_NLR + 
#     theme(
#       plot.margin = unit(c(0,1,0,0), "lines"), 
#       strip.text = element_blank()
#     ),
#   ggdraw() + 
#     draw_label(
#       "NLR",
#       size = 9,
#       x = 0.5,
#       hjust = 0.5,
#       angle = 270
#     ) + 
#     theme(
#       plot.margin = unit(c(0,0.1,0,0), "lines")
#     ),
#   ncol = 2,
#   rel_widths = c(1, 0.01)
# )

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
  # PLOT_NLR + 
  #   theme(
  #     plot.margin = unit(c(2,0.1,0,0), "lines")
  #   ),
  ncol = 1
)

rm(PLOT_RDW, PLOT_HDW)

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
    'Lifespan conditioned on 7-month RDW or HDW',
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
  'The plot above shows the results of the genome-wide variant-association scans for ',
  'lifespan',
  ' conditioned on sex, generation and either RDW or HDW.', 
  ' The y-axis shows the LOD score ', 
  # '(base-10 logarithm of the ratio of the likelihood of a QTL at the given position to the likelihood of no QTL) ',
  '(the strength of the genetic association) ',
  'and the x-axis indicates the position on each chromosome.',
  ' The blue-green lines indicate the LOD score for the genetic association of lifespan conditioned on the corresponding blood trait (indicated to the right of the plot) and the light gray lines indicate the genetic association of lifespan ignoring the blood trait.',
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


tiff(
  'figures/variant_associations/GenomeWide_VarScan_Lifespan_CondSexGenBlood.tiff',
  width = 7, height = 7, units = 'in',
  res = 200,
  type = 'quartz'
)
plot(FULL_PLOT)
dev.off()

rm(FULL_PLOT)

rm(PHENO)

#####


################################################################################
## RDW, covar Sex and Gen  ###########################################

PHENO <- 'rdw_07'

PT <- PERM_THRESH[, PHENO]
SCAN1 <- ADD_SCANS$lod[, PHENO, drop = FALSE]
MAXLOD <- maxlod(SCAN1)
MAXLOD <- ceiling(max(c(MAXLOD, PT))/2.5)*2.5

SCAN_PLOT <- ggplot_scan1(
  x = SCAN1,
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
  )

rm(PT, SCAN1, MAXLOD)


TITLE <- ggdraw() + 
  draw_label(
    'RDW at 7 months',
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
  'The plot above shows the results of the genome-wide variant-association scan for ',
  'RDW measured at 7 months',
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

tiff(
  'figures/variant_associations/GenomeWide_VarScan_RDW_CondSexGen.tiff',
  width = 7, height = 3.5, units = 'in',
  res = 200,
  type = 'quartz'
)
plot(FULL_PLOT)
dev.off()

rm(FULL_PLOT)

rm(PHENO)

#####


################################################################################
## HDW, covar Sex and Gen  ###########################################

PHENO <- 'hdw_07'

PT <- PERM_THRESH[, PHENO]
SCAN1 <- ADD_SCANS$lod[, PHENO, drop = FALSE]
MAXLOD <- maxlod(SCAN1)
MAXLOD <- ceiling(max(c(MAXLOD, PT))/2.5)*2.5

SCAN_PLOT <- ggplot_scan1(
  x = SCAN1,
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
  )

rm(PT, SCAN1, MAXLOD)


TITLE <- ggdraw() + 
  draw_label(
    'HDW at 7 months',
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
  'The plot above shows the results of the genome-wide variant-association scan for ',
  'HDW measured at 7 months',
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

tiff(
  'figures/variant_associations/GenomeWide_VarScan_HDW_CondSexGen.tiff',
  width = 7, height = 3.5, units = 'in',
  res = 200,
  type = 'quartz'
)
plot(FULL_PLOT)
dev.off()

rm(FULL_PLOT)

rm(PHENO)

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
