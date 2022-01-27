# 2020-04-28
################################################################################
#
#   This script creates the CPH coefficient plots for the marginal effects of
#   the phenotypes on survival
# 
#
#   Author:  Andrew Deighan
#   E-mail: andrew.deighan@jax.org
#
################################################################################

options(na.action = 'na.exclude')
options(stringsAsFactors = FALSE)

#####


################################################################################
## libraries etc #########################################################

library(tidyverse)

#####


################################################################################
## helper functions #########################################################


#####


################################################################################
## load data etc #########################################################

load('data/processed/phenotypic/CBC_AnalysisEnvironment.RData')

LS_EFFECTS <- read_csv(
  'tables/lifespan_aging_hematology/hematology_effects_on_lifespan.csv'
) %>% 
  rename(
    Feature = Trait
  ) %>% 
  left_join(
    Info_Table %>% 
      rename(Feature = Phenotype),
    by = 'Feature'
  ) %>% 
  mutate(
    Feature = factor(
      Feature, 
      levels = Info_Table$Phenotype[Info_Table$Phenotype %in% Feature]
    ),
    Timepoint = str_replace(
      Timepoint, 'Months', 'months'
    ),
    Timepoint = factor(
      Timepoint, 
      levels = paste0(c(7, 13, 19), ' months')
    )
  )

#####


################################################################################
## Prepare some arguments for plotting functions #####################

YLAB <- 'Regression coefficient for lifespan\nmonths per standard deviation change in trait'

HLINE_COLOR <- 'grey50'
HLINE_WIDTH <- 0.4

ERROR_BAR_WIDTH <- 0.4
ERROR_BAR_SIZE <- 0.5

SIGNIFICANCE_CUTOFF <- 0.05
SIGNIFICANCE_LEVELS <- c('< 0.05', '> 0.05')
SIGNIFICANCE_COLORS <- c('firebrick', 'grey20')
POINT_SIZE <- 1.5

BASE_SIZE <- 13

X_AXIS_TEXT_ELEMENT <- element_text(angle = 45, hjust = 1)

Y_BREAKS <- seq(-10, 10, 1)
Y_MINOR_BREAKS <- seq(-10, 10, .5)
Y_LIMITS <- c(-3, 3)


LAYOUT_MATRIX <- matrix(c(rep(1, each = 10),
                          rep(2:3, each = 8)), 
                        ncol = 1)
AXIS_TITLE_SIZE <- 14

#####


################################################################################
## Prepare data for plotting functions #####################

PLOT_DATA <- LS_EFFECTS %>%  
  mutate(
    Coefficient = Effect * TraitSD,
    SE = SE * TraitSD,
    Significance = ifelse(
      FWER_PBS < SIGNIFICANCE_CUTOFF, 
      SIGNIFICANCE_LEVELS[1], 
      SIGNIFICANCE_LEVELS[2]
    ),
    Significance = factor(
      Significance,
      levels = SIGNIFICANCE_LEVELS,
      ordered = TRUE
    ),
    Group = factor(
      Group,
      levels = c(
        'RBC', 'Platelet', 'WBC'
      )
    ),
    Abbreviation = str_replace(
      Abbreviation, 'coumt', 'count'
    ),
    Abbreviation = gsub('([[:upper:]]{1}[[:lower:]]+)( #)', '\\1. Count', Abbreviation),
    Abbreviation = gsub('Platelet\\.', 'Platelet', Abbreviation),
    Abbreviation = gsub('Hct', 'Hematocrit', Abbreviation),
    Abbreviation = gsub('Hgb', 'Hemoglobin', Abbreviation),
    Abbreviation = gsub('#', 'Count', Abbreviation),
    Abbreviation = factor(
      Abbreviation,
      levels = c(
        "RBC count", "Hematocrit", "Hemoglobin",  "CH", "CHCM", "HDW", "MCV", "RDW",
        "Platelet count", "MPV", "MPM",
        "WBC count", "NLR", "Lymph. count", "Neut. count", "Mono. count", "Eos. count"
      )
    )
  )

#####


################################################################################
## Create plot #####################

PLOT <- PLOT_DATA %>% 
  ggplot(aes(x = Abbreviation, y = Coefficient)) + 
  # Add horizontal line at no effect
  geom_hline(yintercept = 0, color = HLINE_COLOR, lwd = 1) +
  # Create +/- 1 SE error bars
  geom_errorbar(
    aes(ymin = Coefficient - SE, ymax = Coefficient + SE),
    width = ERROR_BAR_WIDTH, size = ERROR_BAR_SIZE
  ) +
  # Add points for coefficient (log-2 hazard ratio) estimates
  geom_point(aes(color = Significance), size = POINT_SIZE) +
  # facet by timepoint and CBC group (RBC, WBC, and platelets)
  facet_grid(Timepoint ~ Group, scales = 'free_x', space = 'free_x') + 
  # Add color scale for points, pvalue < 0.05 = red, otherwise black
  scale_color_manual(values = SIGNIFICANCE_COLORS, drop = FALSE) +
  # Set base font size
  theme_minimal(base_size = BASE_SIZE) +
  # Set the angle of the X-labels and the size and justification of the caption
  theme(
    axis.text.x = X_AXIS_TEXT_ELEMENT, 
    panel.spacing.y = unit(1.5, "lines")
    # plot.caption = CAPTION_TEXT_ELEMENT
  ) +
  # Set the y-scale minor and major breaks and axis limits
  scale_y_continuous(
    breaks = Y_BREAKS,
    minor_breaks = Y_MINOR_BREAKS,
    limits = Y_LIMITS,
    exp = c(0, 0)
  ) +
  # Remove color bar for red/black p-value
  guides(color = FALSE) +
  # Add plot labels
  labs(
    title = NULL,
    y = YLAB,
    x = NULL
    # caption = CAPTION
  ); PLOT

  
pdf(
  'figures/lifespan_aging_hematology/HemaLS_Coefficients.pdf',
  width = 7, height = 7.75
  )
plot(PLOT)
dev.off()


#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####