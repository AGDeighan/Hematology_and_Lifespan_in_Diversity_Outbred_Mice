# 2021-05-27
################################################################################
#
#   This script estimates the effect of each of the CBC traits, across 
#   timepoints, on lifespan. 
# 
#   This script is designed to by run on a high performance computing cluster 
#   with at least 32 cores due to the multitude of simulations performed during 
#   the parametric bootstraps to estimate statistical significance of effects.
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

library(tidyverse)
library(survival)
library(survminer)
library(cowplot)
library(ggbeeswarm)

#####


################################################################################
## helper functions #########################################################

regress_out_sex <- function(x, sex, timepoint, id, catch_model = FALSE){
  
  old_x <- x
  
  options(na.action = 'na.exclude')
  if(catch_model){
    f <- function(...){'warning'}
    MODEL <- tryCatch(lme4::lmer(x ~ sex + (1|timepoint) + (1|id)),
                      message=f, warning = f)
    if(is.character(MODEL)){
      if(MODEL == 'warning'){
        MODEL <- lme4::lmer(x ~ sex + (1|timepoint) + (1|id))
        saveRDS(MODEL, paste0('SexDOT_Adj_Model_', 
                              floor(mean(old_x, na.rm = TRUE)), 
                              '_',
                              floor(runif(1,1,9999)),
                              '.rds'))
      }
    }
  } else{
    MODEL <- lme4::lmer(x ~ sex + (1|timepoint) + (1|id))
  }
  
  FIXED_EFFECTS <- lme4::fixef(MODEL)
  SEX_EFFECTS <- FIXED_EFFECTS["sex"]
  x <- x - sex*SEX_EFFECTS
  return(x)
}

make_quantiles <- function(values, levels = 5, 
                           min_cp = -1000000000,
                           max_cp = 1000000000,
                           label_precision = 100){
  PROBS <- seq(0 + 1/levels, 1 - 1/levels, 1/levels)
  CUTPOINTS <- c(min_cp, 
                 as.numeric(quantile(values, probs = PROBS, na.rm = TRUE)), 
                 max_cp)
  QUANTILES <- cut(values, breaks = CUTPOINTS)
  
  NEW_NAMES <- paste0(floor(min(values, na.rm = TRUE) * label_precision) / label_precision,
                      ' to ',
                      floor(CUTPOINTS[2] * label_precision)/ label_precision)
  for(i in 2:(length(CUTPOINTS) - 2)){
    NEW_NAMES <- c(NEW_NAMES,
                   paste0(floor(CUTPOINTS[i] * label_precision)/ label_precision,
                          ' to ',
                          floor(CUTPOINTS[i + 1] * label_precision)/ label_precision))
  }
  NEW_NAMES <- c(NEW_NAMES,
                 paste0(floor(CUTPOINTS[length(CUTPOINTS) - 1] * label_precision)/ label_precision,
                        ' to ',
                        ceiling(max(values, na.rm = TRUE) * label_precision) / label_precision))
  
  QUANTILES <- factor(NEW_NAMES[as.numeric(QUANTILES)],
                      levels = NEW_NAMES)
  
  return(QUANTILES)
  
}

convert_level_names <- function(x, new_names){
  if(!is.factor(x)){
    stop('x must be a factor')
  }
  return(factor(
    new_names[as.numeric(x)],
    levels = new_names
  ))
}

#####


################################################################################
## load data etc #########################################################

COLPAL_SURV <- c('royalblue4', 'cadetblue', 'darkorange3', 'firebrick')
COLPAL_SEX <- c('darkorange3', 'cadetblue')

load('data/processed/phenotypic/CBC_AnalysisEnvironment.RData')

AGE_EFFECTS <- read_csv(
  'tables/lifespan_aging_hematology/age_effects_on_hematology_traits.csv'
) %>% 
  mutate(
    Trait = factor(Trait, levels = CBCphenos)
  )

SEX_EFFECTS <- read_csv(
  'tables/lifespan_aging_hematology/sex_effects_on_hematology_traits.csv'
) %>% 
  mutate(
    Trait = factor(Trait, levels = CBCphenos)
  )

AGESEX_INTERACTION_EFFECTS <- read_csv(
  'tables/lifespan_aging_hematology/agesex_interaction_effects_on_hematology_traits.csv'
) %>% 
  mutate(
    Trait = factor(Trait, levels = CBCphenos)
  )

LIFESPAN_EFFECTS <- read_csv(
  'tables/lifespan_aging_hematology/hematology_effects_on_lifespan.csv'
) %>% 
  mutate(
    Trait = factor(Trait, levels = CBCphenos),
    Timepoint = factor(Timepoint, levels = levels(PhenoData_Raw$Timepoint))
  )

#####


################################################################################
## Determine which traits to log transform before plotting #################

par(mfcol = c(2,4))
for(PHENO in CBCphenos){
  UNIT <- Info_Table$Units[Info_Table$Phenotype == PHENO]
  hist(
    PhenoData_ClumpDateAdj[[PHENO]],
    breaks = 100,
    main = PHENO,
    xlab = UNIT
  )
  hist(
    log2(PhenoData_ClumpDateAdj[[PHENO]]),
    breaks = 100,
    main  = '',
    xlab = paste0('log2[', UNIT, ']')
  )
}
rm(PHENO, UNIT)

# The following variables will be log-transformed before plotting
LOG_PHENOS <- c(
  'hdw', 'rdw', 'wbc', 'nlr', 'n.lymph', 'n.neut', 'n.mono', 'n.eos', 'plt', 'mpv', 'mpm'
)

# Check if any of the variables to be transformed have values less than or
# equal to zero (can't log transform x <= 0)
PhenoData_ClumpDateAdj %>% 
  select(
    MouseID, Lifespan, Timepoint, AgeAtTest,
    all_of(LOG_PHENOS)
  ) %>% 
  filter_at(vars(all_of(LOG_PHENOS)), any_vars(. <= 0))
# A tibble: 0 x 15

#####


################################################################################
## Create summary age, sex, and lifespan plots for each CBC trait ###########

for(PHENO in CBCphenos){
  pdf(
    paste0(
      'figures/lifespan_aging_hematology/aging_and_lifespan_pheno_plots/',
      PHENO,
      '.pdf'
    ),
    width = 7, height = 7
  )
  
  ## Age and sex plot
  
  AGE_PLOT_DATA <- PhenoData_ClumpDateAdj
  AGE_PLOT_DATA[['Pheno']] <- AGE_PLOT_DATA[[PHENO]]
  
  if(PHENO %in% LOG_PHENOS){AGE_PLOT_DATA$Pheno <- log2(AGE_PLOT_DATA$Pheno)}
  
  AGE_PLOT_DATA <- AGE_PLOT_DATA %>% 
    select(
      MouseID, Sex, Timepoint, AgeAtTest, Pheno
    ) %>% 
    mutate(
      Timepoint = factor(
        gsub(' Months', '', Timepoint),
        levels = c('7', '13', '19')
      )
    )
  
  AGE_PLOT <- AGE_PLOT_DATA %>% 
    ggplot() +
    theme_minimal(
      base_size = 9
    ) +
    geom_beeswarm(
      aes(x = Timepoint, y = Pheno, color = Sex),
      groupOnX = TRUE, 
      alpha = 2/3,
      size = 1/3,
      cex = 2/3
    ) +
    geom_boxplot(
      aes(x = Timepoint, y = Pheno),
      outlier.shape = NA,
      alpha = 0
    ) +
    facet_grid(
      Sex ~ .,
      scale = 'fixed'
    ) +
    scale_color_manual(
      values = COLPAL_SEX
    ) +
    theme(
      legend.position = 'top',
      panel.spacing = unit(2, "lines")
    ) +
    labs(
      title = '',
      color = NULL,
      x = 'Age (months)',
      y = ifelse(
        PHENO %in% LOG_PHENOS,
        paste0(
          Info_Table$Abbreviation[Info_Table$Phenotype == PHENO],
          ' (log2[', Info_Table$Units[Info_Table$Phenotype == PHENO], '])'
        ),
        paste0(
          Info_Table$Abbreviation[Info_Table$Phenotype == PHENO],
          ' (', Info_Table$Units[Info_Table$Phenotype == PHENO], ')'
        )
      )
    )
  
  rm(AGE_PLOT_DATA)
  
  
  ## Surv curvs
  SURV_PLOT_DATA <- PhenoData_ClumpDateAdj
  SURV_PLOT_DATA[['Pheno']] <- SURV_PLOT_DATA[[PHENO]]
  
  SURV_PLOT_DATA <- SURV_PLOT_DATA %>% 
    select(
      MouseID, Sex, Timepoint, Lifespan, AgeAtTest, Pheno
    ) %>% 
    mutate(
      Sex = ifelse(Sex == 'Female', 0, 1),
      Pheno = regress_out_sex(
        x = Pheno, sex = Sex, timepoint = Timepoint, id = MouseID
      ),
      ResidualLife = Lifespan - AgeAtTest
    ) %>% 
    select(
      MouseID, Timepoint, ResidualLife, Pheno
    ) %>% 
    group_by(Timepoint) %>% 
    mutate(
      Pheno = convert_level_names(
        make_quantiles(Pheno, levels = 4),
        new_names = c('First', 'Second', 'Third', 'Fourth')
      )
    ) %>% 
    ungroup() %>% 
    mutate(
      Died = TRUE
    ) %>% 
    select(
      MouseID, Timepoint, ResidualLife, Died, Pheno
    )
  
  PLOT_DATA <- SURV_PLOT_DATA %>% filter(Timepoint == '7 Months')
  SURV_PLOT <- ggsurvplot(
    fit = survfit(
      Surv(ResidualLife, Died) ~ Pheno, 
      data = PLOT_DATA
    ),
    size = 0.4,
    censor.shape = 124,
    palette = COLPAL_SURV,
    xlim = c(0,51),
    data = PLOT_DATA,
    alpha = 0.75,
    legend.title = 'Quartile: ',
    legend.labs = c('First', 'Second', 'Third', 'Fourth'),
    ggtheme = theme_survminer(
      font.legend = c(6, "plain", "black"),
      font.tickslab = c(6, "plain", "black")
    )
  )
  SURV_CURV_07 <- SURV_PLOT$plot +
    draw_label(
      paste0(
        'p = ', 
        signif(
          LIFESPAN_EFFECTS$PValue_PBS[LIFESPAN_EFFECTS$Trait == PHENO & LIFESPAN_EFFECTS$Timepoint == '7 Months'], 
          digits = 3
        ),
        '\nFWER = ',
        signif(
          LIFESPAN_EFFECTS$FWER_PBS[LIFESPAN_EFFECTS$Trait == PHENO & LIFESPAN_EFFECTS$Timepoint == '7 Months'], 
          digits = 3
        )
      ),
      size = 7,
      x = 36,
      y = 1,
      hjust = 0,
      vjust = 1,
      angle = 0
    ) +
    scale_x_continuous(
      exp = c(0,0),
      breaks = seq(0,60,12),
      minor_breaks = seq(0,60,6)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, .2),
      minor_breaks = seq(0, 1, .1),
      exp = c(0,0),
      limits = c(0, 1.01)
    ) +
    theme(
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7)
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = NULL
    )
  
  PLOT_DATA <- SURV_PLOT_DATA %>% filter(Timepoint == '13 Months')
  SURV_PLOT <- ggsurvplot(
    fit = survfit(
      Surv(ResidualLife, Died) ~ Pheno, 
      data = PLOT_DATA
    ),
    size = 0.4,
    censor.shape = 124,
    palette = COLPAL_SURV,
    xlim = c(0,51),
    data = PLOT_DATA,
    alpha = 0.75,
    legend.title = 'Quartile: ',
    legend.labs = c('First', 'Second', 'Third', 'Fourth'),
    ggtheme = theme_survminer(
      font.legend = c(6, "plain", "black"),
      font.tickslab = c(6, "plain", "black")
    )
  )
  SURV_CURV_13 <- SURV_PLOT$plot +
    draw_label(
      paste0(
        'p = ', 
        signif(
          LIFESPAN_EFFECTS$PValue_PBS[LIFESPAN_EFFECTS$Trait == PHENO & LIFESPAN_EFFECTS$Timepoint == '13 Months'], 
          digits = 3
        ),
        '\nFWER = ',
        signif(
          LIFESPAN_EFFECTS$FWER_PBS[LIFESPAN_EFFECTS$Trait == PHENO & LIFESPAN_EFFECTS$Timepoint == '13 Months'], 
          digits = 3
        )
      ),
      size = 7,
      x = 36,
      y = 1,
      hjust = 0,
      vjust = 1,
      angle = 0
    ) +
    scale_x_continuous(
      exp = c(0,0),
      breaks = seq(0,60,12),
      minor_breaks = seq(0,60,6)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, .2),
      minor_breaks = seq(0, 1, .1),
      exp = c(0,0),
      limits = c(0, 1.01)
    ) +
    theme(
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7)
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = NULL
    )
  
  PLOT_DATA <- SURV_PLOT_DATA %>% filter(Timepoint == '19 Months')
  SURV_PLOT <- ggsurvplot(
    fit = survfit(
      Surv(ResidualLife, Died) ~ Pheno, 
      data = PLOT_DATA
    ),
    size = 0.4,
    censor.shape = 124,
    palette = COLPAL_SURV,
    xlim = c(0,51),
    data = PLOT_DATA,
    alpha = 0.75,
    legend.title = 'Quartile: ',
    legend.labs = c('First', 'Second', 'Third', 'Fourth'),
    ggtheme = theme_survminer(
      font.legend = c(6, "plain", "black"),
      font.tickslab = c(6, "plain", "black")
    )
  )
  SURV_CURV_19 <- SURV_PLOT$plot +
    draw_label(
      paste0(
        'p = ', 
        signif(
          LIFESPAN_EFFECTS$PValue_PBS[LIFESPAN_EFFECTS$Trait == PHENO & LIFESPAN_EFFECTS$Timepoint == '19 Months'], 
          digits = 3
        ),
        '\nFWER = ',
        signif(
          LIFESPAN_EFFECTS$FWER_PBS[LIFESPAN_EFFECTS$Trait == PHENO & LIFESPAN_EFFECTS$Timepoint == '19 Months'], 
          digits = 3
        )
      ),
      size = 7,
      x = 36,
      y = 1,
      hjust = 0,
      vjust = 1,
      angle = 0
    ) +
    scale_x_continuous(
      exp = c(0,0),
      breaks = seq(0,60,12),
      minor_breaks = seq(0,60,6)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, .2),
      minor_breaks = seq(0, 1, .1),
      exp = c(0,0),
      limits = c(0, 1.01)
    ) +
    theme(
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7)
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = NULL
    )
  
  SURV_CURVS <- plot_grid(
    SURV_CURV_07 + 
      theme(
        legend.position = 'none', axis.ticks.x = element_blank(), axis.text.x = element_blank()
      ),
    SURV_CURV_13 + 
      theme(
        legend.position = 'none', axis.ticks.x = element_blank(), axis.text.x = element_blank()
      ),
    SURV_CURV_19 + 
      theme(
        legend.position = 'none'
      ),
    nrow = 3, 
    ncol = 1,
    align = 'h',
    axis = 'lr'
  )
  
  TP_LABELS <- plot_grid(
    ggdraw() + 
      draw_label(
        "7 months",
        size = 7,
        x = 0.5,
        hjust = 0.5,
        angle = 270
      ),
    ggdraw() + 
      draw_label(
        "13 months",
        size = 7,
        x = 0.5,
        hjust = 0.5,
        angle = 270
      ),
    ggdraw() + 
      draw_label(
        "19 months",
        size = 7,
        x = 0.5,
        hjust = 0.5,
        angle = 270
      ),
    nrow = 3, align = 'v'
  ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 7, 0, 0)
    )
  
  LEGEND <- plot_grid(
    get_legend(
      SURV_CURV_07 + 
        scale_color_manual(
          values = c('royalblue4', 'cadetblue', 'darkorange3', 'firebrick'), 
          name = 'Quartile:  ', 
          labels = c('First    ', 'Second    ', 'Third    ', 'Fourth    ')
        ) + 
        theme(
          legend.position = 'top'
        )
    )
  )
  
  YLAB <- ggdraw() + 
    draw_label(
      "Survival probability",
      size = 9,
      x = 0.5, 
      y = 0.5,
      vjust = 0.5,
      hjust = 0.5,
      angle = 90
    )
  
  XLAB <- ggdraw() + 
    draw_label(
      "Survival time\n(months after measurement)",
      size = 9,
      x = .5,
      hjust = 0.5,
      angle = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 7, 0)
    )
  
  COMBINED_SURV_CURVES <- plot_grid(
    YLAB,
    plot_grid(
      LEGEND,
      plot_grid(
        plot_grid(
          SURV_CURVS, 
          TP_LABELS,
          ncol = 2,
          rel_widths = c(1, 0.01)
        ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 14, 0)
          ),
        XLAB  +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 14, 0)
          ),
        nrow = 2,
        rel_heights = c(1,.01)
      ),
      ncol = 1,
      # rel_heights values control vertical title margins
      rel_heights = c(0.05, 1)
    ),
    ncol = 2, rel_widths = c(0.05,1)
  )
  
  rm(
    SURV_PLOT_DATA,
    PLOT_DATA,
    SURV_PLOT, 
    SURV_CURV_07, SURV_CURV_13, SURV_CURV_19, 
    SURV_CURVS, TP_LABELS, LEGEND, XLAB, YLAB
  )
  
  
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
  
  SD_CAPTION <- paste0(
    'The standard deviation (SD) of ', 
    Info_Table$Abbreviation[Info_Table$Phenotype == PHENO],
    ' is ', 
    signif(AGE_EFFECTS$TraitSD[AGE_EFFECTS$Trait == PHENO], digits = 3),
    ' ',
    Info_Table$Units[Info_Table$Phenotype == PHENO],
    '.'
  )
  SEX_CAPTION <- ifelse(
    SEX_EFFECTS$FWER_PBS[SEX_EFFECTS$Trait == PHENO] < 0.05,
    paste0(
      'Mean ', Info_Table$Abbreviation[Info_Table$Phenotype == PHENO],
      ' is significantly (p = ',
      signif(SEX_EFFECTS$PValue_PBS[SEX_EFFECTS$Trait == PHENO], digits = 3),
      ', FWER = ',
      signif(SEX_EFFECTS$FWER_PBS[SEX_EFFECTS$Trait == PHENO], digits = 3),
      ') ',
      ifelse(
        SEX_EFFECTS$Effect[SEX_EFFECTS$Trait == PHENO] < 0, 'lower', 'higher'
      ),
      ' in males by ',
      signif(SEX_EFFECTS$Effect[SEX_EFFECTS$Trait == PHENO]/SEX_EFFECTS$TraitSD[SEX_EFFECTS$Trait == PHENO], digits = 3),
      ' SD.'
    ),
    paste0(
      'Mean ', Info_Table$Abbreviation[Info_Table$Phenotype == PHENO],
      ' is ',
      ifelse(
        SEX_EFFECTS$Effect[SEX_EFFECTS$Trait == PHENO] < 0, 'lower', 'higher'
      ),
      ' in males by ',
      signif(SEX_EFFECTS$Effect[SEX_EFFECTS$Trait == PHENO]/AGE_EFFECTS$TraitSD[AGE_EFFECTS$Trait == PHENO], digits = 3),
      ' SD, but the difference is not significant (p = ',
      signif(SEX_EFFECTS$PValue_PBS[SEX_EFFECTS$Trait == PHENO], digits = 3),
      ', FWER = ',
      signif(SEX_EFFECTS$FWER_PBS[SEX_EFFECTS$Trait == PHENO], digits = 3),
      ').'
    )
  )
  AGE_CAPTION <- ifelse(
    AGE_EFFECTS$FWER_PBS[AGE_EFFECTS$Trait == PHENO] < 0.05,
    paste0(
      Info_Table$Abbreviation[Info_Table$Phenotype == PHENO],
      ' is significantly (p = ',
      signif(AGE_EFFECTS$PValue_PBS[AGE_EFFECTS$Trait == PHENO], digits = 3),
      ', FWER = ',
      signif(AGE_EFFECTS$FWER_PBS[AGE_EFFECTS$Trait == PHENO], digits = 3),
      ') affected by age, changing at a rate of ',
      signif(AGE_EFFECTS$Effect[AGE_EFFECTS$Trait == PHENO]/AGE_EFFECTS$TraitSD[AGE_EFFECTS$Trait == PHENO], digits = 3),
      ' SD per month.'
    ),
    paste0(
      Info_Table$Abbreviation[Info_Table$Phenotype == PHENO],
      ' is not significantly (p = ',
      signif(AGE_EFFECTS$PValue_PBS[AGE_EFFECTS$Trait == PHENO], digits = 3),
      ', FWER = ',
      signif(AGE_EFFECTS$FWER_PBS[AGE_EFFECTS$Trait == PHENO], digits = 3),
      ') affected by age.'
    )
  )
  AGESEX_CAPTION <- ifelse(
    AGESEX_INTERACTION_EFFECTS$FWER_PBS[AGESEX_INTERACTION_EFFECTS$Trait == PHENO] < 0.05,
    paste0(
      'The effect of age on ',
      Info_Table$Abbreviation[Info_Table$Phenotype == PHENO],
      ' differs significantly by sex (p = ',
      signif(AGESEX_INTERACTION_EFFECTS$FWER_PBS[AGESEX_INTERACTION_EFFECTS$Trait == PHENO], digits = 3),
      ', FWER = ',
      signif(AGESEX_INTERACTION_EFFECTS$PValue_PBS[AGESEX_INTERACTION_EFFECTS$Trait == PHENO], digits = 3),
      '). In females ',
      Info_Table$Abbreviation[Info_Table$Phenotype == PHENO],
      ' changes with age at a rate of ',
      signif(AGESEX_INTERACTION_EFFECTS$FemaleEffect[AGESEX_INTERACTION_EFFECTS$Trait == PHENO]/AGESEX_INTERACTION_EFFECTS$TraitSD[AGESEX_INTERACTION_EFFECTS$Trait == PHENO], digits = 3),
      ' SD per month, and in males it changes at a rate of ',
      signif(AGESEX_INTERACTION_EFFECTS$MaleEffect[AGESEX_INTERACTION_EFFECTS$Trait == PHENO]/AGESEX_INTERACTION_EFFECTS$TraitSD[AGESEX_INTERACTION_EFFECTS$Trait == PHENO], digits = 3),
      ' SD per month.'
    ),
    ''
  )
  CAPTION_TEXT <- paste(SD_CAPTION, SEX_CAPTION, AGE_CAPTION, AGESEX_CAPTION)
  
  i <- 0
  for(INDEX in seq(150, nchar(CAPTION_TEXT), by = 150)){
    POSITION <- INDEX + i
    if(substr(CAPTION_TEXT, POSITION, POSITION) != ' '){
      for(j in 1:20){
        if(substr(CAPTION_TEXT, POSITION - j, POSITION - j) == ' '){
          POSITION <- POSITION - j
          break
        }
        if(substr(CAPTION_TEXT, POSITION + j, POSITION + j) == ' '){
          POSITION <- POSITION + j
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
  
  rm(SD_CAPTION, SEX_CAPTION, AGE_CAPTION, AGESEX_CAPTION, CAPTION_TEXT)
  
  PLOT <- plot_grid(
    TITLE,
    plot_grid(
      AGE_PLOT + theme(legend.position = 'none'),
      COMBINED_SURV_CURVES,
      ncol = 2
    ),
    CAPTION,
    nrow = 3,
    rel_heights = c(0.1, 1, 0.15)
  )
  
  rm(AGE_PLOT, COMBINED_SURV_CURVES, TITLE, CAPTION)
  
  plot(PLOT)
  
  rm(PLOT)
  
  dev.off()
}

rm(PHENO)

#####


################################################################################
##  ################################

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