# 2021-01-28
################################################################################
#
#   This script cleans the hematology data, creates distribution plots, and
#   creates the hematology datasets for use in following analyses
#
#
#   Author:  Andrew Deighan
#   E-mail: andrew.deighan@jax.org
#
################################################################################

CREATE_PLOTS = TRUE

################################################################################
## libraries etc #########################################################
options(na.action = 'na.exclude')
options(stringsAsFactors = FALSE)
library(tidyverse)


################################################################################
## helper functions #########################################################

mdy_Date <- function(x){
  
  if(length(x) == 1){
    if(is.na(x)){
      return(as.Date(x))
    }
    test <- str_split(x, '/', simplify = TRUE)[3]
    if(nchar(test) == 4){
      x <- as.Date(x, format = '%m/%d/%Y')
    } else{
      x <- as.Date(x, format = '%m/%d/%y')
    }
    return(x)
  } else{
    f <- function(y){
      if(is.na(y)){
        return(as.character(y))
      }
      test <- str_split(y, '/', simplify = TRUE)[3]
      if(nchar(test) == 4){
        y <- as.character(as.Date(y, format = '%m/%d/%Y'))
      } else{
        y <- as.character(as.Date(y, format = '%m/%d/%y'))
      }
      return(y)
    }
    
    x <- unlist(lapply(x, f))
    
    x <- as.Date(x, format = '%Y-%m-%d')
    
    return(x)
  }
  
  
}

clumping_effects_plots <- function(data,
                                   feature, 
                                   sex = 'Sex',
                                   timepoint = 'Timepoint',
                                   sex_levels = c('F', 'M'),
                                   return_plotlist = FALSE,
                                   rankZ = FALSE){
  require(dplyr)
  require(ggplot2)
  require(grid)
  require(gridExtra)
  
  REVISED_DATA <- data
  FEATURE <- feature
  
  REVISED_DATA[['FEATURE']] <- REVISED_DATA[[FEATURE]]
  REVISED_DATA[['Sex']] <- REVISED_DATA[[sex]]
  REVISED_DATA[['Timepoint']] <- REVISED_DATA[[timepoint]]
  
  REVISED_DATA <- REVISED_DATA %>% 
    filter(!is.na(FEATURE),
           !is.na(pltclm),
           !is.na(clumps))
  
  log_tf <- function(x, adj = 0.01, ...){
    low <- min(x, na.rm = TRUE)
    if(low <= 0){
      x <- log(x - low + adj, ...)
    } else{
      x <- log(x, ...)
    }
    return(x)
  }
  
  if(rankZ){
    rankZ = function(x) {
      x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
      return(qnorm(x))
    }
    
    REVISED_DATA <- REVISED_DATA %>% 
      mutate(FEATURE = rankZ(FEATURE))
  }
  
  FULL_MODEL <- lm(FEATURE ~ Sex + Generation + Timepoint + pltclm, data = REVISED_DATA)
  NULL_MODEL <- lm(FEATURE ~ Sex + Generation + Timepoint, data = REVISED_DATA )
  ANOVA <- anova(NULL_MODEL, FULL_MODEL)
  PVAL <- ANOVA$`Pr(>F)`[2]
  
  
  P1 <- REVISED_DATA %>% 
    mutate(pltclm = factor(pltclm)) %>% 
    filter(!is.na(pltclm)) %>% 
    ggplot(aes(x = pltclm, y = FEATURE, color = Sex)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(x = jitter(as.numeric(pltclm)) + if_else(Sex == sex_levels[1], -.2, .2), 
                   shape = Timepoint),
               alpha = 0.5) +
    labs(title = paste0('Effect of clumping flag on ', FEATURE),
         subtitle = paste0('P-value of effect of clumping flag ~ ', signif(PVAL, 3)),
         y = FEATURE)
  
  
  FULL_MODEL <- lm(FEATURE ~ Sex + Generation + Timepoint + log_tf(clumps), data = REVISED_DATA)
  NULL_MODEL <- lm(FEATURE ~ Sex + Generation + Timepoint, data = REVISED_DATA )
  ANOVA <- anova(NULL_MODEL, FULL_MODEL)
  PVAL <- ANOVA$`Pr(>F)`[2]
  
  P2 <- REVISED_DATA %>% 
    mutate(pltclm = factor(pltclm)) %>% 
    filter(!is.na(pltclm)) %>% 
    ggplot(aes(x = log_tf(clumps), y = FEATURE, color = Sex, shape = Timepoint)) +
    geom_point(alpha = 0.5) +
    labs(title = paste0('Effect of number of platelet clumps on ', FEATURE),
         subtitle = paste0('P-value of effect of number of platelet clumps ~ ', signif(PVAL, 3)),
         y = FEATURE)
  
  
  
  FULL_MODEL <- lm(FEATURE ~ Sex + Generation + Timepoint + pltclm + log_tf(clumps), data = REVISED_DATA)
  NULL_MODEL <- lm(FEATURE ~ Sex + Generation + Timepoint + pltclm, data = REVISED_DATA )
  ANOVA <- anova(NULL_MODEL, FULL_MODEL)
  PVAL_CLUMPS <- ANOVA$`Pr(>F)`[2]
  
  FULL_MODEL <- lm(FEATURE ~ Sex + Generation + Timepoint + log_tf(clumps) + pltclm, data = REVISED_DATA)
  NULL_MODEL <- lm(FEATURE ~ Sex + Generation + Timepoint + log_tf(clumps), data = REVISED_DATA )
  ANOVA <- anova(NULL_MODEL, FULL_MODEL)
  PVAL_FLAG <- ANOVA$`Pr(>F)`[2]
  
  FULL_MODEL <- lm(FEATURE ~ Sex + Generation + Timepoint + log_tf(clumps) + pltclm + pltclm:log_tf(clumps), data = REVISED_DATA)
  NULL_MODEL <- lm(FEATURE ~ Sex + Generation + Timepoint + log_tf(clumps) + pltclm, data = REVISED_DATA )
  ANOVA <- anova(NULL_MODEL, FULL_MODEL)
  PVAL_INT <- ANOVA$`Pr(>F)`[2]
  
  P3 <- REVISED_DATA %>% 
    mutate(
      pltclm = ifelse(
        pltclm == 1, 'Clumping Flag',
        'No Flag'
      ),
      pltclm = factor(
        pltclm,
        levels = c('No Flag', 'Clumping Flag')
      )
    ) %>% 
    filter(!is.na(pltclm)) %>% 
    ggplot(aes(x = log_tf(clumps), y = FEATURE, color = Sex, shape = Timepoint)) +
    geom_point(alpha = 0.5) +
    facet_grid(. ~ pltclm, scales = 'free_x') +
    labs(title = paste0('Effect of platelet clumping on ', FEATURE),
         subtitle = paste0('P-value flag after accounting for clump count ~ ', signif(PVAL_FLAG, 3),
                           '\nP-value clump count after accounting for flag ~ ', signif(PVAL_CLUMPS, 3),
                           '\nP-value of count:flag interaction effect ~ ', signif(PVAL_INT, 3)),
         y = FEATURE)
  
  if(return_plotlist){
    return(list(P1, P2, P3))
  }
  
  grid.arrange(P1, P2, P3,
               layout_matrix = matrix(c(1,2,3,3), nrow = 2, byrow = TRUE))
  
}

standardise <- function(x){
  return(
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  )
}

ab_scale <- function(x, a = 0.01, b = 0.99){
  # This function rescales a vector of values to a specified range
  
  MIN <- min(x, na.rm = TRUE)
  RANGE <- max(x, na.rm = TRUE) - MIN
  
  OUT <- ((x - MIN)/(RANGE)) * (b - a) + a
  return(OUT)
}

flip <- function(x){
  return(max(x, na.rm = TRUE) - x + min(x, na.rm = TRUE))
}

log_tf <- function(x, shift = 0.01){
  if(min(x, na.rm = TRUE) <= 0){
    x <- x - min(x, na.rm = TRUE) + shift
  }
  return(log(x))
}

best_transform <- function(pheno, shift = 0.01){
  SKEW <- e1071::skewness(pheno, na.rm = TRUE)
  KURT <- e1071::kurtosis(pheno, na.rm = TRUE)
  
  if(
    SKEW > 0 &
    abs(e1071::skewness(
      log_tf(x = pheno, shift = shift), na.rm = TRUE
    )) < abs(SKEW) &
    abs(e1071::kurtosis(
      log_tf(x = pheno, shift = shift), na.rm = TRUE
    )) < abs(KURT)
  ){
    TRANSFORM <- 'log'
    pheno <- log_tf(x = pheno, shift = shift)
    TMIN <- min(pheno, na.rm = TRUE); TMAX <- max(pheno, na.rm = TRUE)
  } else if(
    SKEW < 0 &
    abs(e1071::skewness(
      log_tf(x = flip(pheno), shift = shift), na.rm = TRUE
    )) < abs(SKEW) &
    abs(e1071::kurtosis(
      log_tf(x = flip(pheno), shift = shift), na.rm = TRUE
    )) < abs(KURT)
  ){
    TRANSFORM <- 'flipped log'
    pheno <- log_tf(x = flip(pheno), shift = shift)
    TMIN <- min(pheno, na.rm = TRUE); TMAX <- max(pheno, na.rm = TRUE)
  } else{
    TRANSFORM <- TYPE <- 'none'
  }
  
  return(TRANSFORM)
}

adjust_pheno_for_batch <- function(pheno, batch, tp, mouse, sex){
  
  # This function uses a fully random effects model to account estimate batch
  # effects conditioned on individual and timepoint, and then removes the
  # estimated batch effects.
  
  ORIG_PHENO <- pheno
  PHENO_NAME <- deparse(substitute(pheno))
  message(PHENO_NAME)
  
  # We will save the original min and max of the data so that we can 
  # reset the adjusted values to the same range as the original values.
  ORIGINAL_MIN <- min(pheno, na.rm = TRUE)
  ORIGINAL_MAX <- max(pheno, na.rm = TRUE)
  
  
  
  # We will use a log tranform if reduces skewness and kurtosis.
  shift <- 0.01
  SKEW <- e1071::skewness(pheno, na.rm = TRUE)
  KURT <- e1071::kurtosis(pheno, na.rm = TRUE)
  
  if(
    SKEW > 0 &
    abs(e1071::skewness(
      log_tf(x = pheno, shift = shift), na.rm = TRUE
    )) < abs(SKEW) &
    abs(e1071::kurtosis(
      log_tf(x = pheno, shift = shift), na.rm = TRUE
    )) < abs(KURT)
  ){
    TRANSFORM <- 'log'
    message('log transform')
    pheno <- log_tf(x = pheno, shift = shift)
    TMIN <- min(pheno, na.rm = TRUE); TMAX <- max(pheno, na.rm = TRUE)
  } else if(
    SKEW < 0 &
    abs(e1071::skewness(
      log_tf(x = flip(pheno), shift = shift), na.rm = TRUE
    )) < abs(SKEW) &
    abs(e1071::kurtosis(
      log_tf(x = flip(pheno), shift = shift), na.rm = TRUE
    )) < abs(KURT)
  ){
    TRANSFORM <- 'flipped log'
    message('flipped log transform')
    pheno <- log_tf(x = flip(pheno), shift = shift)
    TMIN <- min(pheno, na.rm = TRUE); TMAX <- max(pheno, na.rm = TRUE)
  } else{
    TRANSFORM <- TYPE <- 'none'
    message('No transform')
  }
  
  # Standarised (mean = 0, sd = 1) data for model fitting to reduce risk 
  # of convergence errors. The mean and SD used to standardise the data 
  # must be saved so that the data can be un-standardised afterwards
  MEAN <- mean(pheno, na.rm = TRUE)
  SD <- sd(pheno, na.rm = TRUE)
  pheno <- (pheno - MEAN)/SD
  
  # Convert batch and tp to characters so that they are treated as a categorical 
  # variables
  batch <- as.character(batch)
  tp <- as.character(tp)
  sex <- as.character(sex)
  
  DATA <- data.frame(
    mouse = mouse,
    pheno = pheno,
    tp = tp,
    sex = sex,
    batch = batch
  )
  
  # Fit model to estimate batch effects
  MODEL <- lme4::lmer(
    pheno ~ (1|tp) + (1|sex) + (1|mouse) + (1|batch), 
    control = lme4::lmerControl(                        # If you get convergence errors
      optCtrl = list(                                   # you may be able to prevent them
        xtol_abs = 1e-16,                               # by making these tolerances more
        ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
      )                                                 # 1e-6)
    ),
    data = DATA
  )
  
  VC <- lme4::VarCorr(MODEL)
  if(sum(VC == 0) > 0 | lme4::isSingular(MODEL)){
    message('Estimated variance of one of the REF terms is near zero')
    print(VC)
  }
  
  # Pull estimated batch effects from the model
  REF <- lme4::ranef(MODEL)
  BATCH <- REF$batch
  BATCH <- BATCH[batch,]
  
  # Remove estimated batch effects
  PHENO <- pheno - BATCH
  
  DATA <- DATA %>% 
    mutate(
      adj_transformed = PHENO,
      BatchEffect = BATCH
    )
  
  # Un-standardise data
  PHENO <- PHENO * SD + MEAN
  
  # Un-transform if a transformation was made
  if(TRANSFORM == 'log'){
    PHENO <- ab_scale(PHENO, a = TMIN, b = TMAX)
    PHENO <- exp(PHENO)
  } else if(TRANSFORM == 'flipped log'){
    PHENO <- ab_scale(PHENO, a = TMIN, b = TMAX)
    PHENO <- exp(PHENO)
    PHENO <- flip(PHENO)
  }
  
  # reset to original range
  PHENO <- ab_scale(PHENO, a = ORIGINAL_MIN, b = ORIGINAL_MAX)
  
  DATA <- DATA %>% 
    rename(
      MouseID = mouse,
      raw_transformed = pheno,
      Timepoint = tp,
      Batch = batch
    ) %>% 
    mutate(
      Transform = TRANSFORM,
      adj_untransformed = PHENO,
      raw_untransformed = ORIG_PHENO
    ) %>% 
    select(
      MouseID, Timepoint, Batch, BatchEffect,
      raw_untransformed, raw_transformed,
      adj_untransformed, adj_transformed
    )
  
  FILENAME <- paste0(
    'data/raw/phenotypic/trait_adjustments/',
    PHENO_NAME,
    '.csv'
  )
  
  message(paste0(
    'Saving adjustment info to:\n',
    FILENAME
  ))
  
  write_csv(
    DATA, FILENAME
  )
  
  message('-------------------- #\n\n')
  
  return(PHENO)
}

adjust_pheno_for_batch_and_clumping <- function(pheno, batch, clumps, tp, mouse, sex){
  
  # This function uses a fully random effects model to account estimate batch
  # effects conditioned on individual and timepoint, and then removes the
  # estimated batch effects.
  
  ORIG_PHENO <- pheno
  PHENO_NAME <- deparse(substitute(pheno))
  message(PHENO_NAME)
  
  # We will save the original min and max of the data so that we can 
  # reset the adjusted values to the same range as the original values.
  ORIGINAL_MIN <- min(pheno, na.rm = TRUE)
  ORIGINAL_MAX <- max(pheno, na.rm = TRUE)
  
  
  
  # We will use a log tranform if reduces skewness and kurtosis.
  shift <- 0.01
  SKEW <- e1071::skewness(pheno, na.rm = TRUE)
  KURT <- e1071::kurtosis(pheno, na.rm = TRUE)
  
  if(
    SKEW > 0 &
    abs(e1071::skewness(
      log_tf(x = pheno, shift = shift), na.rm = TRUE
    )) < abs(SKEW) &
    abs(e1071::kurtosis(
      log_tf(x = pheno, shift = shift), na.rm = TRUE
    )) < abs(KURT)
  ){
    TRANSFORM <- 'log'
    message('log transform')
    pheno <- log_tf(x = pheno, shift = shift)
    TMIN <- min(pheno, na.rm = TRUE); TMAX <- max(pheno, na.rm = TRUE)
  } else if(
    SKEW < 0 &
    abs(e1071::skewness(
      log_tf(x = flip(pheno), shift = shift), na.rm = TRUE
    )) < abs(SKEW) &
    abs(e1071::kurtosis(
      log_tf(x = flip(pheno), shift = shift), na.rm = TRUE
    )) < abs(KURT)
  ){
    TRANSFORM <- 'flipped log'
    message('flipped log transform')
    pheno <- log_tf(x = (ORIGINAL_MAX - pheno)/ORIGINAL_MAX, shift = shift)
    TMIN <- min(pheno, na.rm = TRUE); TMAX <- max(pheno, na.rm = TRUE)
  } else{
    TRANSFORM <- TYPE <- 'none'
    message('No transform')
  }
  
  
  # Standarised (mean = 0, sd = 1) data for model fitting to reduce risk 
  # of convergence errors. The mean and SD used to standardise the data 
  # must be saved so that the data can be un-standardised afterwards
  MEAN <- mean(pheno, na.rm = TRUE)
  SD <- sd(pheno, na.rm = TRUE)
  pheno <- (pheno - MEAN)/SD
  
  # Convert batch and tp to characters so that they are treated as a categorical 
  # variables
  batch <- as.character(batch)
  tp <- as.character(tp)
  sex <- as.character(sex)
  
  DATA <- data.frame(
    mouse = mouse,
    pheno = pheno,
    tp = tp,
    sex = sex,
    batch = batch,
    clumps = clumps
  )
  
  # Fit model to estimate batch effects
  MODEL <- lme4::lmer(
    pheno ~ (1|tp) + (1|sex) + (1|mouse) + (1|batch) + clumps,
    control = lme4::lmerControl(                        # If you get convergence errors
      optCtrl = list(                                   # you may be able to prevent them
        xtol_abs = 1e-16,                               # by making these tolerances more
        ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
      )                                                 # 1e-6)
    ),
    data = DATA       
  )
  
  VC <- lme4::VarCorr(MODEL)
  if(sum(VC == 0) > 0 | lme4::isSingular(MODEL)){
    message('Estimated variance of one of the REF terms is near zero')
    print(VC)
  }
  
  # Pull clump coefficient from model
  COEF <- lme4::fixef(MODEL)
  COEF <- COEF['clumps']
  
  # Pull estimated batch effects from the model
  REF <- lme4::ranef(MODEL)
  BATCH <- REF$batch
  BATCH <- BATCH[batch,]
  
  # Remove estimated batch effects
  PHENO <- pheno - BATCH - clumps*COEF
  
  DATA <- DATA %>% 
    mutate(
      adj_transformed = PHENO,
      BatchEffect = BATCH,
      ClumpCoef = COEF
    )
  
  # Un-standardise data
  PHENO <- PHENO * SD + MEAN
  
  # Un-transform if a transformation was made
  if(TRANSFORM == 'log'){
    PHENO <- ab_scale(PHENO, a = TMIN, b = TMAX)
    PHENO <- exp(PHENO)
  } else if(TRANSFORM == 'flipped log'){
    PHENO <- ab_scale(PHENO, a = TMIN, b = TMAX)
    PHENO <- exp(PHENO)
    PHENO <- flip(PHENO)
  }
  
  # reset to original range
  PHENO <- ab_scale(PHENO, a = ORIGINAL_MIN, b = ORIGINAL_MAX)
  
  DATA <- DATA %>% 
    rename(
      MouseID = mouse,
      raw_transformed = pheno,
      Timepoint = tp,
      Batch = batch,
      Clumps = clumps
    ) %>% 
    mutate(
      Transform = TRANSFORM,
      adj_untransformed = PHENO,
      raw_untransformed = ORIG_PHENO
    ) %>% 
    select(
      MouseID, Timepoint, 
      Batch, BatchEffect, Clumps, ClumpCoef,
      raw_untransformed, raw_transformed,
      adj_untransformed, adj_transformed
    )
  
  FILENAME <- paste0(
    'data/raw/phenotypic/trait_adjustments/',
    PHENO_NAME,
    '.csv'
  )
  
  message(paste0(
    'Saving adjustment info to:\n',
    FILENAME
  ))
  
  write_csv(
    DATA, FILENAME
  )
  
  message('-------------------- #\n\n')
  
  return(PHENO)
}

#####


################################################################################
## load data etc #########################################################

GENOPROBS <- readRDS(
  'data/raw/genetic/Long_AlleleProbs.rds'
)
REVISED_DATA <- ORIG_DATA <- read_csv(
  'data/raw/phenotypic/JAC_DO_HemaLS_Phenotypes.csv'
)
DICTIONARY <- read_csv(
  'data/raw/phenotypic/JAC_DO_HemaLS_Phenotypes_descr.csv'
)

#####

RBCphenos <- c(
  'rbc', 'ch', 'chcm', 'hdw',
  'mcv', 'rdw', 
  'hct', 'hgb'
)
WBCphenos <- c(
  'wbc', 'nlr',
  'n.lymph', 
  'n.neut',
  'n.mono',
  'n.eos'
)
Platelets <- c('plt', 'mpv', 'mpm')
CBCPhenos <- c(RBCphenos, WBCphenos, Platelets)




################################################################################
##### Rename Columns #####

REVISED_DATA <- REVISED_DATA %>% 
  rename(
    MouseID = Mouse.ID,
    Cage = cage,
    Coat = Coat.Color,
    DOD = Death.date,
    Lifespan = lifespan,
    DeathType = Death.type
  )

#####


################################################################################
##### Reformat mouse ID names to match those in genoprobs #####

# The mouse IDs in the phenotype data are of the format DO.XXXX while in the
# genoprobs they are DO-XXXX

REVISED_DATA <- REVISED_DATA %>% 
  mutate(
    MouseID = gsub('\\.', '-', MouseID)
  )

#####


################################################################################
##### filtering #####

# We will remove all mice that died before 8 months OR missed their first CBC.
# We will also remove all mice that do not have genotype data. Also, since there
# are only three mice in this subset that are censored (two missing their cause
# of death and one missing its death date), we will remove all censored mice so
# we are focusing only on mice with full lifespan data.

GENOMICE <- rownames(GENOPROBS[[1]])

REVISED_DATA <- REVISED_DATA %>% 
  filter(
    MouseID %in% GENOMICE,
    Lifespan/30.4 > 8,
    !is.na(hema.date.6),
    Died
  ) %>% 
  select(-Died, -DeathType)

rm(GENOMICE, GENOPROBS)

#####


################################################################################
##### Make long by TP and bring in test date info #####

REVISED_DATA <- REVISED_DATA %>% 
  pivot_longer(
    names_to = 'K',
    values_to = 'V',
    cols = c(contains('.'), -contains('date'))
  ) %>% 
  mutate(
    Timepoint = ifelse(
      substr(K, nchar(K), 1000) == '6',  '7 Months',
      ifelse(
        substr(K, nchar(K) - 1, 1000) == '12', '13 Months',
        '19 Months'
        )
    ),
    Timepoint = factor(
      Timepoint, 
      levels = c('7 Months', '13 Months', '19 Months')
    ),
    K = ifelse(
      substr(K, nchar(K), 1000) == '6', 
      substr(K, 1, nchar(K) - 2),
      substr(K, 1, nchar(K) - 3)
    )
  ) %>% 
  pivot_wider(
    names_from = K,
    values_from = V
  ) %>% 
  mutate(
    DOT = if_else(
      Timepoint == '7 Months', hema.date.6,
      if_else(
        Timepoint == '13 Months', hema.date.12,
        hema.date.18
      )
    ),
    DOB = mdy_Date(DOB),
    DOD = mdy_Date(DOD),
    DOT = mdy_Date(DOT),
    AgeAtTest = as.numeric(difftime(DOT, DOB, units = 'days'))/30.4,
    Lifespan = Lifespan/30.4
  ) %>% 
  select(-hema.date.6, -hema.date.12, -hema.date.18) %>% 
  select(
    MouseID, Sex, Coat, Cage, Generation, 
    DOB, DOD, Lifespan,
    Timepoint, DOT, AgeAtTest, everything()
  ) %>% 
  arrange(MouseID) %>% 
  arrange(MouseID, Timepoint)



# No missing dates
REVISED_DATA %>% 
  filter(is.na(DOT)) %>% 
  filter_at(vars(wbc:n.eos), any_vars(!is.na(.))) %>% 
  nrow()
# [1] 0

#####


################################################################################
##### Calculate NLR and convert HDW from SD to CV #####

# The neutrophil to lynphocyte ratio is a biologically relevant measure (it has 
# been associated with various pathological states). The red cell distribution
# width (RDW) is reported as the coefficient of variation (which is the common
# metric) in red blood cell volume, but the hemoglobin distribution width (HDW)
# is reported in standard deviations. To make them more comparable, we will 
# convert HDW to a coefficient of variation in red blood cell hemoglobin 
# concentration

REVISED_DATA <- REVISED_DATA %>% 
  mutate(
    hdw = (hdw/chcm) * 100,
    nlr = n.neut/n.lymph
  ) 

#####


################################################################################
##### Check for impossible test dates #####

# None

REVISED_DATA %>% 
  filter(DOT > DOD) %>% 
  nrow()
# [1] 0

REVISED_DATA %>% 
  filter(AgeAtTest > Lifespan) %>% 
  nrow()
# [1] 0

#####


################################################################################
##### Missing values #####

# # All mice that have their 7-month CBC data have all metrics
# 
# # 2 mice that have their 13-month CBC data are missing the eosinophil count
# 
# # 1 Mouse with 19-month CBC data is missing all the WBC metrics and 2 mice
# # are missing just the eosinophils.

REVISED_DATA %>% 
  filter(!is.na(AgeAtTest)) %>% 
  select(Sex, Generation, Timepoint, all_of(CBCPhenos)) %>% 
  gather(key = Phenotype, value = IsMissing, -Sex, -Generation, -Timepoint) %>% 
  mutate(IsMissing = is.na(IsMissing)) %>% 
  group_by(Phenotype, Timepoint) %>% 
  summarise(
    PercMissing = mean(IsMissing) * 100,
    Missing = sum(IsMissing),
    NotMissing = sum(!IsMissing),
    Total = n()
  ) %>% 
  ungroup() %>% 
  mutate(
    Phenotype = factor(Phenotype, levels = CBCPhenos),
    PercMissing = round(PercMissing, digits = 2)
  ) %>% 
  arrange(Timepoint, Phenotype) %>% 
  data.frame()
#    Phenotype Timepoint PercMissing Missing NotMissing Total
# 1        rbc  7 Months        0.00       0        521   521
# 2         ch  7 Months        0.00       0        521   521
# 3       chcm  7 Months        0.00       0        521   521
# 4        hdw  7 Months        0.00       0        521   521
# 5        mcv  7 Months        0.00       0        521   521
# 6        rdw  7 Months        0.00       0        521   521
# 7        hct  7 Months        0.00       0        521   521
# 8        hgb  7 Months        0.00       0        521   521
# 9        wbc  7 Months        0.00       0        521   521
# 10       nlr  7 Months        0.00       0        521   521
# 11   n.lymph  7 Months        0.00       0        521   521
# 12    n.neut  7 Months        0.00       0        521   521
# 13    n.mono  7 Months        0.00       0        521   521
# 14     n.eos  7 Months        0.00       0        521   521
# 15       plt  7 Months        0.00       0        521   521
# 16       mpv  7 Months        0.00       0        521   521
# 17       mpm  7 Months        0.00       0        521   521
# 18       rbc 13 Months        0.00       0        499   499
# 19        ch 13 Months        0.00       0        499   499
# 20      chcm 13 Months        0.00       0        499   499
# 21       hdw 13 Months        0.00       0        499   499
# 22       mcv 13 Months        0.00       0        499   499
# 23       rdw 13 Months        0.00       0        499   499
# 24       hct 13 Months        0.00       0        499   499
# 25       hgb 13 Months        0.00       0        499   499
# 26       wbc 13 Months        0.00       0        499   499
# 27       nlr 13 Months        0.00       0        499   499
# 28   n.lymph 13 Months        0.00       0        499   499
# 29    n.neut 13 Months        0.00       0        499   499
# 30    n.mono 13 Months        0.00       0        499   499
# 31     n.eos 13 Months        0.40       2        497   499
# 32       plt 13 Months        0.00       0        499   499
# 33       mpv 13 Months        0.00       0        499   499
# 34       mpm 13 Months        0.00       0        499   499
# 35       rbc 19 Months        0.00       0        440   440
# 36        ch 19 Months        0.00       0        440   440
# 37      chcm 19 Months        0.00       0        440   440
# 38       hdw 19 Months        0.00       0        440   440
# 39       mcv 19 Months        0.00       0        440   440
# 40       rdw 19 Months        0.00       0        440   440
# 41       hct 19 Months        0.00       0        440   440
# 42       hgb 19 Months        0.00       0        440   440
# 43       wbc 19 Months        0.00       0        440   440
# 44       nlr 19 Months        0.23       1        439   440
# 45   n.lymph 19 Months        0.23       1        439   440
# 46    n.neut 19 Months        0.23       1        439   440
# 47    n.mono 19 Months        0.23       1        439   440
# 48     n.eos 19 Months        0.68       3        437   440
# 49       plt 19 Months        0.23       1        439   440
# 50       mpv 19 Months        0.23       1        439   440
# 51       mpm 19 Months        0.23       1        439   440

#####


################################################################################
##### Clumping effects #####

# About half the samples have flags for too many platelet clumps. The proportion
# of flagged samples is higher for the later timepoints, especially 19 months.
# The samples are flagged if the clump count is above 300. Eosinophil counts 
# and all platelet metrics are known to be affected by platelet clumping.

REVISED_DATA %>%  filter(pltclm == 0) %>% select(clumps) %>% range()
# [1]  20 300
REVISED_DATA %>%  filter(pltclm == 1) %>% select(clumps) %>% range()
# [1]   301 52327


if(CREATE_PLOTS){
  pdf(
    'figures/data_prep/clumping_effects/Platelet_Clumping_Effects_before_Adjustment.pdf',
    width = 12, height = 12
  )
  for(pheno in c('n.eos', 'plt', 'mpv', 'mpm')){
    clumping_effects_plots(REVISED_DATA, pheno)
  }
  dev.off()
  rm(pheno)
}

if(CREATE_PLOTS){
  pdf(
    'figures/data_prep/clumping_effects/Platelet_Clumping_Effects_before_Adjustment_lnPheno.pdf',
    width = 12, height = 12
  )
  for(pheno in c('n.eos', 'plt', 'mpv', 'mpm')){
    clumping_effects_plots(
      REVISED_DATA %>% 
        mutate_at(vars(all_of(pheno)), log), 
      pheno
    )
  }
  dev.off()
  rm(pheno)
}

#####


################################################################################
##### Adjusting for date and clumping #####

TRANSFORM_TABLE <- REVISED_DATA %>%
  select(
    MouseID, Timepoint, wbc:n.eos, nlr
  ) %>% 
  pivot_longer(
    cols = wbc:nlr,
    names_to = 'Phenotype',
    values_to = 'X'
  ) %>% 
  group_by(Phenotype) %>% 
  summarise(
    Transform = best_transform(X)
  ) %>% 
  ungroup() %>% 
  mutate(
    Phenotype = factor(
      Phenotype,
      levels = c(
        "rbc", "ch", "chcm", "hdw", "mcv", "rdw", "hct", "hgb", 
        "wbc", "nlr", "n.lymph", "n.neut", "n.mono", "n.eos",
        "plt", "mpv", "mpm" 
      )
    )
  ) %>% 
  arrange(Phenotype)

ADJUSTED_DATA <- REVISED_DATA %>% 
  mutate(
    rbc = adjust_pheno_for_batch(pheno = rbc,  batch = DOT, tp = Timepoint, mouse = MouseID, sex = Sex),
    ch = adjust_pheno_for_batch(pheno = ch,  batch = DOT, tp = Timepoint, mouse = MouseID, sex = Sex),
    chcm = adjust_pheno_for_batch(pheno = chcm,  batch = DOT, tp = Timepoint, mouse = MouseID, sex = Sex),
    hdw = adjust_pheno_for_batch(pheno = hdw,  batch = DOT, tp = Timepoint, mouse = MouseID, sex = Sex),
    mcv = adjust_pheno_for_batch(pheno = mcv,  batch = DOT, tp = Timepoint, mouse = MouseID, sex = Sex),
    rdw = adjust_pheno_for_batch(pheno = rdw,  batch = DOT, tp = Timepoint, mouse = MouseID, sex = Sex),
    hct = adjust_pheno_for_batch(pheno = hct,  batch = DOT, tp = Timepoint, mouse = MouseID, sex = Sex),
    hgb = adjust_pheno_for_batch(pheno = hgb,  batch = DOT, tp = Timepoint, mouse = MouseID, sex = Sex),
    wbc = adjust_pheno_for_batch(pheno = wbc,  batch = DOT, tp = Timepoint, mouse = MouseID, sex = Sex),
    nlr = adjust_pheno_for_batch(pheno = nlr,  batch = DOT, tp = Timepoint, mouse = MouseID, sex = Sex),
    n.lymph = adjust_pheno_for_batch(pheno = n.lymph,  batch = DOT, tp = Timepoint, mouse = MouseID, sex = Sex),
    n.neut = adjust_pheno_for_batch(pheno = n.neut,  batch = DOT, tp = Timepoint, mouse = MouseID, sex = Sex),
    n.mono = adjust_pheno_for_batch(pheno = n.mono,  batch = DOT, tp = Timepoint, mouse = MouseID, sex = Sex),
    n.eos = adjust_pheno_for_batch_and_clumping(pheno = n.eos,  batch = DOT, clumps = log(clumps), tp = Timepoint, mouse = MouseID, sex = Sex),
    plt = adjust_pheno_for_batch_and_clumping(pheno = plt,  batch = DOT, clumps = log(clumps), tp = Timepoint, mouse = MouseID, sex = Sex),
    mpv = adjust_pheno_for_batch_and_clumping(pheno = mpv,  batch = DOT, clumps = log(clumps), tp = Timepoint, mouse = MouseID, sex = Sex),
    mpm = adjust_pheno_for_batch_and_clumping(pheno = mpm,  batch = DOT, clumps = log(clumps), tp = Timepoint, mouse = MouseID, sex = Sex)
  ) %>% 
  select(
    MouseID, Sex, Coat, Generation, Lifespan, Timepoint, DOT, AgeAtTest,
    all_of(CBCPhenos), clumps, pltclm,
  )
#   rbc
# No transform
# boundary (singular) fit: see ?isSingular
# Estimated variance of one of the REF terms is near zero
# Groups   Name        Std.Dev.  
# mouse    (Intercept) 0.61123697
# batch    (Intercept) 0.21590170
# tp       (Intercept) 0.14703694
# sex      (Intercept) 0.00001189
# Residual             0.76163180
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/rbc.csv
# -------------------- #
#   
#   
#   ch
# flipped log transform
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/ch.csv
# -------------------- #
#   
#   
#   chcm
# flipped log transform
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/chcm.csv
# -------------------- #
#   
#   
#   hdw
# log transform
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/hdw.csv
# -------------------- #
#   
#   
#   mcv
# log transform
# boundary (singular) fit: see ?isSingular
# Estimated variance of one of the REF terms is near zero
# Groups   Name        Std.Dev.
# mouse    (Intercept) 0.55104 
# batch    (Intercept) 0.50978 
# tp       (Intercept) 0.49495 
# sex      (Intercept) 0.00000 
# Residual             0.49585 
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/mcv.csv
# -------------------- #
#   
#   
#   rdw
# log transform
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/rdw.csv
# -------------------- #
#   
#   
#   hct
# flipped log transform
# boundary (singular) fit: see ?isSingular
# Estimated variance of one of the REF terms is near zero
# Groups   Name        Std.Dev.  
# mouse    (Intercept) 4.4412e-01
# batch    (Intercept) 5.5078e-01
# tp       (Intercept) 5.8212e-01
# sex      (Intercept) 1.3064e-05
# Residual             5.1521e-01
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/hct.csv
# -------------------- #
#   
#   
#   hgb
# flipped log transform
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/hgb.csv
# -------------------- #
#   
#   
#   wbc
# log transform
# boundary (singular) fit: see ?isSingular
# Estimated variance of one of the REF terms is near zero
# Groups   Name        Std.Dev.
# mouse    (Intercept) 0.70948 
# batch    (Intercept) 0.23539 
# tp       (Intercept) 0.00000 
# sex      (Intercept) 0.11614 
# Residual             0.66542 
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/wbc.csv
# -------------------- #
#   
#   
#   nlr
# log transform
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/nlr.csv
# -------------------- #
#   
#   
#   n.lymph
# log transform
# boundary (singular) fit: see ?isSingular
# Estimated variance of one of the REF terms is near zero
# Groups   Name        Std.Dev.
# mouse    (Intercept) 0.66235 
# batch    (Intercept) 0.21523 
# tp       (Intercept) 0.12768 
# sex      (Intercept) 0.00000 
# Residual             0.71413 
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/n.lymph.csv
# -------------------- #
#   
#   
#   n.neut
# log transform
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/n.neut.csv
# -------------------- #
#   
#   
#   n.mono
# log transform
# Warning message:
#   ℹ Model failed to converge with max|grad| = 0.00370358 (tol = 0.002, component 1)
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/n.mono.csv
# -------------------- #
#   
#   
#   n.eos
# log transform
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/n.eos.csv
# -------------------- #
#   
#   
#   plt
# No transform
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/plt.csv
# -------------------- #
#   
#   
#   mpv
# No transform
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/mpv.csv
# -------------------- #
#   
#   
#   mpm
# log transform
# boundary (singular) fit: see ?isSingular
# Estimated variance of one of the REF terms is near zero
# Groups   Name        Std.Dev.
# mouse    (Intercept) 0.68632 
# batch    (Intercept) 0.43650 
# tp       (Intercept) 0.14083 
# sex      (Intercept) 0.00000 
# Residual             0.56956 
# Saving adjustment info to:
#   data/raw/phenotypic/trait_adjustments/mpm.csv
# -------------------- #

ADJUSTED_DATA %>% 
  filter_at(vars(all_of(CBCPhenos)), any_vars(. < 0)) %>% 
  nrow()
# [1] 0

CLUMPING_TRAITS <- c('n.eos', Platelets)

if(CREATE_PLOTS){
  pdf(
    'figures/data_prep/distributions/RawAndAdj_TraitDistributions_Untransformed.pdf',
    width = 14, height = 21
  )
  for(pheno in CBCPhenos){
    
    print(pheno)
    
    DATA_ADJ <- ADJUSTED_DATA 
    DATA_RAW <- REVISED_DATA
    
    DATA_ADJ[['Adj']] <- DATA_ADJ[[pheno]]
    DATA_RAW[['Raw']] <- DATA_RAW[[pheno]]
    
    DATA <- DATA_ADJ %>% 
      filter(!is.na(Adj)) %>% 
      select(
        MouseID, Sex, Generation, DOT, Timepoint, AgeAtTest, clumps, Adj
      ) %>% 
      full_join(
        DATA_RAW %>% 
          filter(!is.na(Raw)) %>% 
          select(
            MouseID, Timepoint, Raw
          ),
        by = c('MouseID', 'Timepoint')
      ) %>% 
      mutate(
        DOT = as.character(DOT),
        AgeAtTest = standardise(AgeAtTest),
        clumps = standardise(log(clumps))
      )
    rm(DATA_ADJ, DATA_RAW)
    
    # Calculate p-values
    if(pheno %in% CLUMPING_TRAITS){
      # Raw data
      RAW_AGE <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Raw ~ AgeAtTest + clumps + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          ),
          lme4::lmer(
            Raw ~ clumps + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      RAW_TP <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Raw ~ (1|Timepoint) + clumps + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          ),
          lme4::lmer(
            Raw ~ clumps + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      RAW_SEX <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Raw ~ AgeAtTest + clumps + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          ),
          lme4::lmer(
            Raw ~ AgeAtTest + clumps + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      RAW_SEXAGE <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Raw ~ AgeAtTest + clumps + (1 + AgeAtTest|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          ),
          lme4::lmer(
            Raw ~ AgeAtTest + clumps + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      RAW_BATCH <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Raw ~ (1|Timepoint) + clumps + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          ),
          lme4::lmer(
            Raw ~ (1|Timepoint) + clumps + (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      RAW_CLUMPS <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Raw ~ (1|Timepoint) + clumps + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          ),
          lme4::lmer(
            Raw ~ (1|Timepoint) + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      
      
      
      # Adjusted data
      ADJ_AGE <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Adj ~ AgeAtTest + clumps + (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          ),
          lme4::lmer(
            Adj ~ clumps + (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      ADJ_TP <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Adj ~ (1|Timepoint) + clumps + (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          ),
          lme4::lmer(
            Adj ~ clumps + (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      ADJ_SEX <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Adj ~ AgeAtTest + clumps + (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          ),
          lme4::lmer(
            Adj ~ AgeAtTest + clumps + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      ADJ_SEXAGE <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Adj ~ AgeAtTest + clumps + (1 + AgeAtTest|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          ),
          lme4::lmer(
            Adj ~ AgeAtTest + clumps + (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      ADJ_BATCH <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Adj ~ (1|Timepoint) + clumps + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          ),
          lme4::lmer(
            Adj ~ (1|Timepoint) + clumps + (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      ADJ_CLUMPS <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Adj ~ (1|Timepoint) + clumps + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          ),
          lme4::lmer(
            Adj ~ (1|Timepoint) + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      
      PVALUE_STRING <- paste0(
        'P-values for unadjusted values:\n    ',
        'Timepoint ~ ', RAW_TP, '; Age ~ ', RAW_AGE, '; Sex ~ ', RAW_SEX, '; Age:Sex ~ ', RAW_SEXAGE, '; Batch ~ ', RAW_BATCH, '; Clumps ~ ', RAW_CLUMPS, '\n',
        'P-values for adjusted values:\n    ',
        'Timepoint ~ ', ADJ_TP, '; Age ~ ', ADJ_AGE, '; Sex ~ ', ADJ_SEX, '; Age:Sex ~ ', ADJ_SEXAGE, '; Batch ~ ', ADJ_BATCH, '; Clumps ~ ', ADJ_CLUMPS
      )
      
      rm(
        RAW_AGE, RAW_BATCH, RAW_SEX, RAW_SEXAGE, RAW_TP, RAW_CLUMPS,
        ADJ_AGE, ADJ_BATCH, ADJ_SEX, ADJ_SEXAGE, ADJ_TP, ADJ_CLUMPS
      )
    } else{
      # Raw data
      RAW_AGE <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Raw ~ AgeAtTest + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          ),
          lme4::lmer(
            Raw ~ (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      RAW_TP <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Raw ~ (1|Timepoint) + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          ),
          lme4::lmer(
            Raw ~ (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      RAW_SEX <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Raw ~ AgeAtTest + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          ),
          lme4::lmer(
            Raw ~ AgeAtTest + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      RAW_SEXAGE <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Raw ~ AgeAtTest + (1 + AgeAtTest|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          ),
          lme4::lmer(
            Raw ~ AgeAtTest + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      RAW_BATCH <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Raw ~ (1|Timepoint) + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          ),
          lme4::lmer(
            Raw ~ (1|Timepoint) + (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Raw = standardise(Raw))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      
      
      
      # Adjusted data
      ADJ_AGE <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Adj ~ AgeAtTest + (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          ),
          lme4::lmer(
            Adj ~ (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      ADJ_TP <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Adj ~ (1|Timepoint) + (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          ),
          lme4::lmer(
            Adj ~ (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      ADJ_SEX <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Adj ~ AgeAtTest + (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          ),
          lme4::lmer(
            Adj ~ AgeAtTest + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      ADJ_SEXAGE <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Adj ~ AgeAtTest + (1 + AgeAtTest|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          ),
          lme4::lmer(
            Adj ~ AgeAtTest + (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      ADJ_BATCH <- signif(
        suppressMessages(anova(
          lme4::lmer(
            Adj ~ (1|Timepoint) + (1|Sex) + (1|Generation) + (1|MouseID) + (1|DOT),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          ),
          lme4::lmer(
            Adj ~ (1|Timepoint) + (1|Sex) + (1|Generation) + (1|MouseID),
            REML = FALSE,
            control = lme4::lmerControl(                        # If you get convergence errors
              optCtrl = list(                                   # you may be able to prevent them
                xtol_abs = 1e-16,                               # by making these tolerances more
                ftol_abs = 1e-16                                # strict (e.g. 1e-8, the defult is
              )                                                 # 1e-6)
            ),
            data = DATA %>% 
              mutate(Adj = standardise(Adj))
          )
        )$`Pr(>Chisq)`[2]),
        digits = 3
      )
      
      PVALUE_STRING <- paste0(
        'P-values for unadjusted values:\n    ',
        'Timepoint ~ ', RAW_TP, '; Age ~ ', RAW_AGE, '; Sex ~ ', RAW_SEX, '; Age:Sex ~ ', RAW_SEXAGE, '; Batch ~ ', RAW_BATCH, '\n',
        'P-values for adjusted values:\n    ',
        'Timepoint ~ ', ADJ_TP, '; Age ~ ', ADJ_AGE, '; Sex ~ ', ADJ_SEX, '; Age:Sex ~ ', ADJ_SEXAGE, '; Batch ~ ', ADJ_BATCH
      )
      
      rm(
        RAW_AGE, RAW_BATCH, RAW_SEX, RAW_SEXAGE, RAW_TP,
        ADJ_AGE, ADJ_BATCH, ADJ_SEX, ADJ_SEXAGE, ADJ_TP
      )
    }
    
    SCATTER_PLOT <- DATA %>% 
      ggplot() +
      theme_minimal() +
      geom_point(
        aes(x = Raw, y = Adj, color = Timepoint, shape = Sex)
      ) +
      theme(
        legend.position = 'none'
      ) +
      labs(
        title = '\n'
      )
    
    HISTOGRAMS <- DATA %>% 
      gather(key = Data, value = X, Adj, Raw) %>% 
      mutate(
        Data = factor(Data, levels = c('Raw', 'Adj'))
      ) %>% 
      ggplot() +
      theme_minimal() +
      geom_histogram(
        aes(x = X, fill = Timepoint),
        bins = 100
      ) +
      facet_grid(
        Data ~ .
      ) +
      theme(
        legend.position = 'top'
      ) +
      labs(
        x = pheno
      )
    
    TP_BOXPLOTS <- DATA %>% 
      gather(key = Data, value = X, Adj, Raw) %>% 
      mutate(
        Data = factor(Data, levels = c('Raw', 'Adj'))
      ) %>% 
      ggplot() +
      theme_minimal() +
      geom_boxplot(
        aes(x = Timepoint, y = X, color = Timepoint),
        outlier.shape = NA
      ) +
      geom_jitter(
        aes(x = Timepoint, y = X, color = Timepoint, shape = Sex),
        alpha = 0.5
      ) +
      facet_grid(
        Data ~ .
      ) +
      guides(
        color = FALSE
      ) +
      theme(
        legend.position = 'top'
      ) +
      labs(
        y = pheno
      )
    
    BATCH_BOXPLOTS <- DATA %>% 
      gather(key = Data, value = X, Adj, Raw) %>% 
      mutate(
        Data = factor(Data, levels = c('Raw', 'Adj')),
        Date = as.character(DOT)
      ) %>% 
      ggplot() +
      theme_minimal() +
      geom_boxplot(
        aes(x = Date, y = X, color = Timepoint),
        outlier.shape = NA
      ) +
      geom_jitter(
        aes(x = Date, y = X, color = Timepoint, shape = Sex),
        alpha = 0.5
      ) +
      facet_wrap(
        Data ~ Timepoint,
        scales = 'free_x',
        dir = 'v',
        ncol = 2
      ) +
      theme(
        axis.text.x = element_text(angle = 45),
        legend.position = 'none'
      ) +
      labs(
        y = pheno
      )
    
    TITLE <- ggplot() + 
      labs(
        title = toupper(pheno),
        subtitle = PVALUE_STRING
      ) + 
      theme_minimal()
    
    rm(PVALUE_STRING)
    
    PLOT <- cowplot::plot_grid(
      TITLE,
      cowplot::plot_grid(
        cowplot::plot_grid(
          SCATTER_PLOT, HISTOGRAMS, TP_BOXPLOTS,
          nrow = 1,
          rel_widths = c(1, 1, 1.25)
        ),
        BATCH_BOXPLOTS,
        ncol = 1,
        rel_heights = c(1, 2)
      ),
      ncol = 1,
      rel_heights = c(0.1, 1)
    )
    
    plot(PLOT)
  }
  dev.off()
  rm(pheno, DATA, SCATTER_PLOT, HISTOGRAMS, TP_BOXPLOTS, BATCH_BOXPLOTS, TITLE, PLOT)
}

rm(CLUMPING_TRAITS)

#######


################################################################################
##### Order columns #####

ADJUSTED_DATA <- ADJUSTED_DATA %>% 
  filter(
    !is.na(AgeAtTest)     # Remove rows missing all data
  ) %>% 
  mutate(
    Sex = ifelse(Sex == 'F', 'Female', 'Male')
  ) %>% 
  select(
    MouseID, Generation, Sex, Coat,
    Lifespan, Timepoint, AgeAtTest, 
    all_of(CBCPhenos)
  )

RAW_DATA <- REVISED_DATA %>% 
  filter(
    !is.na(AgeAtTest)     # Remove rows missing all data
  ) %>% 
  mutate(
    Sex = ifelse(Sex == 'F', 'Female', 'Male')
  ) %>% 
  select(
    all_of(names(ADJUSTED_DATA)),
    DOT, clumps
  )

#####


################################################################################
##### saving #####

TRANSFORM_TABLE %>% 
  write_csv(
    'data/processed/phenotypic/Transform_Table.csv'
  )

ADJUSTED_DATA %>% 
  write_csv(
    'data/processed/phenotypic/ClumpDateAdj_HemaLS.csv'
  )
# # A tibble: 1,460 x 24
#    MouseID Generation Sex    Coat   Lifespan Timepoint AgeAtTest   rbc    ch  chcm   hdw   mcv   rdw   hct   hgb   wbc    nlr n.lymph n.neut n.mono n.eos   plt   mpv   mpm
#    <chr>   <chr>      <chr>  <chr>     <dbl> <fct>         <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>   <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>
#  1 DO-0001 G7         Female agouti     37.9 7 Months       7.14  11.5  14.4  30.9  6.03  47.4  11.9  55.8  16.5  5.60 0.0917    4.82  0.457  0.101 0.371 1222.  7.08  1.16
#  2 DO-0001 G7         Female agouti     37.9 13 Months     13.1   11.4  14.0  32.0  6.03  44.3  12.2  51.8  16.3 12.1  0.131    10.2   1.37   0.319 0.632 1057.  6.09  1.12
#  3 DO-0001 G7         Female agouti     37.9 19 Months     19.1   10.7  14.1  31.3  5.76  45.9  11.8  49.1  16.0  7.66 0.182     6.28  1.10   0.224 0.742 1216.  6.00  1.09
#  4 DO-0002 G7         Female agouti     11.9 7 Months       7.14  11.6  14.6  31.6  5.84  46.7  13.7  55.1  16.8  6.66 0.419     4.40  1.82   0.206 0.612 1537.  6.37  1.14
#  5 DO-0003 G7         Female agouti     29.8 7 Months       7.07  10.9  15.1  31.5  5.66  48.9  13.9  54.3  16.6  8.10 0.111     7.02  0.796  0.111 0.467 1407.  6.91  1.12
#  6 DO-0003 G7         Female agouti     29.8 13 Months     13.1   10.6  15.4  32.8  5.78  47.4  14.2  51.1  16.4  9.79 0.429     6.85  2.90   0.190 0.394 1538.  6.45  1.11
#  7 DO-0003 G7         Female agouti     29.8 19 Months     19.0   10.9  14.9  32.8  5.80  46.3  13.1  50.3  16.5  6.24 0.359     4.62  1.58   0.198 0.469 1218.  6.04  1.08
#  8 DO-0004 G7         Female agouti     15.3 7 Months       7.07  11.0  14.9  29.8  6.11  51.9  13.4  59.1  16.4  6.03 0.0753    5.15  0.403  0.121 0.392 1727.  6.85  1.23
#  9 DO-0004 G7         Female agouti     15.3 13 Months     13.1   10.9  14.3  31.4  6.19  46.4  12.9  51.5  15.6  8.24 0.114     6.89  0.806  0.179 0.644 1613.  5.96  1.21
# 10 DO-0005 G7         Female agouti     30.2 7 Months       7.17  11.3  14.4  30.0  6.34  49.5  13.3  57.8  16.1  8.42 0.108     7.31  0.806  0.195 0.251 1299.  7.00  1.19
# # … with 1,450 more rows

RAW_DATA %>% 
  write_csv(
    'data/processed/phenotypic/Raw_HemaLS.csv'
  )
# # A tibble: 1,460 x 26
#    MouseID Generation Sex    Coat   Lifespan Timepoint AgeAtTest   rbc    ch  chcm   hdw   mcv   rdw   hct   hgb   wbc    nlr n.lymph n.neut n.mono n.eos   plt   mpv   mpm DOT        clumps
#    <chr>   <chr>      <chr>  <chr>     <dbl> <fct>         <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>   <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <date>      <dbl>
#  1 DO-0001 G7         Female agouti     37.9 7 Months       7.14  11.3  14.5  27.3  6.19  53    11.9  59.9  16.6  5.26 0.0946    4.44   0.42   0.1   0.28  1116   7.6  1.18 2012-01-25    900
#  2 DO-0001 G7         Female agouti     37.9 13 Months     13.1   11.0  14.1  29.7  5.86  47.7  12.4  52.3  16.2 13.5  0.132    11.0    1.46   0.28  0.57  1007   8.1  1.21 2012-07-25   1031
#  3 DO-0001 G7         Female agouti     37.9 19 Months     19.1   10.8  14    30.9  5.86  45.4  11.8  49.2  16.2  9.46 0.201     7      1.41   0.27  0.69  1038   4.6  1.13 2013-01-22   1288
#  4 DO-0002 G7         Female agouti     11.9 7 Months       7.14  11.4  14.7  28.1  6.01  52.3  13.7  59.4  16.9  6.24 0.425     4.05   1.72   0.2   0.24  1691   6.9  1.16 2012-01-25     66
#  5 DO-0003 G7         Female agouti     29.8 7 Months       7.07  10.7  15.2  27.9  5.84  54.7  13.9  58.7  16.7  7.57 0.114     6.49   0.74   0.11  0.17  1605   7.4  1.14 2012-01-25     45
#  6 DO-0003 G7         Female agouti     29.8 13 Months     13.1   10.2  15.5  30.5  5.64  50.9  14.4  51.7  16.3 11.0  0.428     7.34   3.14   0.17  0.17  1807   8.4  1.2  2012-07-25     41
#  7 DO-0003 G7         Female agouti     29.8 19 Months     19.0   11.0  14.8  32.4  5.90  45.8  13.1  50.3  16.7  7.73 0.395     5.14   2.03   0.24  0.23  1331   4.6  1.12 2013-01-22     77
#  8 DO-0004 G7         Female agouti     15.3 7 Months       7.07  10.8  15    26    6.27  57.9  13.4  62.5  16.5  5.66 0.0779    4.75   0.37   0.12  0.38  1496   7.4  1.25 2012-01-25   2606
#  9 DO-0004 G7         Female agouti     15.3 13 Months     13.1   10.4  14.4  29    6     49.9  13.1  52.1  15.5  9.27 0.115     7.38   0.85   0.16  0.75  1435   8    1.31 2012-07-25   3046
# 10 DO-0005 G7         Female agouti     30.2 7 Months       7.17  11.1  14.5  26.2  6.49  55.3  13.3  61.5  16.2  7.87 0.111     6.76   0.75   0.19  0.11  1438   7.5  1.21 2012-01-25     82
# # … with 1,450 more rows

#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####