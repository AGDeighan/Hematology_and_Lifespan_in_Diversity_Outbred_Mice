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
  diet <- as.character(diet)
  
  DATA <- data.frame(
    mouse = mouse,
    pheno = pheno,
    tp = tp,
    diet = diet,
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

adjust_pheno_for_batch_and_clumping <- function(pheno, batch, clumps, tp, diet, mouse, outmult){
  
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
  diet <- as.character(diet)
  
  DATA <- data.frame(
    mouse = mouse,
    pheno = pheno,
    tp = tp,
    diet = diet,
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

load('data/original/genetic/JAC_longitudinal_genoprobs_20180618.Rdata')
REVISED_DATA <- ORIG_DATA <- read_csv('data/original/phenotypic/JAC_DO_HemaLS_Phenotypes.csv')
DICTIONARY <- read_csv('data/original/phenotypic/JAC_DO_HemaLS_Phenotypes_descr.csv')

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
    DOD = Death.date,
    Lifespan = lifespan,
    DeathType = Death.type
  )

#####


################################################################################
##### filtering #####

# We will remove all mice that died before 8 months OR missed their first CBC.
# We will also remove all mice that do not have genotype data. Also, since there
# are only three mice in this subset that are censored (two missing their cause
# of death and one missing its death date), we will remove all censored mice so
# we are focusing only on mice with full lifespan data.

GENOMICE <- rownames(probs[[1]])

REVISED_DATA <- REVISED_DATA %>% 
  filter(
    MouseID %in% GENOMICE,
    Lifespan/30.4 > 8,
    !is.na(hema.date.6),
    Died
  ) %>% 
  select(-Died)

#####


################################################################################
##### Make long by TP and bring in test date info #####

REVISED_DATA <- REVISED_DATA %>% 
  gather(key = K, value = V, contains('.'), -contains('date')) %>% 
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
  spread(key = K, value = V) %>% 
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
    MouseID, Sex, Cage, Generation, 
    DOB, DOD, Lifespan,
    Timepoint, DOT, AgeAtTest, everything()
  ) %>% 
  arrange(MouseID) %>% 
  arrange(MouseID, Timepoint)



# No missing dates
REVISED_DATA %>% filter(is.na(DOT), !is.na(ch))
REVISED_DATA %>% filter(is.na(DOT), !is.na(rbc))
REVISED_DATA %>% filter(is.na(DOT), !is.na(wbc))
REVISED_DATA %>% filter(is.na(DOT), !is.na(plt))
#####


################################################################################
##### Calculate NLR and convert HDW from SD to CV #####

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
  filter(DOT > DOD)
REVISED_DATA %>% 
  filter(AgeAtTest > Lifespan)

#####


################################################################################
##### Missing values #####

# # All mice that have their 6-month CBC data have all metrics
# 
# # 2 mice that have their 12-month CBC data are missing the eosinophil count
# 
# # 1 Mouse with 18-month CBC data is missing all the WBC metrics and 2 mice
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
# 1        rbc  7 Months        0.00       0        526   526
# 2         ch  7 Months        0.00       0        526   526
# 3       chcm  7 Months        0.00       0        526   526
# 4        hdw  7 Months        0.00       0        526   526
# 5        mcv  7 Months        0.00       0        526   526
# 6        rdw  7 Months        0.00       0        526   526
# 7        hct  7 Months        0.00       0        526   526
# 8        hgb  7 Months        0.00       0        526   526
# 9        wbc  7 Months        0.00       0        526   526
# 10       lnr  7 Months        0.00       0        526   526
# 11   n.lymph  7 Months        0.00       0        526   526
# 12    n.neut  7 Months        0.00       0        526   526
# 13    n.mono  7 Months        0.00       0        526   526
# 14     n.eos  7 Months        0.00       0        526   526
# 15       plt  7 Months        0.00       0        526   526
# 16       mpv  7 Months        0.00       0        526   526
# 17       mpm  7 Months        0.00       0        526   526
# 18       rbc 13 Months        0.00       0        504   504
# 19        ch 13 Months        0.00       0        504   504
# 20      chcm 13 Months        0.00       0        504   504
# 21       hdw 13 Months        0.00       0        504   504
# 22       mcv 13 Months        0.00       0        504   504
# 23       rdw 13 Months        0.00       0        504   504
# 24       hct 13 Months        0.00       0        504   504
# 25       hgb 13 Months        0.00       0        504   504
# 26       wbc 13 Months        0.00       0        504   504
# 27       lnr 13 Months        0.00       0        504   504
# 28   n.lymph 13 Months        0.00       0        504   504
# 29    n.neut 13 Months        0.00       0        504   504
# 30    n.mono 13 Months        0.00       0        504   504
# 31     n.eos 13 Months        0.40       2        502   504
# 32       plt 13 Months        0.00       0        504   504
# 33       mpv 13 Months        0.00       0        504   504
# 34       mpm 13 Months        0.00       0        504   504
# 35       rbc 19 Months        0.00       0        443   443
# 36        ch 19 Months        0.00       0        443   443
# 37      chcm 19 Months        0.00       0        443   443
# 38       hdw 19 Months        0.00       0        443   443
# 39       mcv 19 Months        0.00       0        443   443
# 40       rdw 19 Months        0.00       0        443   443
# 41       hct 19 Months        0.00       0        443   443
# 42       hgb 19 Months        0.00       0        443   443
# 43       wbc 19 Months        0.00       0        443   443
# 44       lnr 19 Months        0.23       1        442   443
# 45   n.lymph 19 Months        0.23       1        442   443
# 46    n.neut 19 Months        0.23       1        442   443
# 47    n.mono 19 Months        0.23       1        442   443
# 48     n.eos 19 Months        0.68       3        440   443
# 49       plt 19 Months        0.23       1        442   443
# 50       mpv 19 Months        0.23       1        442   443
# 51       mpm 19 Months        0.23       1        442   443

#####


################################################################################
##### Clumping effects #####

# About half the samples have flags for too many platelet clumps. The proportion
# of flagged samples is higher for the later timepoints, especially 18 months.
# The samples are flagged if the clump count is above 300. Some variables are 
# associated with the clumping variables
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
  mutate(
    # hct = (rbc * mcv) / 10,
    # wbc = n.lymph + n.neut + n.mono + n.eos,
    # nlr = n.neut/n.lymph
  ) %>% 
  select(
    MouseID, Sex, Generation, Lifespan, Timepoint, DOT, AgeAtTest,
    all_of(CBCPhenos), clumps, pltclm,
  )
#   rbc
# Power transform with power = 2
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/rbc.csv
# -------------------- #
#   
#   
#   ch
# Power transform with power = 2
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/ch.csv
# -------------------- #
#   
#   
#   chcm
# Power transform with power = 3
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/chcm.csv
# -------------------- #
#   
#   
#   hdw
# Power transform with power = 0
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/hdw.csv
# -------------------- #
#   
#   
#   mcv
# Power transform with power = 0
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/mcv.csv
# -------------------- #
#   
#   
#   rdw
# Power transform with power = 0
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/rdw.csv
# -------------------- #
#   
#   
#   hct
# Power transform with power = 2
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/hct.csv
# -------------------- #
#   
#   
#   hgb
# Power transform with power = 3
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/hgb.csv
# -------------------- #
#   
#   
#   wbc
# Power transform with power = 0
# boundary (singular) fit: see ?isSingular
# Estimated variance of one of the REF terms is near zero
# Groups   Name        Std.Dev.
# mouse    (Intercept) 0.70805 
# batch    (Intercept) 0.24478 
# tp       (Intercept) 0.00000 
# Residual             0.55469 
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/wbc.csv
# -------------------- #
#   
#   
#   nlr
# Power transform with power = 0
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/nlr.csv
# -------------------- #
#   
#   
#   n.lymph
# Power transform with power = 0
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/n.lymph.csv
# -------------------- #
#   
#   
#   n.neut
# Power transform with power = 0
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/n.neut.csv
# -------------------- #
#   
#   
#   n.mono
# Power transform with power = 0
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/n.mono.csv
# -------------------- #
#   
#   
#   n.eos
# Power transform with power = 0
# boundary (singular) fit: see ?isSingular
# Estimated variance of one of the REF terms is near zero
# Groups   Name        Std.Dev.
# mouse    (Intercept) 0.56070 
# batch    (Intercept) 0.14654 
# tp       (Intercept) 0.00000 
# Residual             0.58869 
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/n.eos.csv
# -------------------- #
#   
#   
#   plt
# Power transform with power = 0.5
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/plt.csv
# -------------------- #
#   
#   
#   mpv
# Power transform with power = 0
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/mpv.csv
# -------------------- #
#   
#   
#   mpm
# Power transform with power = 0
# Saving adjustment info to:
#   data/original/phenotypic/trait_adjustments/mpm.csv
# -------------------- #

ADJUSTED_DATA %>% filter_at(vars(all_of(CBCPhenos)), any_vars(. < 0))
# A tibble: 0 x 26

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
        MouseID, Sex, Generation, DOT, Timepoint, Adj
      ) %>% 
      full_join(
        DATA_RAW %>% 
          filter(!is.na(Raw)) %>% 
          select(
            MouseID, Timepoint, Raw
          ),
        by = c('MouseID', 'Timepoint')
      )
    rm(DATA_ADJ, DATA_RAW)
    
    DATA %>% filter(is.na(Adj) | is.na(Raw)) %>% print()
    
    SCATTER_PLOT <- DATA %>% 
      ggplot() +
      theme_minimal() +
      geom_point(
        aes(x = Raw, y = Adj, color = Timepoint)
      ) +
      theme(
        legend.position = 'top'
      ) +
      labs(
        title = toupper(pheno)
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
        legend.position = 'none'
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
    
    
    PLOT <- cowplot::plot_grid(
      cowplot::plot_grid(
        SCATTER_PLOT, HISTOGRAMS, TP_BOXPLOTS,
        nrow = 1,
        rel_widths = c(1, 1.25, 1.5)
      ),
      BATCH_BOXPLOTS,
      ncol = 1,
      rel_heights = c(1, 2)
    )
    
    plot(PLOT)
  }
  dev.off()
  rm(pheno, DATA, SCATTER_PLOT, HISTOGRAMS, TP_BOXPLOTS, BATCH_BOXPLOTS, PLOT)
  
  
  pdf(
    'figures/data_prep/distributions/RawAndAdj_TraitDistributions_BestPowerTransform.pdf',
    width = 14, height = 21
  )
  for(pheno in CBCPhenos){
    
    print(pheno)
    
    PT <- optimal_pt_minimize_skewness(REVISED_DATA[[pheno]])
    print(paste0("Power transform of ", PT))
    
    
    DATA_ADJ <- ADJUSTED_DATA 
    DATA_RAW <- REVISED_DATA
    
    DATA_ADJ[['Adj']] <- DATA_ADJ[[pheno]] %>% power_transform(power = PT)
    DATA_RAW[['Raw']] <- DATA_RAW[[pheno]] %>% power_transform(power = PT)
    
    DATA <- DATA_ADJ %>% 
      filter(!is.na(Adj)) %>% 
      select(
        MouseID, Sex, Generation, DOT, Timepoint, Adj
      ) %>% 
      full_join(
        DATA_RAW %>% 
          filter(!is.na(Raw)) %>% 
          select(
            MouseID, Timepoint, Raw
          ),
        by = c('MouseID', 'Timepoint')
      )
    rm(DATA_ADJ, DATA_RAW)
    
    DATA %>% filter(is.na(Adj) | is.na(Raw)) %>% print()
    
    SCATTER_PLOT <- DATA %>% 
      ggplot() +
      theme_minimal() +
      geom_point(
        aes(x = Raw, y = Adj, color = Timepoint)
      ) +
      theme(
        legend.position = 'top'
      ) +
      labs(
        title = paste0(toupper(pheno), ' with power transform of ', PT)
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
        legend.position = 'none'
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
    
    
    PLOT <- cowplot::plot_grid(
      cowplot::plot_grid(
        SCATTER_PLOT, HISTOGRAMS, TP_BOXPLOTS,
        nrow = 1,
        rel_widths = c(1, 1.25, 1.5)
      ),
      BATCH_BOXPLOTS,
      ncol = 1,
      rel_heights = c(1, 2)
    )
    
    plot(PLOT)
  }
  dev.off()
  rm(pheno, DATA, SCATTER_PLOT, HISTOGRAMS, TP_BOXPLOTS, BATCH_BOXPLOTS, PLOT)
}

#######


################################################################################
##### Clumping effects after regressing out batch effects #####

if(CREATE_PLOTS){
  pdf(
    'figures/data_prep/clumping_effects/Platelet_Clumping_Effects_after_Adjustment.pdf',
    width = 12, height = 12
  )
  for(pheno in c('n.eos', 'plt', 'mpv', 'mpm')){
    clumping_effects_plots(ADJUSTED_DATA, pheno)
  }
  dev.off()
  rm(pheno)
}

if(CREATE_PLOTS){
  pdf(
    'figures/data_prep/clumping_effects/Platelet_Clumping_Effects_after_Adjustment_lnPheno.pdf',
    width = 12, height = 12
  )
  for(pheno in c('n.eos', 'plt', 'mpv', 'mpm')){
    clumping_effects_plots(
      ADJUSTED_DATA %>% 
        mutate_at(vars(all_of(pheno)), log), 
      pheno
    )
  }
  dev.off()
  rm(pheno)
}

#####



################################################################################
##### Order columns #####

ADJUSTED_DATA <- ADJUSTED_DATA %>% 
  select(
    MouseID, Generation, Sex, 
    Lifespan, Timepoint, AgeAtTest, 
    all_of(CBCPhenos)
  )

#####



################################################################################
##### saving #####

ADJUSTED_DATA %>% 
  filter(
    !is.na(AgeAtTest)     # Remove rows missing all data
  ) %>% 
  mutate(
    Sex = ifelse(Sex == 'F', 'Female', 'Male')
  ) %>% 
  write_csv(
    'data/working/phenotypic/ClumpDateAdj_HemaLS.csv'
  )

#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####