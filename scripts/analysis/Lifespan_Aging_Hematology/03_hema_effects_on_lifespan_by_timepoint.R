# 2021-05-27
################################################################################
#
#   This script estimates the effect of each of the CBC traits, across 
#   timepoints, on lifespan. 
# 
#   This script is designed to by run on a high performance computing cluster 
#   with at least 50 cores due to the multitude of simulations performed during 
#   the parametric bootstraps to estimate statistical significance of effects.
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
library(lme4)
library(pbkrtest)
library(parallel)

#####


################################################################################
## helper functions #########################################################

rankZ <- function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

standardise <- function(x){
  return(
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  )
}

war <- function(...){'warning'}

seqPBmodcomp <- function(largeModel, smallModel, h = 1, minsim = 10000, maxsim = 10^17, chunk.size = 10000, seed, ...) {
  # The seqPBmodcomp function in the pbkrtest package is not actually defined,
  # so we copied and modified the function from the original paper, "A 
  # Kenward-Roger Approximation and Parametric Bootstrap Methods for Tests in 
  # Linear Mixed Models The R Package pbkrtest"  
  set.seed(seed = seed)
  t.start <- proc.time()
  nchunk <- (maxsim %/% chunk.size) + 2
  LRTstat <- pbkrtest:::getLRT.merMod(largeModel, smallModel)
  
  chunk_seed <- sample.int(99999:1000000, size = 1)
  ref <- c(PBrefdist(largeModel, smallModel, nsim = minsim * 1.05, seed = chunk_seed, ...))
  n.extreme <- sum(ref > LRTstat["tobs"])
  
  if (n.extreme >= h & length(ref > 0) >= minsim){
    
    ans <- PBmodcomp(largeModel, smallModel, ref = ref)
    ans$ctime <- (proc.time() - t.start)[3]
    return(ans)
    
  } else{
    
    chunk_seed <- sample.int(99999:1000000, size = 1)
    ref <- c(ref, PBrefdist(largeModel, smallModel, nsim = minsim * 1.05, seed = chunk_seed, ...))
    n.extreme <- sum(ref > LRTstat["tobs"])
    
    if (n.extreme >= h & length(ref > 0) >= minsim){
      
      ans <- PBmodcomp(largeModel, smallModel, ref = ref)
      ans$ctime <- (proc.time() - t.start)[3]
      return(ans)
      
    } else{
      
      for (ii in 1:nchunk) {
        chunk_seed <- sample.int(99999:1000000, size = 1)
        ref <- c(ref, PBrefdist(largeModel, smallModel, nsim = chunk.size, seed = chunk_seed, ...))
        n.extreme <- sum(ref > LRTstat["tobs"])
        if (n.extreme >= h & length(ref > 0) >= minsim)
          break
      }
      
      ans <- PBmodcomp(largeModel, smallModel, ref = ref)
      ans$ctime <- (proc.time() - t.start)[3]
      return(ans)
    }
  }
}

#####


################################################################################
## load data etc #########################################################

load('data/processed/phenotypic/CBC_AnalysisEnvironment.RData')

TRANSFORM_TABLE <- read_csv(
  'data/processed/phenotypic/Transform_Table.csv'
)

#####

MIN_EXTREME <- 2
MIN_SIM <- 10000
MAX_SIM <- 10^7
CHUNK_SIZE <- 10^6
SEED <- 19940418
CLUSTER <- makeForkCluster(50)


################################################################################
## LS effects at 7 months #############################

LS_EFFECTS_07MO <- data.frame(
  Trait = vector(),
  Timepoint = vector(),
  PValue_LRT = vector(),
  PValue_PBS = vector(),
  # Total number of succesfull simulations (LRT statistic > 0)
  NSim = vector(),
  # Simulations in which the LRT statistic was greater than that observed
  NExtreme = vector(),
  BootstrapTime = vector(),
  Effect = vector(),
  SE = vector(),
  TraitSD = vector()
)

for(TRAIT in CBCphenos){
  cat(paste0(
    TRAIT,'\n'
  ))
  
  # Create dataset for model fitting. We will standardise (mean 0, sd = 1) all
  # continuous variables to improve model fitting
  ADJ_DATA <- PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == '7 Months')
  
  if(TRAIT %in% TRANSFORM_TABLE$Phenotype[TRANSFORM_TABLE$Transform == 'log']){
    ADJ_DATA[['X']] <- log(ADJ_DATA[[TRAIT]])
  } else{
    ADJ_DATA[['X']] <- ADJ_DATA[[TRAIT]]
  }
  
  
  TRAIT_SD <- sd(ADJ_DATA[['X']], na.rm = TRUE)
  
  ADJ_DATA <- ADJ_DATA %>%
    mutate(
      Lifespan = Lifespan - AgeAtTest,
      Lifespan = standardise(Lifespan),
      X = standardise(X)
    ) %>%
    select(
      MouseID, Generation, Sex, X, Lifespan
    )
  
  
  # Fit full model of raw data, but with the trait rank-Z transformed
  FULL_MODEL <- tryCatch(
    suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex + X,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa"),
      data = ADJ_DATA %>% mutate(X = rankZ(X))
    )),
    warning = war
  )
  # Fit null model of raw data, but with the trait rank-Z transformed
  NULL_MODEL <- tryCatch(
    suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa"),
      data = ADJ_DATA %>% mutate(X = rankZ(X))
    )),
    warning = war
  )
  
  if(is.character(FULL_MODEL) | is.character(NULL_MODEL)){
    cat('WARNING: Full and/or null model throws a warning, refitting with Nelder_Mead optimization\n')
    if(is.character(FULL_MODEL)){
      cat('Full model throws a warning\n')
    }
    if(is.character(NULL_MODEL)){
      cat('Null model throws a warning\n')
    }
    # Refit full and null models using Nelder_Mead optimizer
    FULL_MODEL <- suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex + X,
      REML = FALSE,
      control = lmerControl(optimizer = "Nelder_Mead"),
      data = ADJ_DATA %>% mutate(X = rankZ(X))
    ))
    NULL_MODEL <- suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex,
      REML = FALSE,
      control = lmerControl(optimizer = "Nelder_Mead"),
      data = ADJ_DATA %>% mutate(X = rankZ(X))
    ))
  }
  
  # Fit full model of raw data, without rank-Z transform of trait
  FULL_MODEL_NRZ <- tryCatch(
    suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex + X,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa"),
      data = ADJ_DATA
    )),
    warning = war
  )
  if(is.character(FULL_MODEL_NRZ)){
    cat('Non-rankZ full model throws a warning, refitting with Nelder_Mead optimization\n')
    FULL_MODEL_NRZ <- suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex + X,
      REML = FALSE,
      control = lmerControl(optimizer = "Nelder_Mead"),
      data = ADJ_DATA
    ))
    
  }
  
  
  TEST_RESULTS <- seqPBmodcomp(
    FULL_MODEL,
    NULL_MODEL,
    h = MIN_EXTREME,
    minsim = MIN_SIM,
    maxsim = MAX_SIM,
    chunk.size = CHUNK_SIZE,
    seed = SEED,
    cl = CLUSTER
  )
  
  LS_EFFECTS_07MO <- rbind(
    LS_EFFECTS_07MO,
    data.frame(
      Trait = TRAIT,
      Timepoint = '7 Months',
      PValue_LRT = as.numeric(TEST_RESULTS$test['LRT', 'p.value']),
      PValue_PBS = as.numeric(TEST_RESULTS$test['PBtest', 'p.value']),
      # Total number of simulations
      NSim = as.numeric(TEST_RESULTS$samples[2]),
      # Simulations in which the LRT statistic was greater than that observed
      NExtreme = as.numeric(TEST_RESULTS$n.extreme),
      BootstrapTime = as.numeric(TEST_RESULTS$ctime),
      Effect = as.numeric(fixef(FULL_MODEL_NRZ)['X']),
      SE = as.numeric(sqrt(vcov(FULL_MODEL_NRZ)['X', 'X'])),
      TraitSD = TRAIT_SD
    )
  )
  
  print(VarCorr(FULL_MODEL))
  
  cat('\n -------------------- # \n\n')
}
rm(TRAIT, TRAIT_SD, ADJ_DATA, FULL_MODEL, NULL_MODEL, FULL_MODEL_NRZ, TEST_RESULTS)


sd(
  PhenoData_ClumpDateAdj$Lifespan[PhenoData_ClumpDateAdj$Timepoint == '7 Months'] -
    PhenoData_ClumpDateAdj$AgeAtTest[PhenoData_ClumpDateAdj$Timepoint == '7 Months']
)
# 8.194393

LS_EFFECTS_07MO <- LS_EFFECTS_07MO %>%
  mutate(
    FWER_PBS = p.adjust(PValue_PBS, method = 'holm'),
    # The phenotypes and lifespan were standerdised during model fitting, so to 
    # get the effects and standard errors in units of lifespan (months) per 
    # trait, we need to  multiply them by the standard deviation of the lifespan
    # in months and divide by the standard deviation of the trait
    Effect = Effect * 8.194393 / TraitSD,
    SE = SE * 8.194393 / TraitSD
  ) %>%
  select(
    Trait, PValue_PBS, FWER_PBS, PValue_LRT,
    everything()
  )


#####


################################################################################
## LS effects at 13 months #############################

LS_EFFECTS_13MO <- data.frame(
  Trait = vector(),
  Timepoint = vector(),
  PValue_LRT = vector(),
  PValue_PBS = vector(),
  # Total number of succesfull simulations (LRT statistic > 0)
  NSim = vector(),
  # Simulations in which the LRT statistic was greater than that observed
  NExtreme = vector(),
  BootstrapTime = vector(),
  Effect = vector(),
  SE = vector(),
  TraitSD = vector()
)

for(TRAIT in CBCphenos){
  cat(paste0(
    TRAIT,'\n'
  ))
  
  # Create dataset for model fitting. We will standardise (mean 0, sd = 1) all
  # continuous variables to improve model fitting
  ADJ_DATA <- PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == '13 Months')
  
  if(TRAIT %in% TRANSFORM_TABLE$Phenotype[TRANSFORM_TABLE$Transform == 'log']){
    ADJ_DATA[['X']] <- log(ADJ_DATA[[TRAIT]])
  } else{
    ADJ_DATA[['X']] <- ADJ_DATA[[TRAIT]]
  }
  
  
  TRAIT_SD <- sd(ADJ_DATA[['X']], na.rm = TRUE)
  ADJ_DATA <- ADJ_DATA %>%
    mutate(
      Lifespan = Lifespan - AgeAtTest,
      Lifespan = standardise(Lifespan),
      X = standardise(X)
    ) %>%
    select(
      MouseID, Generation, Sex, X, Lifespan
    )
  
  
  # Fit full model of raw data, but with the trait rank-Z transformed
  FULL_MODEL <- tryCatch(
    suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex + X,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa"),
      data = ADJ_DATA %>% mutate(X = rankZ(X))
    )),
    warning = war
  )
  # Fit null model of raw data, but with the trait rank-Z transformed
  NULL_MODEL <- tryCatch(
    suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa"),
      data = ADJ_DATA %>% mutate(X = rankZ(X))
    )),
    warning = war
  )
  
  if(is.character(FULL_MODEL) | is.character(NULL_MODEL)){
    cat('WARNING: Full and/or null model throws a warning, refitting with Nelder_Mead optimization\n')
    if(is.character(FULL_MODEL)){
      cat('Full model throws a warning\n')
    }
    if(is.character(NULL_MODEL)){
      cat('Null model throws a warning\n')
    }
    # Refit full and null models using Nelder_Mead optimizer
    FULL_MODEL <- suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex + X,
      REML = FALSE,
      control = lmerControl(optimizer = "Nelder_Mead"),
      data = ADJ_DATA %>% mutate(X = rankZ(X))
    ))
    NULL_MODEL <- suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex,
      REML = FALSE,
      control = lmerControl(optimizer = "Nelder_Mead"),
      data = ADJ_DATA %>% mutate(X = rankZ(X))
    ))
  }
  
  # Fit full model of raw data, without rank-Z transform of trait
  FULL_MODEL_NRZ <- tryCatch(
    suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex + X,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa"),
      data = ADJ_DATA
    )),
    warning = war
  )
  if(is.character(FULL_MODEL_NRZ)){
    cat('Non-rankZ full model throws a warning, refitting with Nelder_Mead optimization\n')
    FULL_MODEL_NRZ <- suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex + X,
      REML = FALSE,
      control = lmerControl(optimizer = "Nelder_Mead"),
      data = ADJ_DATA
    ))
    
  }
  
  
  TEST_RESULTS <- seqPBmodcomp(
    FULL_MODEL,
    NULL_MODEL,
    h = MIN_EXTREME,
    minsim = MIN_SIM,
    maxsim = MAX_SIM,
    chunk.size = CHUNK_SIZE,
    seed = SEED,
    cl = CLUSTER
  )
  
  LS_EFFECTS_13MO <- rbind(
    LS_EFFECTS_13MO,
    data.frame(
      Trait = TRAIT,
      Timepoint = '13 Months',
      PValue_LRT = as.numeric(TEST_RESULTS$test['LRT', 'p.value']),
      PValue_PBS = as.numeric(TEST_RESULTS$test['PBtest', 'p.value']),
      # Total number of simulations
      NSim = as.numeric(TEST_RESULTS$samples[2]),
      # Simulations in which the LRT statistic was greater than that observed
      NExtreme = as.numeric(TEST_RESULTS$n.extreme),
      BootstrapTime = as.numeric(TEST_RESULTS$ctime),
      Effect = as.numeric(fixef(FULL_MODEL_NRZ)['X']),
      SE = as.numeric(sqrt(vcov(FULL_MODEL_NRZ)['X', 'X'])),
      TraitSD = TRAIT_SD
    )
  )
  
  print(VarCorr(FULL_MODEL))
  
  cat('\n -------------------- # \n\n')
}
rm(TRAIT, TRAIT_SD, ADJ_DATA, FULL_MODEL, NULL_MODEL, FULL_MODEL_NRZ, TEST_RESULTS)


sd(
  PhenoData_ClumpDateAdj$Lifespan[PhenoData_ClumpDateAdj$Timepoint == '13 Months'] -
    PhenoData_ClumpDateAdj$AgeAtTest[PhenoData_ClumpDateAdj$Timepoint == '13 Months']
)
# 7.650548

LS_EFFECTS_13MO <- LS_EFFECTS_13MO %>%
  mutate(
    FWER_PBS = p.adjust(PValue_PBS, method = 'holm'),
    # The phenotypes and lifespan were standerdised during model fitting, so to 
    # get the effects and standard errors in units of lifespan (months) per 
    # trait, we need to  multiply them by the standard deviation of the lifespan
    # in months and divide by the standard deviation of the trait
    Effect = Effect * 7.650548 / TraitSD,
    SE = SE * 7.650548 / TraitSD
  ) %>%
  select(
    Trait, PValue_PBS, FWER_PBS, PValue_LRT,
    everything()
  )


#####


################################################################################
## LS effects at 19 months #############################

LS_EFFECTS_19MO <- data.frame(
  Trait = vector(),
  Timepoint = vector(),
  PValue_LRT = vector(),
  PValue_PBS = vector(),
  # Total number of succesfull simulations (LRT statistic > 0)
  NSim = vector(),
  # Simulations in which the LRT statistic was greater than that observed
  NExtreme = vector(),
  BootstrapTime = vector(),
  Effect = vector(),
  SE = vector(),
  TraitSD = vector()
)

for(TRAIT in CBCphenos){
  cat(paste0(
    TRAIT,'\n'
  ))
  
  # Create dataset for model fitting. We will standardise (mean 0, sd = 1) all
  # continuous variables to improve model fitting
  ADJ_DATA <- PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == '19 Months')
  
  if(TRAIT %in% TRANSFORM_TABLE$Phenotype[TRANSFORM_TABLE$Transform == 'log']){
    ADJ_DATA[['X']] <- log(ADJ_DATA[[TRAIT]])
  } else{
    ADJ_DATA[['X']] <- ADJ_DATA[[TRAIT]]
  }
  
  
  TRAIT_SD <- sd(ADJ_DATA[['X']], na.rm = TRUE)
  ADJ_DATA <- ADJ_DATA %>%
    mutate(
      Lifespan = Lifespan - AgeAtTest,
      Lifespan = standardise(Lifespan),
      X = standardise(X)
    ) %>%
    select(
      MouseID, Generation, Sex, X, Lifespan
    )
  
  
  # Fit full model of raw data, but with the trait rank-Z transformed
  FULL_MODEL <- tryCatch(
    suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex + X,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa"),
      data = ADJ_DATA %>% mutate(X = rankZ(X))
    )),
    warning = war
  )
  # Fit null model of raw data, but with the trait rank-Z transformed
  NULL_MODEL <- tryCatch(
    suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa"),
      data = ADJ_DATA %>% mutate(X = rankZ(X))
    )),
    warning = war
  )
  
  if(is.character(FULL_MODEL) | is.character(NULL_MODEL)){
    cat('WARNING: Full and/or null model throws a warning, refitting with Nelder_Mead optimization\n')
    if(is.character(FULL_MODEL)){
      cat('Full model throws a warning\n')
    }
    if(is.character(NULL_MODEL)){
      cat('Null model throws a warning\n')
    }
    # Refit full and null models using Nelder_Mead optimizer
    FULL_MODEL <- suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex + X,
      REML = FALSE,
      control = lmerControl(optimizer = "Nelder_Mead"),
      data = ADJ_DATA %>% mutate(X = rankZ(X))
    ))
    NULL_MODEL <- suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex,
      REML = FALSE,
      control = lmerControl(optimizer = "Nelder_Mead"),
      data = ADJ_DATA %>% mutate(X = rankZ(X))
    ))
  }
  
  # Fit full model of raw data, without rank-Z transform of trait
  FULL_MODEL_NRZ <- tryCatch(
    suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex + X,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa"),
      data = ADJ_DATA
    )),
    warning = war
  )
  if(is.character(FULL_MODEL_NRZ)){
    cat('Non-rankZ full model throws a warning, refitting with Nelder_Mead optimization\n')
    FULL_MODEL_NRZ <- suppressMessages(lmer(
      Lifespan ~ (1|Generation) + Sex + X,
      REML = FALSE,
      control = lmerControl(optimizer = "Nelder_Mead"),
      data = ADJ_DATA
    ))
    
  }
  
  
  TEST_RESULTS <- seqPBmodcomp(
    FULL_MODEL,
    NULL_MODEL,
    h = MIN_EXTREME,
    minsim = MIN_SIM,
    maxsim = MAX_SIM,
    chunk.size = CHUNK_SIZE,
    seed = SEED,
    cl = CLUSTER
  )
  
  LS_EFFECTS_19MO <- rbind(
    LS_EFFECTS_19MO,
    data.frame(
      Trait = TRAIT,
      Timepoint = '19 Months',
      PValue_LRT = as.numeric(TEST_RESULTS$test['LRT', 'p.value']),
      PValue_PBS = as.numeric(TEST_RESULTS$test['PBtest', 'p.value']),
      # Total number of simulations
      NSim = as.numeric(TEST_RESULTS$samples[2]),
      # Simulations in which the LRT statistic was greater than that observed
      NExtreme = as.numeric(TEST_RESULTS$n.extreme),
      BootstrapTime = as.numeric(TEST_RESULTS$ctime),
      Effect = as.numeric(fixef(FULL_MODEL_NRZ)['X']),
      SE = as.numeric(sqrt(vcov(FULL_MODEL_NRZ)['X', 'X'])),
      TraitSD = TRAIT_SD
    )
  )
  
  print(VarCorr(FULL_MODEL))
  
  cat('\n -------------------- # \n\n')
}
rm(TRAIT, TRAIT_SD, ADJ_DATA, FULL_MODEL, NULL_MODEL, FULL_MODEL_NRZ, TEST_RESULTS)


sd(
  PhenoData_ClumpDateAdj$Lifespan[PhenoData_ClumpDateAdj$Timepoint == '19 Months'] -
    PhenoData_ClumpDateAdj$AgeAtTest[PhenoData_ClumpDateAdj$Timepoint == '19 Months']
)
# 6.607949

LS_EFFECTS_19MO <- LS_EFFECTS_19MO %>%
  mutate(
    FWER_PBS = p.adjust(PValue_PBS, method = 'holm'),
    # The phenotypes and lifespan were standerdised during model fitting, so to 
    # get the effects and standard errors in units of lifespan (months) per 
    # trait, we need to  multiply them by the standard deviation of the lifespan
    # in months and divide by the standard deviation of the trait
    Effect = Effect * 6.607949 / TraitSD,
    SE = SE * 6.607949 / TraitSD
  ) %>%
  select(
    Trait, PValue_PBS, FWER_PBS, PValue_LRT,
    everything()
  )


#####


# write_csv(
#   rbind(
#     LS_EFFECTS_07MO,
#     LS_EFFECTS_13MO,
#     LS_EFFECTS_19MO
#   ) %>% 
#     left_join(
#       TRANSFORM_TABLE %>% 
#         mutate(
#           Transform = ifelse(Transform == 'flipped log', 'none', Transform)
#         ) %>% 
#         rename(
#           Trait = Phenotype,
#           TraitTransform = Transform
#         ),
#       by = 'Trait'
#     ),
#   'tables/lifespan_aging_hematology/hematology_effects_on_lifespan.csv'
# )

################################################################################
##  ################################

stopCluster(CLUSTER)

rm(MIN_EXTREME, MIN_SIM, MAX_SIM, CHUNK_SIZE, SEED, CLUSTER)

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
