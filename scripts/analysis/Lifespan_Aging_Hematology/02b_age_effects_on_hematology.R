# 2021-01-28
################################################################################
#
#   This script estimates the effect of age on each of the CBC traits. 
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
  nchunk <- maxsim %/% chunk.size
  LRTstat <- pbkrtest:::getLRT.merMod(largeModel, smallModel)
  ref <- NULL
  for (ii in 1:nchunk) {
    chunk_seed <- sample.int(99999:1000000, size = 1)
    ref <- c(ref, PBrefdist(largeModel, smallModel, nsim = chunk.size, seed = chunk_seed, ...))
    n.extreme <- sum(ref > LRTstat["tobs"])
    if (n.extreme >= h & length(ref > 0) >= minsim)
      break
  }
  ans <- PBmodcomp(largeModel, smallModel, ref = ref)
  ans$ctime <- (proc.time() - t.start)[3]
  ans
}

#####


################################################################################
## load data etc #########################################################

load('data/processed/phenotypic/CBC_AnalysisEnvironment.RData')

#####

MIN_EXTREME <- 1
MIN_SIM <- 10000
MAX_SIM <- 10^17
CHUNK_SIZE <- MIN_SIM + 10
SEED <- 19940418
CLUSTER <- makeCluster(rep("localhost", 32))

################################################################################
## Age effects ################################

AGE_EFFECTS <- data.frame(
  Trait = vector(),
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
  
  TRAIT_SD <- sd(PhenoData_Raw[[TRAIT]], na.rm = TRUE)
  
  # Create dataset for model fitting. We will standardise (mean 0, sd = 1) all 
  # continuous variables to improve model fitting
  RAW_DATA <- PhenoData_Raw
  RAW_DATA[['Y']] <- RAW_DATA[[TRAIT]]
  RAW_DATA <- RAW_DATA %>% 
    mutate(
      DOT = as.character(DOT),
      AgeAtTest = standardise(AgeAtTest),
      clumps = standardise(log(clumps)),
      Y = standardise(Y)
    ) %>% 
    select(
      MouseID, Generation, DOT, clumps, AgeAtTest, Timepoint, Sex, Y
    )
  
  
  if(TRAIT %in% Clump_Affected_Traits){
    
    # Fit full model of raw data, but with the trait rank-Z transformed
    FULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + clumps + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    # Fit null model of raw data, but with the trait rank-Z transformed
    NULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + clumps + Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
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
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + clumps + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "Nelder_Mead"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
      NULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + clumps + Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "Nelder_Mead"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
    }
    
    # Fit full model of raw data, without rank-Z transform of trait
    FULL_MODEL_NRZ <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + clumps + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA
      )),
      warning = war
    )
    if(is.character(FULL_MODEL_NRZ)){
      cat('WARNING: Non-rankZ full model throws a warning, refitting with Nelder_Mead optimization\n')
      FULL_MODEL_NRZ <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + clumps + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "Nelder_Mead"),
        data = RAW_DATA
      ))
    }
    
  } else{
    
    # Fit full model of raw data, but with the trait rank-Z transformed
    FULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    # Fit null model of raw data, but with the trait rank-Z transformed
    NULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
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
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "Nelder_Mead"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
      NULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "Nelder_Mead"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
    }
    
    # Fit full model of raw data, without rank-Z transform of trait
    FULL_MODEL_NRZ <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA
      )),
      warning = war
    )
    if(is.character(FULL_MODEL_NRZ)){
      cat('Non-rankZ full model throws a warning, refitting with Nelder_Mead optimization\n')
      FULL_MODEL_NRZ <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "Nelder_Mead"),
        data = RAW_DATA
      ))
      
    }
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
  
  AGE_EFFECTS <- rbind(
    AGE_EFFECTS,
    data.frame(
      Trait = TRAIT,
      PValue_LRT = as.numeric(TEST_RESULTS$test['LRT', 'p.value']),
      PValue_PBS = as.numeric(TEST_RESULTS$test['PBtest', 'p.value']),
      # Total number of simulations
      NSim = as.numeric(TEST_RESULTS$samples[2]), 
      # Simulations in which the LRT statistic was greater than that observed 
      NExtreme = as.numeric(TEST_RESULTS$n.extreme),
      BootstrapTime = as.numeric(TEST_RESULTS$ctime),
      Effect = as.numeric(fixef(FULL_MODEL_NRZ)['AgeAtTest']),
      SE = as.numeric(sqrt(vcov(FULL_MODEL_NRZ)['AgeAtTest', 'AgeAtTest'])),
      TraitSD = TRAIT_SD
    )
  )
  
  print(VarCorr(FULL_MODEL))
  
  cat('\n -------------------- # \n\n')
}
rm(TRAIT, TRAIT_SD, RAW_DATA, FULL_MODEL, NULL_MODEL, FULL_MODEL_NRZ, TEST_RESULTS)
#   rbc
# Groups         Name        Std.Dev.  
# MouseID        (Intercept) 7.4952e-01
# Generation:DOT (Intercept) 2.2942e-01
# Generation     (Intercept) 1.0810e-07
# Residual                   6.1482e-01
# 
# -------------------- # 
#   
#   ch
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.80912 
# Generation:DOT (Intercept) 0.21665 
# Generation     (Intercept) 0.18279 
# Residual                   0.46507 
# 
# -------------------- # 
#   
#   chcm
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.52508 
# Generation:DOT (Intercept) 0.62464 
# Generation     (Intercept) 0.29440 
# Residual                   0.38865 
# 
# -------------------- # 
#   
#   hdw
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.63880 
# Generation:DOT (Intercept) 0.42694 
# Generation     (Intercept) 0.17111 
# Residual                   0.58153 
# 
# -------------------- # 
#   
#   mcv
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.60344 
# Generation:DOT (Intercept) 0.44531 
# Generation     (Intercept) 0.21030 
# Residual                   0.47486 
# 
# -------------------- # 
#   
#   rdw
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.72652 
# Generation:DOT (Intercept) 0.18938 
# Generation     (Intercept) 0.10171 
# Residual                   0.62868 
# 
# -------------------- # 
#   
#   hct
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.51764 
# Generation:DOT (Intercept) 0.43802 
# Generation     (Intercept) 0.23759 
# Residual                   0.52897 
# 
# -------------------- # 
#   
#   hgb
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.70238 
# Generation:DOT (Intercept) 0.17940 
# Generation     (Intercept) 0.12967 
# Residual                   0.62893 
# 
# -------------------- # 
#   
#   wbc
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.74370 
# Generation:DOT (Intercept) 0.24391 
# Generation     (Intercept) 0.00000 
# Residual                   0.61379 
# 
# -------------------- # 
#   
#   nlr
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.67807 
# Generation:DOT (Intercept) 0.11793 
# Generation     (Intercept) 0.00000 
# Residual                   0.60625 
# 
# -------------------- # 
#   
#   n.lymph
# Groups         Name        Std.Dev.  
# MouseID        (Intercept) 7.3568e-01
# Generation:DOT (Intercept) 2.3185e-01
# Generation     (Intercept) 5.4943e-08
# Residual                   6.2582e-01
# 
# -------------------- # 
#   
#   n.neut
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.643061
# Generation:DOT (Intercept) 0.193570
# Generation     (Intercept) 0.054061
# Residual                   0.628099
# 
# -------------------- # 
#   
#   n.mono
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.62359 
# Generation:DOT (Intercept) 0.27434 
# Generation     (Intercept) 0.00000 
# Residual                   0.64165 
# 
# -------------------- # 
#   
#   n.eos
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.60949 
# Generation:DOT (Intercept) 0.16025 
# Generation     (Intercept) 0.00000 
# Residual                   0.61728 
# 
# -------------------- # 
#   
#   plt
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.631218
# Generation:DOT (Intercept) 0.118605
# Generation     (Intercept) 0.071786
# Residual                   0.592795
# 
# -------------------- # 
#   
#   mpv
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.30503 
# Generation:DOT (Intercept) 0.64213 
# Generation     (Intercept) 0.41537 
# Residual                   0.38133 
# 
# -------------------- # 
#   
#   mpm
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.69923 
# Generation:DOT (Intercept) 0.42099 
# Generation     (Intercept) 0.29660 
# Residual                   0.46373 
# 
# -------------------- # 

sd(PhenoData_Raw$AgeAtTest)
# 4.852629

AGE_EFFECTS <- AGE_EFFECTS %>% 
  mutate(
    FWER_PBS = p.adjust(PValue_PBS, method = 'holm'),
    # Age and the phenotypes were standerdised during model fitting, so to get 
    # the age effects and standard errors in units of trait/month, we need to 
    # divide them by the standard deviation of age in months (4.852629 months) 
    # and multiply them by the standard deviation of the trait
    Effect = Effect / 4.852629 * TraitSD,
    SE = SE / 4.852629 * TraitSD
  ) %>% 
  select(
    Trait, PValue_PBS, FWER_PBS, PValue_LRT,
    everything()
  )

# write_csv(
#   AGE_EFFECTS,
#   'tables/lifespan_aging_hematology/age_effects_on_hematology_traits.csv'
# )

write_csv(
  AGE_EFFECTS,
  'out/age_effects_on_hematology_traits.csv'
)

#####


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