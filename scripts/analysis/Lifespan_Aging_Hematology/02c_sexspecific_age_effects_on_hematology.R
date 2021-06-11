# 2021-01-28
################################################################################
#
#   This script estimates the sex-specfic effect of age on each of the CBC traits. 
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

#####

MIN_EXTREME <- 2
MIN_SIM <- 10000
MAX_SIM <- 10^7
CHUNK_SIZE <- 10^6
SEED <- 19940418
CLUSTER <- makeForkCluster(50)

################################################################################
## Age-Sex interaction effects ################################

AGESEX_INTERACTION_EFFECTS <- data.frame(
  Trait = vector(),
  PValue_LRT = vector(),
  PValue_PBS = vector(),
  # Total number of succesfull simulations (LRT statistic > 0)
  NSim = vector(),
  # Simulations in which the LRT statistic was greater than that observed 
  NExtreme = vector(),
  BootstrapTime = vector(),
  FemaleEffect = vector(),
  FemaleSE = vector(),
  MaleEffect = vector(),
  MaleSE = vector(),
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
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + clumps + Sex + AgeAtTest + AgeAtTest:Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    # Fit null model of raw data, but with the trait rank-Z transformed
    NULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + clumps + Sex + AgeAtTest,
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
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + clumps + Sex + AgeAtTest + AgeAtTest:Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "Nelder_Mead"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
      NULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + clumps + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "Nelder_Mead"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
    }
    
    # Fit full model of raw data, without rank-Z transform of trait
    FULL_MODEL_NRZ <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + clumps + Sex + AgeAtTest + AgeAtTest:Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA
      )),
      warning = war
    )
    if(is.character(FULL_MODEL_NRZ)){
      cat('WARNING: Non-rankZ full model throws a warning, refitting with Nelder_Mead optimization\n')
      FULL_MODEL_NRZ <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + clumps + Sex + AgeAtTest + AgeAtTest:Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "Nelder_Mead"),
        data = RAW_DATA
      ))
    }
    
  } else{
    
    # Fit full model of raw data, but with the trait rank-Z transformed
    FULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + Sex + AgeAtTest + AgeAtTest:Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    # Fit null model of raw data, but with the trait rank-Z transformed
    NULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + Sex + AgeAtTest,
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
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + Sex + AgeAtTest + AgeAtTest:Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "Nelder_Mead"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
      NULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "Nelder_Mead"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
    }
    
    # Fit full model of raw data, without rank-Z transform of trait
    FULL_MODEL_NRZ <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|DOT) + (1|MouseID) + Sex + AgeAtTest + AgeAtTest:Sex,
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
  
  AGESEX_INTERACTION_EFFECTS <- rbind(
    AGESEX_INTERACTION_EFFECTS,
    data.frame(
      Trait = TRAIT,
      PValue_LRT = as.numeric(TEST_RESULTS$test['LRT', 'p.value']),
      PValue_PBS = as.numeric(TEST_RESULTS$test['PBtest', 'p.value']),
      # Total number of simulations
      NSim = as.numeric(TEST_RESULTS$samples[2]), 
      # Simulations in which the LRT statistic was greater than that observed 
      NExtreme = as.numeric(TEST_RESULTS$n.extreme),
      BootstrapTime = as.numeric(TEST_RESULTS$ctime),
      FemaleEffect = as.numeric(fixef(FULL_MODEL_NRZ)['AgeAtTest']),
      FemaleSE = as.numeric(sqrt(vcov(FULL_MODEL_NRZ)['AgeAtTest', 'AgeAtTest'])),
      MaleEffect = as.numeric(
        fixef(FULL_MODEL_NRZ)['SexMale:AgeAtTest'] +
          fixef(FULL_MODEL_NRZ)['AgeAtTest']
      ),
      MaleSE = as.numeric(sqrt(
        vcov(FULL_MODEL_NRZ)['AgeAtTest', 'AgeAtTest'] +
          vcov(FULL_MODEL_NRZ)['SexMale:AgeAtTest', 'SexMale:AgeAtTest'] +
          2*vcov(FULL_MODEL_NRZ)['SexMale:AgeAtTest', 'AgeAtTest']
      )),
      TraitSD = TRAIT_SD
    )
  )

  write_csv(
    AGESEX_INTERACTION_EFFECTS,
    'out/agesex_interaction_effects_on_hematology_traits_e7_bigchunks.csv'
  )
  
  print(VarCorr(FULL_MODEL))
  
  cat('\n -------------------- # \n\n')
}
rm(TRAIT, TRAIT_SD, RAW_DATA, FULL_MODEL, NULL_MODEL, FULL_MODEL_NRZ, TEST_RESULTS)
#   rbc
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.75075 
# Generation:DOT (Intercept) 0.22860 
# Generation     (Intercept) 0.00000 
# Residual                   0.61315 
# 
# -------------------- # 
#   
#   ch
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.81288 
# Generation:DOT (Intercept) 0.21413 
# Generation     (Intercept) 0.17980 
# Residual                   0.45414 
# 
# -------------------- # 
#   
#   chcm
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.52512 
# Generation:DOT (Intercept) 0.62472 
# Generation     (Intercept) 0.29487 
# Residual                   0.38845 
# 
# -------------------- # 
#   
#   hdw
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.63875 
# Generation:DOT (Intercept) 0.42747 
# Generation     (Intercept) 0.17117 
# Residual                   0.58130 
# 
# -------------------- # 
#   
#   mcv
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.60847 
# Generation:DOT (Intercept) 0.44583 
# Generation     (Intercept) 0.21226 
# Residual                   0.46095 
# 
# -------------------- # 
#   
#   rdw
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.72704 
# Generation:DOT (Intercept) 0.18751 
# Generation     (Intercept) 0.10242 
# Residual                   0.62780 
# 
# -------------------- # 
#   
#   hct
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.51788 
# Generation:DOT (Intercept) 0.43878 
# Generation     (Intercept) 0.23873 
# Residual                   0.52639 
# 
# -------------------- # 
#   
#   hgb
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.70210 
# Generation:DOT (Intercept) 0.17945 
# Generation     (Intercept) 0.13018 
# Residual                   0.62767 
# 
# -------------------- # 
#   
#   wbc
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.74361 
# Generation:DOT (Intercept) 0.24426 
# Generation     (Intercept) 0.00000 
# Residual                   0.61354 
# 
# -------------------- # 
#   
#   nlr
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.67948 
# Generation:DOT (Intercept) 0.11698 
# Generation     (Intercept) 0.00000 
# Residual                   0.60382 
# 
# -------------------- # 
#   
#   n.lymph
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.73556 
# Generation:DOT (Intercept) 0.23211 
# Generation     (Intercept) 0.00000 
# Residual                   0.62575 
# 
# -------------------- # 
#   
#   n.neut
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.644460
# Generation:DOT (Intercept) 0.193868
# Generation     (Intercept) 0.052633
# Residual                   0.624761
# 
# -------------------- # 
#   
#   n.mono
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.62391 
# Generation:DOT (Intercept) 0.27479 
# Generation     (Intercept) 0.00000 
# Residual                   0.64125 
# 
# -------------------- # 
#   
#   n.eos
# Groups         Name        Std.Dev.  
# MouseID        (Intercept) 6.0943e-01
# Generation:DOT (Intercept) 1.5978e-01
# Generation     (Intercept) 2.0475e-08
# Residual                   6.1723e-01
# 
# -------------------- # 
#   
#   plt
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.63582 
# Generation:DOT (Intercept) 0.12255 
# Generation     (Intercept) 0.07016 
# Residual                   0.57965 
# 
# -------------------- # 
#   
#   mpv
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.30493 
# Generation:DOT (Intercept) 0.64238 
# Generation     (Intercept) 0.41517 
# Residual                   0.38130 
# 
# -------------------- # 
#   
#   mpm
# Groups         Name        Std.Dev.
# MouseID        (Intercept) 0.70165 
# Generation:DOT (Intercept) 0.42038 
# Generation     (Intercept) 0.29594 
# Residual                   0.45944 
# 
# -------------------- # 

sd(PhenoData_Raw$AgeAtTest)
# 4.852629

AGESEX_INTERACTION_EFFECTS <- AGESEX_INTERACTION_EFFECTS %>% 
  mutate(
    FWER_PBS = p.adjust(PValue_PBS, method = 'holm'),
    # Age and the phenotypes were standerdised during model fitting, so to get 
    # the age effects and standard errors in units of trait/month, we need to 
    # divide them by the standard deviation of age in months (4.852629 months) 
    # and multiply them by the standard deviation of the trait
    FemaleEffect = FemaleEffect / 4.852629 * TraitSD,
    FemaleSE = FemaleSE / 4.852629 * TraitSD,
    MaleEffect = MaleEffect / 4.852629 * TraitSD,
    MaleSE = MaleSE / 4.852629 * TraitSD,
  ) %>% 
  select(
    Trait, PValue_PBS, FWER_PBS, PValue_LRT,
    everything()
  )

# write_csv(
#   AGESEX_INTERACTION_EFFECTS,
#   'tables/lifespan_aging_hematology/agesex_interaction_effects_on_hematology_traits.csv'
# )

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
