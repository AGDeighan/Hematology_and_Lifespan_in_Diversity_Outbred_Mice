# 2021-01-28
################################################################################
#
#   This script estimates the marginal effects of each CBC component on survival
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
# library(RLRsim)

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

seqPBmodcomp <- function(largeModel, smallModel, h = 20, nsim = 1000, seed, ...) {
  # The seqPBmodcomp function in the pbkrtest package is not actually defined,
  # so we copied the following function from the original paper, "A Kenward-Roger 
  # Approximation and Parametric Bootstrap Methods for Tests in Linear Mixed 
  # Models The R Package pbkrtest"  
  set.seed(seed = seed)
  t.start <- proc.time()
  chunk.size <- 50
  nchunk <- nsim %/% chunk.size
  LRTstat <- pbkrtest:::getLRT.merMod(largeModel, smallModel)
  ref <- NULL
  for (ii in 1:nchunk) {
    chunk_seed <- sample.int(99999:1000000, size = 1)
    ref <- c(ref, PBrefdist(largeModel, smallModel, nsim = chunk.size, seed = chunk_seed, ...))
    n.extreme <- sum(ref > LRTstat["tobs"])
    if (n.extreme >= h)
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

NUMBER_OF_SIM_TO_PERFORM <- 10000

################################################################################
## Sex effects ################################

SEX_EFFECTS <- data.frame(
  Trait = vector(),
  PValue_LRT = vector(),
  PValut_PBS = vector(),
  # Total number of simulations
  NSim = vector(),
  # Simulations in which the LRT statistic was greater than that observed 
  NExtreme = vector(),
  BootstrapTime = vector(),
  Effect = vector(),
  SE = vector(),
  TraitSD = vector()
)

# Sex as fixed effect
CLUSTER <- makeCluster(rep("localhost", 6))
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
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + AgeAtTest + Sex,
        REML = FALSE,
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    # Fit null model of raw data, but with the trait rank-Z transformed
    NULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + AgeAtTest,
        REML = FALSE,
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    
    if(is.character(FULL_MODEL) | is.character(NULL_MODEL)){
      cat('WARNING: Full and/or null model throws a warning, refitting with BOBYQA optimization\n')
      if(is.character(FULL_MODEL)){
        cat('Full model throws a warning\n')
      }
      if(is.character(NULL_MODEL)){
        cat('Null model throws a warning\n')
      }
      # Refit full and null models using BOBYQA optimizer
      FULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + AgeAtTest + Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
      NULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
    }
    
    
    
    # Fit full model of raw data, without rank-Z transform of trait
    FULL_MODEL_NRZ <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + AgeAtTest + Sex,
        REML = FALSE,
        data = RAW_DATA
      )),
      warning = war
    )
    if(is.character(FULL_MODEL_NRZ)){
      cat('WARNING: Non-rankZ full model throws a warning, refitting with BOBYQA optimization\n')
      FULL_MODEL_NRZ <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + AgeAtTest + Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA
      ))
    }
    
  } else{
    
    # Fit full model of raw data, but with the trait rank-Z transformed
    FULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + AgeAtTest + Sex,
        REML = FALSE,
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    # Fit null model of raw data, but with the trait rank-Z transformed
    NULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + AgeAtTest,
        REML = FALSE,
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    
    if(is.character(FULL_MODEL) | is.character(NULL_MODEL)){
      cat('WARNING: Full and/or null model throws a warning, refitting with BOBYQA optimization\n')
      if(is.character(FULL_MODEL)){
        cat('Full model throws a warning\n')
      }
      if(is.character(NULL_MODEL)){
        cat('Null model throws a warning\n')
      }
      # Refit full and null models using BOBYQA optimizer
      FULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + AgeAtTest + Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
      NULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
    }
    
    # Fit full model of raw data, without rank-Z transform of trait
    FULL_MODEL_NRZ <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + AgeAtTest + Sex,
        REML = FALSE,
        data = RAW_DATA
      )),
      warning = war
    )
    if(is.character(FULL_MODEL_NRZ)){
      cat('Non-rankZ full model throws a warning, refitting with BOBYQA optimization\n')
      FULL_MODEL_NRZ <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + AgeAtTest + Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA
      ))
      
    }
  }
  
  # Perform 1000 parametric bootstraps or enough until 20 simulated LRT 
  # statistics greater than the observed statistic are generated. This is
  # done in batches of 50, so it is possible to actually generate more than
  # 20 extreme LRT statistics
  TEST_RESULTS <- seqPBmodcomp(
    FULL_MODEL,
    NULL_MODEL,
    h = 20,
    nsim = NUMBER_OF_SIM_TO_PERFORM,
    seed = 19940418,
    cl = CLUSTER
  )
  
  SEX_EFFECTS <- rbind(
    SEX_EFFECTS,
    data.frame(
      Trait = TRAIT,
      PValue_LRT = as.numeric(TEST_RESULTS$test['LRT', 'p.value']),
      PValut_PBS = as.numeric(TEST_RESULTS$test['PBtest', 'p.value']),
      # Total number of simulations
      NSim = as.numeric(TEST_RESULTS$samples[2]), 
      # Simulations in which the LRT statistic was greater than that observed 
      NExtreme = as.numeric(TEST_RESULTS$n.extreme),
      BootstrapTime = as.numeric(TEST_RESULTS$ctime),
      Effect = as.numeric(fixef(FULL_MODEL_NRZ)['SexMale']),
      SE = as.numeric(sqrt(vcov(FULL_MODEL_NRZ)['SexMale', 'SexMale'])),
      TraitSD = TRAIT_SD
    )
  )
  
  print(VarCorr(FULL_MODEL))
  
  cat('\n -------------------- # \n\n')
}
rm(TRAIT, TRAIT_SD, RAW_DATA, FULL_MODEL, NULL_MODEL, FULL_MODEL_NRZ, TEST_RESULTS)
#   rbc
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.74952 
# DOT        (Intercept) 0.22942 
# Generation (Intercept) 0.00000 
# Residual               0.61482 
# 
# -------------------- # 
#   
#   ch
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.80912 
# DOT        (Intercept) 0.21666 
# Generation (Intercept) 0.18279 
# Residual               0.46507 
# 
# -------------------- # 
#   
#   chcm
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.52508 
# DOT        (Intercept) 0.62464 
# Generation (Intercept) 0.29440 
# Residual               0.38865 
# 
# -------------------- # 
#   
#   hdw
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.63880 
# DOT        (Intercept) 0.42694 
# Generation (Intercept) 0.17109 
# Residual               0.58153 
# 
# -------------------- # 
#   
#   mcv
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.60344 
# DOT        (Intercept) 0.44530 
# Generation (Intercept) 0.21030 
# Residual               0.47486 
# 
# -------------------- # 
#   
#   rdw
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.72652 
# DOT        (Intercept) 0.18938 
# Generation (Intercept) 0.10172 
# Residual               0.62868 
# 
# -------------------- # 
#   
#   hct
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.51764 
# DOT        (Intercept) 0.43803 
# Generation (Intercept) 0.23758 
# Residual               0.52896 
# 
# -------------------- # 
#   
#   hgb
# Non-rankZ full model throws a warning, refitting with BOBYQA optimization
#   In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#     Model failed to converge with max|grad| = 0.0020389 (tol = 0.002, component 1)
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.70238 
# DOT        (Intercept) 0.17940 
# Generation (Intercept) 0.12968 
# Residual               0.62893 
# 
# -------------------- # 
#   
#   wbc
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.74369 
# DOT        (Intercept) 0.24391 
# Generation (Intercept) 0.00000 
# Residual               0.61379 
# 
# -------------------- # 
#   
#   nlr
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.67807 
# DOT        (Intercept) 0.11794 
# Generation (Intercept) 0.00000 
# Residual               0.60625 
# 
# -------------------- # 
#   
#   n.lymph
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.73568 
# DOT        (Intercept) 0.23185 
# Generation (Intercept) 0.00000 
# Residual               0.62582 
# 
# -------------------- # 
#   
#   n.neut
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.64306 
# DOT        (Intercept) 0.19357 
# Generation (Intercept) 0.05406 
# Residual               0.62810 
# 
# -------------------- # 
#   
#   n.mono
# WARNING: Full and/or null model throws a warning, refitting with BOBYQA optimization
#   Full model throws a warning
#     In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#       Model failed to converge with max|grad| = 0.00548718 (tol = 0.002, component 1)
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.62359 
# DOT        (Intercept) 0.27435 
# Generation (Intercept) 0.00000 
# Residual               0.64165 
# 
# -------------------- # 
#   
#   n.eos
# Groups     Name        Std.Dev.  
# MouseID    (Intercept) 0.60949371
# DOT        (Intercept) 0.16025528
# Generation (Intercept) 0.00003816
# Residual               0.61728397
# 
# -------------------- # 
#   
#   plt
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.63121 
# DOT        (Intercept) 0.11861 
# Generation (Intercept) 0.07179 
# Residual               0.59280 
# 
# -------------------- # 
#   
#   mpv
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.30503 
# DOT        (Intercept) 0.64213 
# Generation (Intercept) 0.41538 
# Residual               0.38133 
# 
# -------------------- # 
#   
#   mpm
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.69923 
# DOT        (Intercept) 0.42099 
# Generation (Intercept) 0.29660 
# Residual               0.46373 
# 
# -------------------- # 




# ------------------------- Checking hgb -------------------------------------- #

# Check that BOBYQA fitting of hgb doesn't also throw convergence warnings
TRAIT <- 'hgb'

RAW_DATA <- PhenoData_Raw
RAW_DATA[['Y']] <- RAW_DATA[[TRAIT]]
RAW_DATA <- RAW_DATA %>% 
  mutate(
    AgeAtTest = standardise(AgeAtTest),
    clumps = standardise(log(clumps)),
    Y = standardise(Y)
  ) %>% 
  select(
    MouseID, Generation, DOT, clumps, AgeAtTest, Timepoint, Sex, Y
  )

FULL_MODEL_NRZ <- lmer(
  Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + AgeAtTest + Sex,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa"),
  data = RAW_DATA
)

VarCorr(FULL_MODEL_NRZ)
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.605352
# DOT        (Intercept) 0.156035
# Generation (Intercept) 0.074616
# Residual               0.741361

rm(TRAIT, RAW_DATA, FULL_MODEL_NRZ)

# ---------------------------------------------------------------------------- #


# --------------------------- Checking n.mono -------------------------------- #

# Check that BOBYQA fittings of n.mono don't also throw convergence 
# warnings,  they don't
TRAIT <- 'n.mono'

RAW_DATA <- PhenoData_Raw
RAW_DATA[['Y']] <- RAW_DATA[[TRAIT]]
RAW_DATA <- RAW_DATA %>% 
  mutate(
    AgeAtTest = standardise(AgeAtTest),
    clumps = standardise(log(clumps)),
    Y = standardise(Y)
  ) %>% 
  select(
    MouseID, Generation, DOT, clumps, AgeAtTest, Timepoint, Sex, Y
  )

FULL_MODEL <- lmer(
  Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + AgeAtTest + Sex,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa"),
  data = RAW_DATA %>% mutate(Y = rankZ(Y))
)
NULL_MODEL <- lmer(
  Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + AgeAtTest + Sex,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa"),
  data = RAW_DATA %>% mutate(Y = rankZ(Y))
)
VarCorr(FULL_MODEL)
# Groups     Name        Std.Dev.
#  MouseID    (Intercept) 0.62359 
#  DOT        (Intercept) 0.27434 
#  Generation (Intercept) 0.00000 
#  Residual               0.64165 

rm(TRAIT, RAW_DATA, FULL_MODEL, NULL_MODEL)

# ---------------------------------------------------------------------------- #


stopCluster(CLUSTER); rm(CLUSTER)

write_csv(
  SEX_EFFECTS,
  'tables/lifespan_aging_hematology/sex_effects_on_hematology_traits.csv'
)

#####


################################################################################
## Age effects ################################

AGE_EFFECTS <- data.frame(
  Trait = vector(),
  PValue_LRT = vector(),
  PValut_PBS = vector(),
  # Total number of simulations
  NSim = vector(),
  # Simulations in which the LRT statistic was greater than that observed 
  NExtreme = vector(),
  BootstrapTime = vector(),
  Effect = vector(),
  SE = vector(),
  TraitSD = vector()
)

# Sex as fixed effect
CLUSTER <- makeCluster(rep("localhost", 6))
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
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + Sex + AgeAtTest,
        REML = FALSE,
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    # Fit null model of raw data, but with the trait rank-Z transformed
    NULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + Sex,
        REML = FALSE,
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    if(is.character(FULL_MODEL) | is.character(NULL_MODEL)){
      cat('WARNING: Full and/or null model throws a warning, refitting with BOBYQA optimization\n')
      if(is.character(FULL_MODEL)){
        cat('Full model throws a warning\n')
      }
      if(is.character(NULL_MODEL)){
        cat('Null model throws a warning\n')
      }
      # Refit full and null models using BOBYQA optimizer
      FULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
      NULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
    }
    
    # Fit full model of raw data, without rank-Z transform of trait
    FULL_MODEL_NRZ <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + Sex + AgeAtTest,
        REML = FALSE,
        data = RAW_DATA
      )),
      warning = war
    )
    if(is.character(FULL_MODEL_NRZ)){
      cat('WARNING: Non-rankZ full model throws a warning, refitting with BOBYQA optimization\n')
      FULL_MODEL_NRZ <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA
      ))
    }
    
  } else{
    
    # Fit full model of raw data, but with the trait rank-Z transformed
    FULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest,
        REML = FALSE,
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    # Fit null model of raw data, but with the trait rank-Z transformed
    NULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex,
        REML = FALSE,
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    if(is.character(FULL_MODEL) | is.character(NULL_MODEL)){
      cat('WARNING: Full and/or null model throws a warning, refitting with BOBYQA optimization\n')
      if(is.character(FULL_MODEL)){
        cat('Full model throws a warning\n')
      }
      if(is.character(NULL_MODEL)){
        cat('Null model throws a warning\n')
      }
      # Refit full and null models using BOBYQA optimizer
      FULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
      NULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
    }
    
    # Fit full model of raw data, without rank-Z transform of trait
    FULL_MODEL_NRZ <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest,
        REML = FALSE,
        data = RAW_DATA
      )),
      warning = war
    )
    if(is.character(FULL_MODEL_NRZ)){
      cat('Non-rankZ full model throws a warning, refitting with BOBYQA optimization\n')
      FULL_MODEL_NRZ <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA
      ))
      
    }
  }
  
  # Perform 1000 parametric bootstraps or enough until 20 simulated LRT 
  # statistics greater than the observed statistic are generated. This is
  # done in batches of 50, so it is possible to actually generate more than
  # 20 extreme LRT statistics
  TEST_RESULTS <- seqPBmodcomp(
    FULL_MODEL,
    NULL_MODEL,
    h = 20,
    nsim = NUMBER_OF_SIM_TO_PERFORM,
    seed = 19940418,
    cl = CLUSTER
  )
  
  AGE_EFFECTS <- rbind(
    AGE_EFFECTS,
    data.frame(
      Trait = TRAIT,
      PValue_LRT = as.numeric(TEST_RESULTS$test['LRT', 'p.value']),
      PValut_PBS = as.numeric(TEST_RESULTS$test['PBtest', 'p.value']),
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
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.74952 
# DOT        (Intercept) 0.22942 
# Generation (Intercept) 0.00000 
# Residual               0.61482 
# 
# -------------------- # 
#   
#   ch
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.80912 
# DOT        (Intercept) 0.21666 
# Generation (Intercept) 0.18279 
# Residual               0.46507 
# 
# -------------------- # 
#   
#   chcm
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.52508 
# DOT        (Intercept) 0.62464 
# Generation (Intercept) 0.29440 
# Residual               0.38865 
# 
# -------------------- # 
#   
#   hdw
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.63880 
# DOT        (Intercept) 0.42694 
# Generation (Intercept) 0.17109 
# Residual               0.58153 
# 
# -------------------- # 
#   
#   mcv
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.60344 
# DOT        (Intercept) 0.44530 
# Generation (Intercept) 0.21030 
# Residual               0.47486 
# 
# -------------------- # 
#   
#   rdw
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.72652 
# DOT        (Intercept) 0.18938 
# Generation (Intercept) 0.10172 
# Residual               0.62868 
# 
# -------------------- # 
#   
#   hct
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.51764 
# DOT        (Intercept) 0.43803 
# Generation (Intercept) 0.23758 
# Residual               0.52896 
# 
# -------------------- # 
#   
#   hgb
# WARNING: Full and/or null model throws a warning, refitting with BOBYQA optimization
#   Null model throws a warning
#     In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#       Model failed to converge with max|grad| = 0.00832725 (tol = 0.002, component 1)
# 
# Non-rankZ full model throws a warning, refitting with BOBYQA optimization
#   In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#     Model failed to converge with max|grad| = 0.00203891 (tol = 0.002, component 1)
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.70238 
# DOT        (Intercept) 0.17940 
# Generation (Intercept) 0.12968 
# Residual               0.62893 
# 
# -------------------- # 
#   
#   wbc
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.74369 
# DOT        (Intercept) 0.24391 
# Generation (Intercept) 0.00000 
# Residual               0.61379 
# 
# -------------------- # 
#   
#   nlr
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.67807 
# DOT        (Intercept) 0.11794 
# Generation (Intercept) 0.00000 
# Residual               0.60625 
# 
# -------------------- # 
#   
#   n.lymph
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.73568 
# DOT        (Intercept) 0.23185 
# Generation (Intercept) 0.00000 
# Residual               0.62582 
# 
# -------------------- # 
#   
#   n.neut
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.64306 
# DOT        (Intercept) 0.19357 
# Generation (Intercept) 0.05406 
# Residual               0.62810 
# 
# -------------------- # 
#   
#   n.mono
# WARNING: Full and/or null model throws a warning, refitting with BOBYQA optimization
#   Full model throws a warning
#     In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#       Model failed to converge with max|grad| = 0.00548717 (tol = 0.002, component 1)
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.62359 
# DOT        (Intercept) 0.27435 
# Generation (Intercept) 0.00000 
# Residual               0.64165 
# 
# -------------------- # 
#   
#   n.eos
# Groups     Name        Std.Dev.  
# MouseID    (Intercept) 0.60949371
# DOT        (Intercept) 0.16025528
# Generation (Intercept) 0.00003816
# Residual               0.61728397
# 
# -------------------- # 
#   
#   plt
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.631218
# DOT        (Intercept) 0.118606
# Generation (Intercept) 0.071791
# Residual               0.592795
# 
# -------------------- # 
#   
#   mpv
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.30503 
# DOT        (Intercept) 0.64213 
# Generation (Intercept) 0.41538 
# Residual               0.38133 
# 
# -------------------- # 
#   
#   mpm
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.69923 
# DOT        (Intercept) 0.42099 
# Generation (Intercept) 0.29660 
# Residual               0.46373 
# 
# -------------------- # 



# ------------------------- Checking hgb -------------------------------------- #

# Check that BOBYQA fittings of hgb don't also throw convergence warnings,
# they do not.

TRAIT <- 'hgb'

RAW_DATA <- PhenoData_Raw
RAW_DATA[['Y']] <- RAW_DATA[[TRAIT]]
RAW_DATA <- RAW_DATA %>% 
  mutate(
    AgeAtTest = standardise(AgeAtTest),
    clumps = standardise(log(clumps)),
    Y = standardise(Y)
  ) %>% 
  select(
    MouseID, Generation, DOT, clumps, AgeAtTest, Timepoint, Sex, Y
  )

FULL_MODEL <- lmer(
  Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa"),
  data = RAW_DATA %>% mutate(Y = rankZ(Y))
)
NULL_MODEL <- lmer(
  Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa"),
  data = RAW_DATA %>% mutate(Y = rankZ(Y))
)
VarCorr(FULL_MODEL)
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.70238 
# DOT        (Intercept) 0.17940 
# Generation (Intercept) 0.12967 
# Residual               0.62893 

FULL_MODEL_NRZ <- suppressMessages(lmer(
  Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa"),
  data = RAW_DATA
))

rm(TRAIT, RAW_DATA, FULL_MODEL, NULL_MODEL, FULL_MODEL_NRZ)

# ---------------------------------------------------------------------------- #


# --------------------------- Checking n.mono -------------------------------- #

# Check that BOBYQA fittings of n.mono don't also throw convergence 
# warnings, they do not.

TRAIT <- 'n.mono'

RAW_DATA <- PhenoData_Raw
RAW_DATA[['Y']] <- RAW_DATA[[TRAIT]]
RAW_DATA <- RAW_DATA %>% 
  mutate(
    AgeAtTest = standardise(AgeAtTest),
    clumps = standardise(log(clumps)),
    Y = standardise(Y)
  ) %>% 
  select(
    MouseID, Generation, DOT, clumps, AgeAtTest, Timepoint, Sex, Y
  )

FULL_MODEL <- lmer(
  Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa"),
  data = RAW_DATA %>% mutate(Y = rankZ(Y))
)
NULL_MODEL <- lmer(
  Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa"),
  data = RAW_DATA %>% mutate(Y = rankZ(Y))
)
VarCorr(FULL_MODEL)
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.62359 
# DOT        (Intercept) 0.27434 
# Generation (Intercept) 0.00000 
# Residual               0.64165 

rm(TRAIT, RAW_DATA, FULL_MODEL, NULL_MODEL)

# ---------------------------------------------------------------------------- #


stopCluster(CLUSTER); rm(CLUSTER)

write_csv(
  AGE_EFFECTS,
  'tables/lifespan_aging_hematology/age_effects_on_hematology_traits.csv'
)

#####


################################################################################
## Age-Sex interaction effects ################################

AGESEX_INTERACTION_EFFECTS <- data.frame(
  Trait = vector(),
  PValue_LRT = vector(),
  PValut_PBS = vector(),
  # Total number of simulations
  NSim = vector(),
  # Simulations in which the LRT statistic was greater than that observed 
  NExtreme = vector(),
  BootstrapTime = vector(),
  Effect = vector(),
  SE = vector(),
  TraitSD = vector()
)

# Sex as fixed effect
CLUSTER <- makeCluster(rep("localhost", 6))
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
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + Sex + AgeAtTest + AgeAtTest:Sex,
        REML = FALSE,
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    # Fit null model of raw data, but with the trait rank-Z transformed
    NULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + Sex + AgeAtTest,
        REML = FALSE,
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    if(is.character(FULL_MODEL) | is.character(NULL_MODEL)){
      cat('WARNING: Full and/or null model throws a warning, refitting with BOBYQA optimization\n')
      if(is.character(FULL_MODEL)){
        cat('Full model throws a warning\n')
      }
      if(is.character(NULL_MODEL)){
        cat('Null model throws a warning\n')
      }
      # Refit full and null models using BOBYQA optimizer
      FULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + Sex + AgeAtTest + AgeAtTest:Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
      NULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
    }
    
    # Fit full model of raw data, without rank-Z transform of trait
    FULL_MODEL_NRZ <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + Sex + AgeAtTest + AgeAtTest:Sex,
        REML = FALSE,
        data = RAW_DATA
      )),
      warning = war
    )
    if(is.character(FULL_MODEL_NRZ)){
      cat('WARNING: Non-rankZ full model throws a warning, refitting with BOBYQA optimization\n')
      FULL_MODEL_NRZ <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + clumps + Sex + AgeAtTest + AgeAtTest:Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA
      ))
    }
    
  } else{
    
    # Fit full model of raw data, but with the trait rank-Z transformed
    FULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest + AgeAtTest:Sex,
        REML = FALSE,
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    # Fit null model of raw data, but with the trait rank-Z transformed
    NULL_MODEL <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest,
        REML = FALSE,
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      )),
      warning = war
    )
    if(is.character(FULL_MODEL) | is.character(NULL_MODEL)){
      cat('WARNING: Full and/or null model throws a warning, refitting with BOBYQA optimization\n')
      if(is.character(FULL_MODEL)){
        cat('Full model throws a warning\n')
      }
      if(is.character(NULL_MODEL)){
        cat('Null model throws a warning\n')
      }
      # Refit full and null models using BOBYQA optimizer
      FULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest + AgeAtTest:Sex,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
      NULL_MODEL <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA %>% mutate(Y = rankZ(Y))
      ))
    }
    
    # Fit full model of raw data, without rank-Z transform of trait
    FULL_MODEL_NRZ <- tryCatch(
      suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest + AgeAtTest:Sex,
        REML = FALSE,
        data = RAW_DATA
      )),
      warning = war
    )
    if(is.character(FULL_MODEL_NRZ)){
      cat('Non-rankZ full model throws a warning, refitting with BOBYQA optimization\n')
      FULL_MODEL_NRZ <- suppressMessages(lmer(
        Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa"),
        data = RAW_DATA
      ))
      
    }
  }
  
  # Perform 1000 parametric bootstraps or enough until 20 simulated LRT 
  # statistics greater than the observed statistic are generated. This is
  # done in batches of 50, so it is possible to actually generate more than
  # 20 extreme LRT statistics
  TEST_RESULTS <- seqPBmodcomp(
    FULL_MODEL,
    NULL_MODEL,
    h = 20,
    nsim = NUMBER_OF_SIM_TO_PERFORM,
    seed = 19940418,
    cl = CLUSTER
  )
  
  AGESEX_INTERACTION_EFFECTS <- rbind(
    AGESEX_INTERACTION_EFFECTS,
    data.frame(
      Trait = TRAIT,
      PValue_LRT = as.numeric(TEST_RESULTS$test['LRT', 'p.value']),
      PValut_PBS = as.numeric(TEST_RESULTS$test['PBtest', 'p.value']),
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
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.75075 
# DOT        (Intercept) 0.22860 
# Generation (Intercept) 0.00000 
# Residual               0.61315 
# 
# -------------------- # 
#   
#   ch
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.81288 
# DOT        (Intercept) 0.21414 
# Generation (Intercept) 0.17980 
# Residual               0.45414 
# 
# -------------------- # 
#   
#   chcm
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.52512 
# DOT        (Intercept) 0.62472 
# Generation (Intercept) 0.29487 
# Residual               0.38845 
# 
# -------------------- # 
#   
#   hdw
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.63876 
# DOT        (Intercept) 0.42747 
# Generation (Intercept) 0.17114 
# Residual               0.58130 
# 
# -------------------- # 
#   
#   mcv
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.60847 
# DOT        (Intercept) 0.44583 
# Generation (Intercept) 0.21226 
# Residual               0.46095 
# 
# -------------------- # 
#   
#   rdw
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.72704 
# DOT        (Intercept) 0.18751 
# Generation (Intercept) 0.10241 
# Residual               0.62780 
# 
# -------------------- # 
#   
#   hct
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.51788 
# DOT        (Intercept) 0.43878 
# Generation (Intercept) 0.23873 
# Residual               0.52639 
# 
# -------------------- # 
#   
#   hgb
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.70210 
# DOT        (Intercept) 0.17945 
# Generation (Intercept) 0.13019 
# Residual               0.62767 
# 
# -------------------- # 
#   
#   wbc
# Groups     Name        Std.Dev.  
# MouseID    (Intercept) 7.4362e-01
# DOT        (Intercept) 2.4427e-01
# Generation (Intercept) 1.8365e-05
# Residual               6.1353e-01
# 
# -------------------- # 
#   
#   nlr
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.67948 
# DOT        (Intercept) 0.11698 
# Generation (Intercept) 0.00000 
# Residual               0.60382 
# 
# -------------------- # 
#   
#   n.lymph
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.73556 
# DOT        (Intercept) 0.23211 
# Generation (Intercept) 0.00000 
# Residual               0.62575 
# 
# -------------------- # 
#   
#   n.neut
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.64446 
# DOT        (Intercept) 0.19387 
# Generation (Intercept) 0.05263 
# Residual               0.62476 
# 
# -------------------- # 
#   
#   n.mono
# WARNING: Full and/or null model throws a warning, refitting with BOBYQA optimization
#   Null model throws a warning
#     In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#       Model failed to converge with max|grad| = 0.00548717 (tol = 0.002, component 1)
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.62391 
# DOT        (Intercept) 0.27479 
# Generation (Intercept) 0.00000 
# Residual               0.64125 
# 
# -------------------- # 
#   
#   n.eos
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.60943 
# DOT        (Intercept) 0.15978 
# Generation (Intercept) 0.00000 
# Residual               0.61723 
# 
# -------------------- # 
#   
#   plt
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.635822
# DOT        (Intercept) 0.122544
# Generation (Intercept) 0.070172
# Residual               0.579655
# 
# -------------------- # 
#   
#   mpv
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.30493 
# DOT        (Intercept) 0.64238 
# Generation (Intercept) 0.41517 
# Residual               0.38130 
# 
# -------------------- # 
#   
#   mpm
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.70165 
# DOT        (Intercept) 0.42038 
# Generation (Intercept) 0.29594 
# Residual               0.45944 
# 
# -------------------- # 



# --------------------------- Checking n.mono -------------------------------- #

# Check that BOBYQA fittings of n.mono don't also throw convergence 
# warnings, they do not.

TRAIT <- 'n.mono'

RAW_DATA <- PhenoData_Raw
RAW_DATA[['Y']] <- RAW_DATA[[TRAIT]]
RAW_DATA <- RAW_DATA %>% 
  mutate(
    AgeAtTest = standardise(AgeAtTest),
    clumps = standardise(log(clumps)),
    Y = standardise(Y)
  ) %>% 
  select(
    MouseID, Generation, DOT, clumps, AgeAtTest, Timepoint, Sex, Y
  )

FULL_MODEL <- suppressMessages(lmer(
  Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest + AgeAtTest:Sex,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa"),
  data = RAW_DATA %>% mutate(Y = rankZ(Y))
))
NULL_MODEL <- suppressMessages(lmer(
  Y ~ (1|Generation) + (1|MouseID) + (1|DOT) + Sex + AgeAtTest,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa"),
  data = RAW_DATA %>% mutate(Y = rankZ(Y))
))
VarCorr(FULL_MODEL)
# Groups     Name        Std.Dev.
# MouseID    (Intercept) 0.62391 
# DOT        (Intercept) 0.27479 
# Generation (Intercept) 0.00000 
# Residual               0.64125 

rm(TRAIT, RAW_DATA, FULL_MODEL, NULL_MODEL)

# ---------------------------------------------------------------------------- #


stopCluster(CLUSTER); rm(CLUSTER)

write_csv(
  AGESEX_INTERACTION_EFFECTS,
  'tables/lifespan_aging_hematology/agesex_interaction_effects_on_hematology_traits.csv'
)

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
