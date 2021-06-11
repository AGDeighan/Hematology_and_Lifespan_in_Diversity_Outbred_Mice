# Creating QTL analysis datasets for lifespan scans
# 2021-06-11

################################################################################
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
#### Libraries etc ###################

library(tidyverse)

#####

################################################################################
#### Helper functions ###################

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

create_QTL_dataset <- function(
  dataframe, 
  id_col,
  resp_col,
  factor_cov_col = NULL,
  numeric_cov_col = NULL,
  int_cov_col = NULL,
  type = 'phenotype',
  name = NA
){
  
  require(dplyr)
  
  dataframe <- dataframe %>% 
    filter_at(vars(all_of(id_col), all_of(factor_cov_col), all_of(numeric_cov_col)), all_vars(!is.na(.)))
  
  ################################################################################
  # Create data dictionary
  
  DATA_DICTIONARY <- tibble(data_name = c(id_col, resp_col),
                            short_name = c(id_col, resp_col),
                            R_name = c(id_col, resp_col),
                            description = as.character(NA),
                            units = c('N/A', rep(NA, length(resp_col))),
                            is_id = c(TRUE, rep(FALSE, length(resp_col))),
                            category = c('ID', rep(NA, length(resp_col))),
                            R_category = c('ID', rep(NA, length(resp_col))),
                            is_numeric = c(FALSE, rep(TRUE, length(resp_col))),
                            is_date = FALSE,
                            is_factor = FALSE,
                            factor_levels = as.character(NA),
                            is_covar = FALSE,
                            is_pheno = c(FALSE, rep(TRUE, length(resp_col))),
                            is_derived = FALSE,
                            omit = FALSE,
                            use_covar = c(NA, rep(paste(c(factor_cov_col, numeric_cov_col), collapse = ':'), length(resp_col)))
  )
  
  ################################################################################
  # Create covariate matrix
  
  if(!is.null(factor_cov_col) | !is.null(numeric_cov_col)){
    # Convert Diet to continuous
    COV_MATRIX <- dataframe %>% 
      select(all_of(factor_cov_col), all_of(numeric_cov_col))
    
    
    FORMULA <- paste(c(factor_cov_col, numeric_cov_col), collapse = ' + ')
    FORMULA <- paste0('~ ', FORMULA)
    FORMULA <- as.formula(FORMULA)
    COV_MATRIX <- model.matrix(FORMULA,
                               data = COV_MATRIX)[,-1,drop = FALSE]
    rownames(COV_MATRIX) <- dataframe[[id_col]]
  } else{
    COV_MATRIX <- NULL
  }
  
  
  
  ################################################################################
  # Create covariate factor dataframe
  
  COV_FACTORS <- data.frame(column_name = vector('character'),    
                            display_name = vector('character'),   
                            int.covar = vector('character'),      
                            lod.peaks = vector('character'),                         
                            covar.name = vector('character'))
  
  for(covar in unique(c(factor_cov_col, numeric_cov_col, int_cov_col))){
    NAMES <- colnames(COV_MATRIX)[which(substr(colnames(COV_MATRIX), 1, nchar(covar)) == covar)]
    NAMES <- paste(NAMES, collapse = ',')
    
    NEW_ROW <- data.frame(column_name = covar,    
                          display_name = covar,   
                          int.covar = ifelse(covar %in% int_cov_col,
                                             ifelse(covar %in% factor_cov_col,
                                                    'factor', 'numeric'),
                                             as.character(NA)),      
                          lod.peaks = ifelse(covar %in% int_cov_col,
                                             paste0(tolower(covar), '_int'),
                                             as.character(NA)),                         
                          covar.name = NAMES)
    
    COV_FACTORS <- rbind(COV_FACTORS, NEW_ROW)
  }; rm(covar, NAMES, NEW_ROW)
  
  
  ################################################################################
  # Specify datatype
  
  TYPE <- type
  
  
  ################################################################################
  # Specify display name
  
  NAME <- name
  
  
  ################################################################################
  # Create list of LOD peaks
  
  PEAKS <- list(additive = NA)
  
  for(covar in COV_FACTORS$lod.peaks[!is.na(COV_FACTORS$lod.peaks)]){
    PEAKS <- append(PEAKS, NA)
    names(PEAKS)[length(PEAKS)] <- covar
  }; rm(covar)
  
  
  ################################################################################
  # Create phenotype dataframe
  
  PHENOS <- dataframe %>% 
    select(all_of(id_col), all_of(resp_col)) %>% 
    as.data.frame()
  
  rownames(PHENOS) <- dataframe[[id_col]]
  
  
  ################################################################################
  # Create sample/mouse description dataframe
  
  SAMPLES_DF <- data.frame(mouse.id = dataframe[[id_col]],
                           stringsAsFactors = FALSE)
  
  for(covar in unique(c(factor_cov_col, numeric_cov_col))){
    SAMPLES_DF[[covar]] <- dataframe[[covar]]
  }; rm(covar)
  
  rownames(SAMPLES_DF) <- dataframe[[id_col]]
  
  ################################################################################
  # Create dataset list
  
  DATASET <- list(annots = as.data.frame(DATA_DICTIONARY),
                  covar = COV_MATRIX,
                  covar.factors = COV_FACTORS,
                  datatype = TYPE,
                  display.name = NAME,
                  lod.peaks = PEAKS,
                  pheno = PHENOS,
                  samples = SAMPLES_DF)
  
  return(DATASET)
}

#####

################################################################################
## load data etc #########################################################

load('data/processed/phenotypic/CBC_AnalysisEnvironment.RData')

#####


################################################################################
## Calculate sample sizes ###########################

(
  SAMPLE_SIZES <- PhenoData_ClumpDateAdj %>% 
    filter(Timepoint == '7 Months') %>% 
    select(
      MouseID, Lifespan, hdw, rdw, nlr
    ) %>% 
    gather(key = Trait, value = V, Lifespan, hdw, rdw, nlr) %>% 
    group_by(Trait) %>% 
    summarise(
      N = sum(!is.na(V))
    ) %>% 
    ungroup() %>% 
    data.frame()
)
#      Trait   N
# 1      hdw 526
# 2 Lifespan 526
# 3      nlr 526
# 4      rdw 526

#####


################################################################################
## Widen dataframe and rank-Z transform CBC traits ##########################

DATA <- PhenoData_ClumpDateAdj %>% 
  select(
    MouseID, Generation, Sex, Lifespan, Timepoint, rbc:mpm
  ) %>% 
  pivot_longer(
    cols = rbc:mpm,
    names_to = 'K',
    values_to = 'V'
  ) %>% 
  mutate(
    Timepoint = gsub(' Months', '', Timepoint),
    Timepoint = ifelse(nchar(Timepoint) < 2, paste0(0, Timepoint), Timepoint),
    K = paste0(K, '_', Timepoint)
  ) %>% 
  select(-Timepoint) %>% 
  pivot_wider(
    names_from = 'K',
    values_from = 'V'
  ) %>% 
  select(
    MouseID, Generation, Sex, 
    Lifespan, 
    all_of(paste0(rep(CBCphenos, each = 3), '_', rep(c('07', '13', '19'), length(CBCphenos))))
  ) %>% 
  mutate_at(
    vars(rbc_07:mpm_19), rankZ
  )

#####


################################################################################
## Create unconditioned dataset for qtl2 analysis ############################

DS_UNDCOND <- create_QTL_dataset(
  dataframe = DATA,
  id_col = 'MouseID',
  resp_col = names(DATA)[-c(1:3)],
  factor_cov_col = c('Sex', 'Generation'),
  numeric_cov_col = NULL,
  int_cov_col = NULL,
  type = 'phenotype',
  name = 'na'
)

DS_UNDCOND %>% 
  saveRDS('data/processed/qtl_analysis_data/datasets/LifespanHema_Dataset_CovSexGen.rds')

#####

################################################################################
## Conditioned on rdw ############################

DS_RDW_COND <- create_QTL_dataset(
  dataframe = DATA,
  id_col = 'MouseID',
  resp_col = c('Lifespan'),
  factor_cov_col = c('Sex', 'Generation'),
  numeric_cov_col = c('rdw_07'),
  int_cov_col = NULL,
  type = 'phenotype',
  name = 'na'
)

DS_RDW_COND %>% 
  saveRDS(
    'data/processed/qtl_analysis_data/datasets/Lifespan_Dataset_CovSexGenRDW.rds'
  )

#####

################################################################################
## Conditioned on hdw ############################

DS_HDW_COND <- create_QTL_dataset(
  dataframe = DATA,
  id_col = 'MouseID',
  resp_col = c('Lifespan'),
  factor_cov_col = c('Sex', 'Generation'),
  numeric_cov_col = c('hdw_07'),
  int_cov_col = NULL,
  type = 'phenotype',
  name = 'na'
)

DS_HDW_COND %>% 
  saveRDS(
    'data/processed/qtl_analysis_data/datasets/Lifespan_Dataset_CovSexGenHDW.rds'
  )

#####

################################################################################
## Conditioned on nlr ############################

DS_NLR_COND <- create_QTL_dataset(
  dataframe = DATA,
  id_col = 'MouseID',
  resp_col = c('Lifespan'),
  factor_cov_col = c('Sex', 'Generation'),
  numeric_cov_col = c('nlr_07'),
  int_cov_col = NULL,
  type = 'phenotype',
  name = 'na'
)

DS_NLR_COND %>% 
  saveRDS(
    'data/processed/qtl_analysis_data/datasets/Lifespan_Dataset_CovSexGenNLR.rds'
  )

#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####
