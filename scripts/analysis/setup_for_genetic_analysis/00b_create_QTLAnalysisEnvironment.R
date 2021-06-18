# 2021-06-11

################################################################################
## Create environment for QTL analysis #####################################

# Load the "datasets"
dataset.LifespanHemaCovSexGen.phenotype <- readRDS('data/processed/qtl_analysis_data/datasets/LifespanHema_Dataset_CovSexGen.rds')
dataset.LifespanCovSexGenRDW.phenotype <- readRDS('data/processed/qtl_analysis_data/datasets/Lifespan_Dataset_CovSexGenRDW.rds')
dataset.LifespanCovSexGenHDW.phenotype <- readRDS('data/processed/qtl_analysis_data/datasets/Lifespan_Dataset_CovSexGenHDW.rds')
dataset.LifespanCovSexGenNLR.phenotype <- readRDS('data/processed/qtl_analysis_data/datasets/Lifespan_Dataset_CovSexGenNLR.rds')

# Load the genoprobs
genoprobs <- readRDS('data/processed/genetic/AlleleProbs_8State.rds')

# Calculate leave-one-chromosome-out the kinship matrix
K <- qtl2::calc_kinship(probs = genoprobs, type = 'loco', cores = 6)

# Load marker map as list
map <- readRDS('data/processed/genetic/PhysicalMap_List.rds')

# Load marker map as dataframe
markers <- readRDS('data/processed/genetic/PhysicalMap_Dataframe.rds')

# Save environment
save(
  list = ls(), 
  file = 'data/processed/qtl_analysis_data/environments/LifespanHema_Environment.RData'
)

#####

################################################################################
## clear workspace ###########################################

rm(list = ls())
graphics.off()

#####