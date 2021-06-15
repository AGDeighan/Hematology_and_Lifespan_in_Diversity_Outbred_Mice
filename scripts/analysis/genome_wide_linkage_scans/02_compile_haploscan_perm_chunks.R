# 2021-01-15



################################################################################
## load libraries ###########################################

options(stringsAsFactors = FALSE)
library(qtl2)

#####


################################################################################
## accessory functions #############

read_perm_file <- function(file, path){
  return(readRDS(paste0(path, file)))
}

#####


################################################################################
## Set up #############

FILES <- list.files(
  path = 'results/genome_wide_linkage_scans/haploscan_perm_chunks/', 
  pattern = 'rds$'
)

#####



################################################################################
## Unconditioned #############

PERM_SCANS <- lapply(
  FILES,
  read_perm_file,
  path = 'results/genome_wide_linkage_scans/haploscan_perm_chunks/'
)

COMBINED_PERM_SCANS <- c(
  PERM_SCANS[[1]], PERM_SCANS[[2]]
)
for(i in 3:length(PERM_SCANS)){
  COMBINED_PERM_SCANS <- c(
    COMBINED_PERM_SCANS, PERM_SCANS[[i]]
  )
}

PERMS <- nrow(COMBINED_PERM_SCANS)

FILENAME <- paste0(
  'results/genome_wide_linkage_scans/LifespanHema_',
  PERMS, 'Perms.rds'
)

saveRDS(
  COMBINED_PERM_SCANS,
  file = FILENAME
)

#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####
