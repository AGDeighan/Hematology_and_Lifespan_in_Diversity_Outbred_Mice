# 2021-01-07
################################################################################
#
#   This script creates survival curves by sex
#
#   Author:  Andrew Deighan
#   E-mail: andrew.deighan@jax.org
#
################################################################################


options(na.action = 'na.exclude')
options(stringsAsFactors = FALSE)

################################################################################
## libraries etc #########################################################
library(tidyverse)
library(survival)
library(survminer)
library(lme4)
library(pbkrtest)

#####


################################################################################
## helper functions #########################################################

#####

################################################################################
## load data etc #########################################################

load('data/processed/phenotypic/CBC_AnalysisEnvironment.RData')

#####

################################################################################
## Sex Effects ####

DATA <- PhenoData_Raw %>% 
  select(MouseID, Generation, Sex, Lifespan) %>% 
  group_by(MouseID, Generation, Sex, Lifespan) %>% 
  summarise() %>% 
  ungroup()

(
  POOLED_SUMMARY_STATS <- DATA %>% 
    summarise(
      Mean = round(mean(Lifespan), digits = 2),
      Median = round(median(Lifespan), digits = 2), 
      p25 = round(quantile(Lifespan, 0.25), digits = 2),
      p75 = round(quantile(Lifespan, 0.75), digits = 2),
      Range = paste0(round(range(Lifespan), digits = 2), collapse = ' - '),
      Variance = round(var(Lifespan), digits = 2),
      SD = round(sd(Lifespan), digits = 2),
      MAD = round(mad(Lifespan), digits = 2),
      IQR = round(IQR(Lifespan), digits = 2)
    ) %>% 
    ungroup() %>% 
    data.frame()
)
#    Mean Median   p25   p75        Range Variance  SD  MAD   IQR
# 1 27.51  27.93 21.64 33.42 8.26 - 55.59    67.21 8.2 8.97 11.78

pdf(
  'figures/lifespan_aging_and_hematology/Distribution_of_lifespan.pdf',
  width = 6.5, height = 4
)
DATA %>% 
  ggplot() +
  theme_minimal(
    base_size = 13
  ) +
  geom_histogram(
    aes(x = Lifespan, fill = Sex),
    binwidth = 2
  ) +
  geom_vline(
    xintercept = POOLED_SUMMARY_STATS$Mean,
    linetype = 1
  ) +
  geom_vline(
    xintercept = POOLED_SUMMARY_STATS$Median,
    linetype = 2
  ) +
  geom_vline(
    xintercept = c(POOLED_SUMMARY_STATS$p25, POOLED_SUMMARY_STATS$p75),
    linetype = 3
  ) +
  scale_fill_manual(
    values = c('darkorange3', 'cadetblue')
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 10),
    minor_breaks = seq(0, 100, 2),
    limits = c(0, 56),
    exp = c(0, 0)
  ) +
  theme(
    legend.position = 'top'
  ) +
  labs(
    x = 'Lifespan (months)',
    y = 'Frequency'
  )
dev.off()  

(
  SEX_SPECIFIC_SUMMARY_STATS <- DATA %>% 
    group_by(Sex) %>% 
    summarise(
      Mean = round(mean(Lifespan), digits = 2),
      Median = round(median(Lifespan), digits = 2), 
      p25 = round(quantile(Lifespan, 0.25), digits = 2),
      p75 = round(quantile(Lifespan, 0.75), digits = 2),
      Range = paste0(round(range(Lifespan), digits = 2), collapse = ' - '),
      Variance = round(var(Lifespan), digits = 2),
      SD = round(sd(Lifespan), digits = 2),
      MAD = round(mad(Lifespan), digits = 2),
      IQR = round(IQR(Lifespan), digits = 2)
    ) %>% 
    ungroup() %>% 
    data.frame()
)
#      Sex  Mean Median   p25   p75        Range Variance   SD  MAD   IQR
# 1 Female 28.03  28.75 22.40 33.78 8.85 - 46.22    62.26 7.89 8.63 11.38
# 2   Male 27.00  26.78 21.01 33.31 8.26 - 55.59    71.78 8.47 8.78 12.29

# Test of equal variances (they don't look significantly different)
var.test(Lifespan ~ Sex, data = DATA)
#   F test to compare two variances
# 
# data:  Lifespan by Sex
# F = 0.86737, num df = 256, denom df = 263, p-value = 0.2529                   # Variance of lifespan does not differ significantly by sex
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.6797852 1.1071759
# sample estimates:
# ratio of variances 
# 0.867373 

(T_TEST <- t.test(Lifespan ~ Sex, var.equal = TRUE, data = DATA))
# data:  Lifespan by Sex
# t = 1.427, df = 519, p-value = 0.1542                                         # Sex does not have a significant effect on lifespan
# alternative hypothesis: true difference in means is not equal to 0            
# 95 percent confidence interval:
#   -0.3858387  2.4341227
# sample estimates:
#   mean in group Female   mean in group Male 
#               28.02798             27.00384 
(PVAL_TTEST <- T_TEST$p.value)
# [1] 0.1541953

F_TEST <- pbkrtest::KRmodcomp(
  lme4::lmer(
    Lifespan ~ Sex + (1|Generation),
    data = DATA
  ),
  lme4::lmer(
    Lifespan ~ (1|Generation),
    data = DATA
  )
)
(PVAL_FTEST <- F_TEST$test$p.value[1])
# [1] 0.1545729                                                                 # F-test gives same result as t-test

SURV_MOD <- survfit(
  Surv(Lifespan - 7, Died) ~ Sex, 
  data = DATA %>% 
    mutate(Died = TRUE)
)
SURV_PLOT <- ggsurvplot(
  fit = SURV_MOD,
  size = 1.25,
  censor.shape = 124,
  conf.int = FALSE,
  pval = FALSE,
  palette = c('darkorange3', 'cadetblue'),
  alpha = 0.9,
  legend = c(0.8, 0.8),
  title = NULL,
  legend.title = '',
  legend.labs = levels(DATA$Sex),
  xlab = 'Age (months)',
  surv.median.line = 'n',
  xlim = c(6.75, 57) - 7,
  font.x = c(11, "plain", "black"),
  font.y = c(11, "plain", "black"),
  font.tickslab = c(9, "plain", "black"),
  data = DATA %>% 
    mutate(Died = TRUE)
)
X <- SURV_PLOT$plot +
  annotate(
    'rect',
    xmin = 21.64 - 7, xmax = 33.42 - 7,
    ymin = -0.01, ymax = 1.01,
    fill = 'grey90',
    alpha = 0.5
  ) +
  geom_vline(xintercept = 7 - 7, linetype = 3) +
  geom_vline(xintercept = 13 - 7, linetype = 3) +
  geom_vline(xintercept = 19 - 7, linetype = 3) +
  geom_vline(xintercept = 27.93 - 7, linetype = 1, size = 1, color = 'grey50') +
  annotate(
    'text',
    label = paste0('P-Value = ', formatC(PVAL_TTEST, digits = 3, format = 'f')),
    x = 40,
    y = 0.25,
    hjust = 0
  ) +
  scale_x_continuous(
    breaks = seq(0,60,6) - 7,
    labels = c('', '', 12, '', 24, '', 36, '', 48, '', 60),
    exp = c(0,0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    limits = c(-0.01, 1.01),
    exp = c(0,0)
  ) +
  theme(
    legend.position = 'top'
    ) +
  labs(
    y = 'Survival Probability'
  )

Y <- X$layers
Y[[7]] <- X$layers[[1]]
Y[[8]] <- X$layers[[2]]
Y[[1]] <- X$layers[[3]]
Y[[2]] <- X$layers[[4]]
Y[[3]] <- X$layers[[5]]
Y[[4]] <- X$layers[[6]]
Y[[5]] <- X$layers[[7]]
Y[[6]] <- X$layers[[8]]
X$layers <- Y
X

pdf(
  'figures/lifespan_aging_and_hematology/Sex_surv_curve.pdf',
  width = 6.5, height = 4
)
plot(X)
dev.off()

#####


################################################################################
##   ##########################################

#####

################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####
