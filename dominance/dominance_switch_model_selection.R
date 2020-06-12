# This script performs model selection on variables that affect dominant species switching
#
###########################################################################################
###########################################################################################
# Load libraries
library(tidyverse)
library(glmnet)
library(ggplot2)

## -- REPOSITORY DIR -- ##
setwd("/Users/avahoffman/Dropbox/Research/Grazing_consortium/2020/GEx_SEV2019")


dat <-
  read.csv("community_difference_allmetrics_siteavg_12June2020b.csv")

dat_sub <- 
  dat %>%
  select(c(diff_sp,
           precip,
           CAM,
           bio1,
           NMDS1,
           N.deposition1993,
           PhotoMix,
           ALLC3)) %>%
  drop_na()

dat_y <-
  as.numeric((dat_sub %>%
                select(diff_sp))[, 1])

dat_x <- dat_sub %>%
  select(-c(diff_sp))

#find mean, std in dat_x
mean_x <- sapply(dat_x, mean)
sd_x <- sapply(dat_x, sd)

x_scaled <- scale(dat_x)


lambda_seq <- 10 ^ seq(2,-2, by = -.1)

cv_output <- cv.glmnet(x_scaled,
                       dat_y,
                       alpha = 1,
                       lambda = lambda_seq,
                       nfolds = 5)

best_lam <- cv_output$lambda.min
best_lam

lasso_out <-
  glmnet(x_scaled,
         dat_y,
         family = "gaussian",
         alpha = 1,
         lambda = best_lam)
