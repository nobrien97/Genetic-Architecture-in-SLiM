# Stats: 
# locisigma + pleiorate + rwide + delmu (mu) + model type on variance, covariance, distance
# also the above on allelic effect distributions: means, variance, kurtosis
# kurtosis gives an idea of relative likelihood of getting a rare allele



# Tables
# 1) Regression tables for Va, cov(x,y), sqrt(distance)
# 2) Regression tables for kurtosis, mean, variance of allele frequencies

# Set constants etc. import libraries


source("../../AIM1/AIM-1_R/src_G_mat.R")
source("../../AIM1/AIM-1_R/src_plot.R")

# Set the seed

set.seed(873662137) # sampled using sample(1:2147483647, 1)

library(plyr)
library(tidyverse)


# Import data

d_eucdist_c <- readRDS("d_eucdist_c.RDS")

d_raw_c <- readRDS("d_raw_c.RDS")

d_mean_eucdist <- readRDS("d_mean_eucdist.RDS")

d_mean_var <- readRDS("d_mean_var.RDS")

d_raw_end <- readRDS("d_raw_end.RDS")

levels(d_eucdist_c$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")
levels(d_eucdist_fingen$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")
levels(d_raw_c$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")
levels(d_raw_end$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")



levels(d_eucdist_c$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_eucdist_fingen$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_raw_c$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_raw_end$tau.cat) <- c("Null", "Strong", "Medium", "Weak")

# Comparing between models: Add a new column which is the interaction between selection strength and delmu (actually just mu)

d_eucdist_c$COA.cat <- interaction(d_eucdist_c$delmu.cat, d_eucdist_c$tau.cat)

levels(d_eucdist_c$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                                 "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                                 "Other", "Other")

d_eucdist_fingen <- d_eucdist_c[d_eucdist_c$gen == 150000,]


d_raw_end$COA.cat <- interaction(d_raw_end$delmu.cat, d_raw_end$tau.cat)

levels(d_raw_end$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                               "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                               "Other", "Other")

d_raw_c$COA.cat <- interaction(d_raw_c$delmu.cat, d_raw_c$tau.cat)

levels(d_raw_c$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                             "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                             "Other", "Other")

d_mean_eucdist$COA.cat <- interaction(d_mean_eucdist$delmu.cat, d_mean_eucdist$tau.cat)

levels(d_mean_eucdist$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                                    "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                                    "Other", "Other")

d_mean_var$COA.cat <- interaction(d_mean_var$delmu.cat, d_mean_var$tau.cat)

levels(d_mean_var$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                                "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                                "Other", "Other")


d_muts_stats <- read.csv("d_muts_stats.csv")

# Write these data frames to .csv

d_eucdist_c <- write.csv(d_eucdist_c, "d_eucdist_c.csv", row.names = F)

d_raw_c <- write.csv(d_raw_c, "d_raw_c.csv", row.names = F)

d_mean_eucdist <- write.csv(d_mean_eucdist, "d_mean_eucdist.csv", row.names = F)

d_mean_var <- write.csv(d_mean_var, "d_mean_var.csv", row.names = F)

d_raw_end <- write.csv(d_raw_end, "d_raw_end.csv", row.names = F)

d_eucdist_fingen <- write.csv(d_eucdist_fingen, "d_eucdist_fingen.csv", row.names = F)

###################################################################
# Mean variance at generation 100,000 (150,000 including burn-in)


lm_var_end <- lm(sqrtvar ~ (rwide + locisigma + COA.cat + pleiorate)^2,
                       data = d_raw_end[d_raw_end$COA.cat != "Other",]) # Not looking at the in-between models for now

plot(lm_var_end)


# Use lm_robust to adjust for non-uniform errors, non-normality accounted for by sample size
library(estimatr)
library(xtable)
d_raw_end$sqrtvar <- sqrt(d_raw_end$varmean)

lm_var_end <- lm_robust(sqrtvar ~ (COA.cat + locisigma)^2,
                       data = d_raw_end[d_raw_end$COA.cat != "Other",])

summary(lm_var_end)

# Prepare table for output

var_end_table <- lm_var_end %>% tidy %>% xtable()

print.xtable(var_end_table, type = "latex", file = "var_end_table.tex")



library(emmeans)

emmeans(lm_var_end, pairwise ~ COA.cat | pleiorate)

###################################################################
# Mean covariance at generation 100,000 (150,000 including burn-in)


lm_cov_end <- lm(covmean ~ (rwide + locisigma + COA.cat + pleiorate)^2,
                 data = d_raw_end[d_raw_end$COA.cat != "Other",]) # Not looking at the in-between models for now

plot(lm_cov_end)


# Use lm_robust to adjust for non-uniform errors, non-normality accounted for by sample size
library(estimatr)
library(xtable)

lm_cov_end <- lm_robust(covmean ~ (rwide + locisigma + COA.cat + pleiorate)^2,
                        data = d_raw_end[d_raw_end$COA.cat != "Other",])

# Prepare table for output

cov_end_table <- lm_cov_end %>% tidy %>% xtable()

print.xtable(cov_end_table, type = "latex", file = "cov_end_table.tex")

summary(lm_cov_end)

############################################################
# Distance from optimum

###################################################################
# Mean variance at generation 100,000 (150,000 including burn-in)


lm_dist_end <- lm(distance ~ (rwide + locisigma + COA.cat + pleiorate)^2,
                 data = d_eucdist_fingen[d_eucdist_fingen$COA.cat != "Other",]) # Not looking at the in-between models for now

plot(lm_dist_end)


# Use lm_robust to adjust for non-uniform errors, non-normality accounted for by sample size
library(estimatr)
library(xtable)

lm_dist_end <- lm_robust(distance ~ (rwide + locisigma + COA.cat + pleiorate)^2,
                        data = d_eucdist_fingen[d_eucdist_fingen$COA.cat != "Other",])

summary(lm_dist_end)

# Prepare table for output

dist_end_table <- lm_dist_end %>% tidy %>% xtable()

print.xtable(dist_end_table, type = "latex", file = "dist_end_table.tex")


########################################################################
# Allele frequencies: means, variance, kurtosis

# Means

lm_mean_muts <- lm(mean ~ (rwide + locisigma + COA.cat + pleiorate)^2,
                  data = d_muts_stats[d_muts_stats$COA.cat != "Other",]) # Not looking at the in-between models for now

plot(lm_mean_muts)


# Use lm_robust to adjust for non-uniform errors, non-normality accounted for by sample size
library(estimatr)
library(xtable)

lm_mean_muts <- lm_robust(mean ~ (rwide + locisigma + COA.cat + pleiorate)^2,
                         data = d_muts_stats[d_muts_stats$COA.cat != "Other",])

summary(lm_mean_muts)

# Prepare table for output

mean_muts_table <- lm_mean_muts %>% tidy %>% xtable()

print.xtable(mean_muts_table, type = "latex", file = "mean_muts_table.tex")


###########
# Variance

lm_var_muts <- lm(var ~ (rwide + locisigma + COA.cat + pleiorate)^2,
                   data = d_muts_stats[d_muts_stats$COA.cat != "Other",]) # Not looking at the in-between models for now

plot(lm_var_muts)


# Use lm_robust to adjust for non-uniform errors, non-normality accounted for by sample size
library(estimatr)
library(xtable)

lm_var_muts <- lm_robust(var ~ (rwide + locisigma + COA.cat + pleiorate)^2,
                          data = d_muts_stats[d_muts_stats$COA.cat != "Other",])

summary(lm_var_muts)

# Prepare table for output

var_muts_table <- lm_var_muts %>% tidy %>% xtable()

print.xtable(var_muts_table, type = "latex", file = "var_muts_table.tex")


############################
# Kurtosis

lm_kurt_muts <- lm(kurt ~ (rwide + locisigma + COA.cat + pleiorate)^2,
                   data = d_muts_stats[d_muts_stats$COA.cat != "Other",]) # Not looking at the in-between models for now

plot(lm_kurt_muts)


# Use lm_robust to adjust for non-uniform errors, non-normality accounted for by sample size
library(estimatr)
library(xtable)


lm_kurt_muts <- lm_robust(kurt ~ (rwide + locisigma + COA.cat + pleiorate)^2,
                          data = d_muts_stats[d_muts_stats$COA.cat != "Other",])

summary(lm_kurt_muts)

lm_kurt_muts_nonull <- lm_robust(kurt ~ (rwide + locisigma + COA.cat + pleiorate)^2,
                          data = d_muts_stats[d_muts_stats$COA.cat != "Other" & d_muts_stats$COA.cat != "Null",])

summary(lm_kurt_muts_nonull)

# Prepare table for output

kurt_muts_table <- lm_kurt_muts_nonull %>% tidy %>% xtable()

print.xtable(kurt_muts_table, type = "latex", file = "kurt_muts_table.tex")

