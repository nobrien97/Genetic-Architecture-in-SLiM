
#############################################################################################################################

# Set constants etc. import libraries


source("../../AIM1/AIM-1_R/src_G_mat.R")
source("../../AIM1/AIM-1_R/src_plot.R")

# Set the seed

set.seed(873662137) # sampled using sample(1:2147483647, 1)

library(plyr)
library(tidyverse)


# Import data

# Data over time
d_eucdist_c <- readRDS("d_eucdist_c.RDS")

d_raw_c <- readRDS("d_raw_c.RDS")

# End point 'equilibrium' data
d_combined <- read.csv("d_combined.csv")

# Allelic effects data
d_muts_stats <- read.csv("d_muts_stats.csv", stringsAsFactors = T)

########################

# Fix up some factor labels d_eucdist_c and d_raw_c

levels(d_eucdist_c$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")
levels(d_raw_c$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")
levels(d_muts_stats$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")
levels(d_eucdist_c$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_raw_c$tau.cat) <- c("Null", "Strong", "Medium", "Weak")

d_eucdist_c$locisigma.cat <- factor(d_eucdist_c$locisigma.cat, levels = c("Low variance", "Medium variance", "High variance"))
d_raw_c$locisigma.cat <- factor(d_raw_c$locisigma.cat, levels = c("Low variance", "Medium variance", "High variance"))
d_combined$locisigma.cat <- factor(d_combined$locisigma.cat, levels = c("Low variance", "Medium variance", "High variance"))
d_muts_stats$locisigma.cat <- factor(d_muts_stats$locisigma.cat, levels = c("Small variance", "Medium variance", "Large variance"))

d_combined$rwide.cat <- factor(d_combined$rwide.cat, levels = c("Low", "Medium", "High"))
d_combined$pleiorate.cat <- factor(d_combined$pleiorate.cat, levels = c("Low", "Medium", "High"))
d_combined$pleiocov.cat <- factor(d_combined$pleiocov.cat, levels = c("Low", "Medium", "High"))



d_eucdist_c$COA.cat <- interaction(d_eucdist_c$delmu.cat, d_eucdist_c$tau.cat)

levels(d_eucdist_c$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                                 "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                                 "Other", "Other")

d_raw_c$COA.cat <- interaction(d_raw_c$delmu.cat, d_raw_c$tau.cat)

levels(d_raw_c$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                             "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                             "Other", "Other")

######################

# Add adapted/maladapted variable to end point data

d_combined$sqrtdist <- sqrt(d_combined$distance) 
d_combined$atopt <- d_combined$distance
d_combined[d_combined$sqrtdist < 4.0,]$atopt <- 1
d_combined[d_combined$sqrtdist >= 4.0,]$atopt <- 0
d_combined$atopt <- as.factor(d_combined$atopt)

levels(d_combined$atopt) <- c("Maladapted", "Adapted")



# Chi sq analysis

library(chisq.posthoc.test)

conttab_dist <- table(d_combined[d_combined$COA.cat != "Other",]$atopt, d_combined[d_combined$COA.cat != "Other",]$COA.cat)
chisq.test(conttab_dist)

chisq.posthoc.test(conttab_dist)


# Variance linear model

d_combined_stat <- d_combined

# Using robust SEs
library(estimatr)
library(emmeans)

lm_var_end <- lm_robust(varmean ~ rwide.cat + locisigma.cat + COA.cat + pleiorate.cat + pleiocov.cat +
                   COA.cat*rwide.cat + COA.cat*locisigma.cat + COA.cat*pleiorate.cat + COA.cat*pleiocov.cat,
                 data = d_combined_stat[d_combined_stat$COA.cat != "Other" & d_combined_stat$COA.cat != "Null" & d_combined$atopt == "Adapted",]) # Not looking at the in-between models for now

summary(lm_var_end)

# Relative importance of predictors based on R^2
lm_var_end_OLS <- lm(varmean ~ rwide.cat + locisigma.cat + COA.cat + pleiorate.cat + pleiocov.cat +
                          COA.cat*rwide.cat + COA.cat*locisigma.cat + COA.cat*pleiorate.cat + COA.cat*pleiocov.cat,
                        data = d_combined_stat[d_combined_stat$COA.cat != "Other" & d_combined_stat$COA.cat != "Null" & d_combined$atopt == "Adapted",]) # Not looking at the in-between models for now


# Relative importance based on R^2
library(relaimpo)

calc.relimp(lm_var_end_OLS)



# Estimated marginal means
es_var.m <- emmeans(lm_var_end, pairwise ~  COA.cat)
es_var.m

es_var_l.m <- emmeans(lm_var_end, pairwise ~ locisigma.cat | COA.cat)
es_var_l.m

emm_var_contr_l.m <- pairs(pairs(emmeans(lm_var_end, ~ locisigma.cat | COA.cat,
                                          at = list(locisigma.cat = c("Low variance", "High variance"),
                                                    COA.cat = c("Gaussian", "House-of-Cards")))),
                            by = NULL)

emm_var_contr_l.m

es_var_pr.m <- emmeans(lm_var_end, pairwise ~ pleiorate.cat | COA.cat)
es_var_pr.m

emm_var_contr_pr.m <- pairs(pairs(emmeans(lm_var_end, ~ pleiorate.cat | COA.cat,
                                           at = list(pleiorate.cat = c("Low", "High"),
                                                     COA.cat = c("Gaussian", "House-of-Cards")))),
                             by = NULL)

emm_var_contr_pr.m


es_var_r.m <- emmeans(lm_var_end, pairwise ~ rwide.cat | COA.cat)
es_var_r.m

emm_var_contr_r.m <- pairs(pairs(emmeans(lm_var_end, ~ rwide.cat | COA.cat,
                                          at = list(rwide.cat = c("Low", "High"),
                                                    COA.cat = c("Gaussian", "House-of-Cards")))),
                            by = NULL)

emm_var_contr_r.m


es_var_pc.m <- emmeans(lm_var_end, pairwise ~ pleiocov.cat | COA.cat)
es_var_pc.m

emm_var_contr_pc.m <- pairs(pairs(emmeans(lm_var_end, ~ pleiocov.cat | COA.cat,
                                           at = list(pleiocov.cat = c("Low", "High"),
                                                     COA.cat = c("Gaussian", "House-of-Cards")))),
                             by = NULL)

emm_var_contr_pc.m

################################################################################################################

# Covariance 


lm_cov_end <- lm_robust(covmean ~ rwide.cat + locisigma.cat + COA.cat + pleiorate.cat + pleiocov.cat +
                          COA.cat*rwide.cat + COA.cat*locisigma.cat + COA.cat*pleiorate.cat + COA.cat*pleiocov.cat,
                        data = d_combined_stat[d_combined_stat$COA.cat != "Other" & d_combined_stat$COA.cat != "Null" & d_combined$atopt == "Adapted",]) # Not looking at the in-between models for now

summary(lm_cov_end)


lm_cov_end_OLS <- lm(covmean ~ rwide.cat + locisigma.cat + COA.cat + pleiorate.cat + pleiocov.cat +
                       COA.cat*rwide.cat + COA.cat*locisigma.cat + COA.cat*pleiorate.cat + COA.cat*pleiocov.cat,
                     data = d_combined_stat[d_combined_stat$COA.cat != "Other" & d_combined_stat$COA.cat != "Null" & d_combined$atopt == "Adapted",]) # Not looking at the in-between models for now



library(relaimpo)

calc.relimp(lm_cov_end_OLS)



es_cov.m <- emmeans(lm_cov_end, pairwise ~  COA.cat)
es_cov.m

es_cov_l.m <- emmeans(lm_cov_end, pairwise ~ locisigma.cat | COA.cat)
es_cov_l.m

emm_cov_contr_l.m <- pairs(pairs(emmeans(lm_cov_end, ~ locisigma.cat | COA.cat,
                                         at = list(locisigma.cat = c("Low variance", "High variance"),
                                                   COA.cat = c("Gaussian", "House-of-Cards")))),
                           by = NULL)

emm_cov_contr_l.m

es_cov_pr.m <- emmeans(lm_cov_end, pairwise ~ pleiorate.cat | COA.cat)
es_cov_pr.m

emm_cov_contr_pr.m <- pairs(pairs(emmeans(lm_cov_end, ~ pleiorate.cat | COA.cat,
                                          at = list(pleiorate.cat = c("Low", "High"),
                                                    COA.cat = c("Gaussian", "House-of-Cards")))),
                            by = NULL)

emm_cov_contr_pr.m


es_cov_r.m <- emmeans(lm_cov_end, pairwise ~ rwide.cat) # No effects of interaction 
es_cov_r.m

emm_cov_contr_r.m <- pairs(pairs(emmeans(lm_cov_end, ~ rwide.cat | COA.cat,
                                         at = list(rwide.cat = c("Low", "High"),
                                                   COA.cat = c("Gaussian", "House-of-Cards")))),
                           by = NULL)

emm_cov_contr_r.m


es_cov_pc.m <- emmeans(lm_cov_end, pairwise ~ pleiocov.cat | COA.cat)
es_cov_pc.m

emm_cov_contr_pc.m <- pairs(pairs(emmeans(lm_cov_end, ~ pleiocov.cat | COA.cat,
                                          at = list(pleiocov.cat = c("Low", "High"),
                                                    COA.cat = c("Gaussian", "House-of-Cards")))),
                            by = NULL)

emm_cov_contr_pc.m

#####################################################################################################################

# Distance 


lm_dist_end <- lm_robust(distance ~ rwide.cat + locisigma.cat + COA.cat + pleiorate.cat + pleiocov.cat +
                                         COA.cat*rwide.cat + COA.cat*locisigma.cat + COA.cat*pleiorate.cat + COA.cat*pleiocov.cat,
                                       data = d_combined_stat[d_combined_stat$COA.cat != "Other" & d_combined_stat$COA.cat != "Null" & d_combined$atopt == "Adapted",]) # Not looking at the in-between models for now

summary(lm_dist_end)

lm_dist_end_OLS <- lm(distance ~ rwide.cat + locisigma.cat + COA.cat + pleiorate.cat + pleiocov.cat +
                       COA.cat*rwide.cat + COA.cat*locisigma.cat + COA.cat*pleiorate.cat + COA.cat*pleiocov.cat,
                     data = d_combined_stat[d_combined_stat$COA.cat != "Other" & d_combined_stat$COA.cat != "Null" & d_combined$atopt == "Adapted",]) # Not looking at the in-between models for now



library(relaimpo)

calc.relimp(lm_dist_end_OLS, rela = T)


es_dist.m <- emmeans(lm_dist_end, pairwise ~  COA.cat)
es_dist.m


es_dist_l.m <- emmeans(lm_dist_end, pairwise ~ locisigma.cat * COA.cat)
es_dist_l.m

emm_dist_contr_l.m <- pairs(pairs(emmeans(lm_dist_end, ~ locisigma.cat | COA.cat,
                                            at = list(locisigma.cat = c("Low variance", "High variance"),
                                                      COA.cat = c("Gaussian", "House-of-Cards")))),
                                                      by = NULL)

emm_dist_contr_l.m

es_dist_pr.m <- emmeans(lm_dist_end, pairwise ~ pleiorate.cat)
es_dist_pr.m

emm_dist_contr_pr.m <- pairs(pairs(emmeans(lm_dist_end, ~ pleiorate.cat | COA.cat,
                                          at = list(pleiorate.cat = c("Low", "High"),
                                                    COA.cat = c("Gaussian", "House-of-Cards")))),
                            by = NULL)

emm_dist_contr_pr.m


es_dist_r.m <- emmeans(lm_dist_end, pairwise ~ rwide.cat | COA.cat)
es_dist_r.m

emm_dist_contr_r.m <- pairs(pairs(emmeans(lm_dist_end, ~ rwide.cat | COA.cat,
                                           at = list(rwide.cat = c("Low", "High"),
                                                     COA.cat = c("Gaussian", "House-of-Cards")))),
                             by = NULL)

emm_dist_contr_r.m


es_dist_pc.m <- emmeans(lm_dist_end, pairwise ~ pleiocov.cat | COA.cat)
es_dist_pc.m

emm_dist_contr_pc.m <- pairs(pairs(emmeans(lm_dist_end, ~ pleiocov.cat | COA.cat,
                                          at = list(pleiocov.cat = c("Low", "High"),
                                                    COA.cat = c("Gaussian", "House-of-Cards")))),
                            by = NULL)

emm_dist_contr_pc.m

################################################################################################################################################

# Mutation stats: comparing means, variance, kurtosis, count



lm_mean_muts <- lm_robust(mean ~ rwide.cat + locisigma.cat + COA.cat + pleiorate.cat + pleiocov.cat +
                          COA.cat*rwide.cat + COA.cat*locisigma.cat + COA.cat*pleiorate.cat + COA.cat*pleiocov.cat,
                        data = d_muts_stats[d_muts_stats$COA.cat != "Other" & d_muts_stats$COA.cat != "Null" & d_muts_stats$Po == "Adapted",],
                        se_type = "HC3") # Not looking at the in-between models for now

summary(lm_mean_muts)

lm_mean_muts_OLS <- lm(mean ~ rwide.cat + locisigma.cat + COA.cat + pleiorate.cat + pleiocov.cat +
                            COA.cat*rwide.cat + COA.cat*locisigma.cat + COA.cat*pleiorate.cat + COA.cat*pleiocov.cat,
                          data = d_muts_stats[d_muts_stats$COA.cat != "Other" & d_muts_stats$COA.cat != "Null" & d_muts_stats$Po == "Adapted",]) # Not looking at the in-between models for now


calc.relimp(lm_mean_muts_OLS)

# Not significant, so no post hoc

##############################################################################################################

# Variance


lm_var_muts <- lm_robust(var ~ rwide.cat + locisigma.cat + COA.cat + pleiorate.cat + pleiocov.cat +
                      COA.cat*rwide.cat + COA.cat*locisigma.cat + COA.cat*pleiorate.cat + COA.cat*pleiocov.cat,
                    data = d_muts_stats[d_muts_stats$COA.cat != "Other" & d_muts_stats$COA.cat != "Null" & d_muts_stats$Po == "Adapted",],
                    se_type = "HC3") # Not looking at the in-between models for now

summary(lm_var_muts)


lm_var_muts_OLS <- lm(var ~ rwide.cat + locisigma.cat + COA.cat + pleiorate.cat + pleiocov.cat +
                         COA.cat*rwide.cat + COA.cat*locisigma.cat + COA.cat*pleiorate.cat + COA.cat*pleiocov.cat,
                       data = d_muts_stats[d_muts_stats$COA.cat != "Other" & d_muts_stats$COA.cat != "Null" & d_muts_stats$Po == "Adapted",]) # Not looking at the in-between models for now


calc.relimp(lm_var_muts_OLS)


es_varmuts.m <- emmeans(lm_var_muts, pairwise ~  COA.cat)
es_varmuts.m

es_varmuts_l.m <- emmeans(lm_var_muts, pairwise ~ locisigma.cat | COA.cat)
es_varmuts_l.m


es_varmuts_pr.m <- emmeans(lm_var_muts, pairwise ~ pleiorate.cat | COA.cat)
es_varmuts_pr.m

emm_varmuts_contr_pr.m <- pairs(pairs(emmeans(lm_var_muts, ~ pleiorate.cat | COA.cat,
                                          at = list(pleiorate.cat = c("Low", "Medium", "High"),
                                                    COA.cat = c("Gaussian", "House-of-Cards")))),
                            by = NULL)

emm_varmuts_contr_pr.m


es_varmuts_r.m <- emmeans(lm_var_muts, pairwise ~ rwide.cat | COA.cat)
es_varmuts_r.m




es_varmuts_pc.m <- emmeans(lm_var_muts, pairwise ~ pleiocov.cat | COA.cat)
es_varmuts_pc.m

emm_varmuts_contr_pc.m <- pairs(pairs(emmeans(lm_var_muts, ~ pleiocov.cat | COA.cat,
                                          at = list(pleiocov.cat = c("Low", "High"),
                                                    COA.cat = c("Gaussian", "House-of-Cards")))),
                            by = NULL)

emm_varmuts_contr_pc.m

################################################################################################################

# Kurtosis

lm_kurt_muts <- lm_robust(kurt ~ rwide.cat + locisigma.cat + COA.cat + pleiorate.cat + pleiocov.cat +
                            COA.cat*rwide.cat + COA.cat*locisigma.cat + COA.cat*pleiorate.cat + COA.cat*pleiocov.cat,
                          data = d_muts_stats[d_muts_stats$COA.cat != "Other" & d_muts_stats$COA.cat != "Null" & d_muts_stats$Po == "Adapted",],
                          se_type = "HC3") # Not looking at the in-between models for now

summary(lm_kurt_muts)

lm_kurt_muts_OLS <- lm(kurt ~ rwide.cat + locisigma.cat + COA.cat + pleiorate.cat + pleiocov.cat +
                         COA.cat*rwide.cat + COA.cat*locisigma.cat + COA.cat*pleiorate.cat + COA.cat*pleiocov.cat,
                       data = d_muts_stats[d_muts_stats$COA.cat != "Other" & d_muts_stats$COA.cat != "Null" & d_muts_stats$Po == "Adapted",]) # Not looking at the in-between models for now


calc.relimp(lm_kurt_muts_OLS)


es_kurtmuts.m <- emmeans(lm_kurt_muts, pairwise ~  COA.cat)
es_kurtmuts.m

es_kurt_l.m <- emmeans(lm_kurt_muts, pairwise ~ locisigma.cat | COA.cat)
es_kurt_l.m

emm_kurt_contr_l.m <- pairs(pairs(emmeans(lm_kurt_muts, ~ locisigma.cat | COA.cat,
                                         at = list(locisigma.cat = c("Low variance", "High variance"),
                                                   COA.cat = c("Gaussian", "House-of-Cards")))),
                           by = NULL)

emm_kurt_contr_l.m

es_kurt_pr.m <- emmeans(lm_kurt_muts, pairwise ~ pleiorate.cat | COA.cat)
es_kurt_pr.m

emm_kurt_contr_pr.m <- pairs(pairs(emmeans(lm_kurt_muts, ~ pleiorate.cat | COA.cat,
                                          at = list(pleiorate.cat = c("Low", "High"),
                                                    COA.cat = c("Gaussian", "House-of-Cards")))),
                            by = NULL)

emm_kurt_contr_pr.m


es_kurt_r.m <- emmeans(lm_kurt_muts, pairwise ~ rwide.cat | COA.cat)
es_kurt_r.m

emm_kurt_contr_r.m <- pairs(pairs(emmeans(lm_kurt_muts, ~ rwide.cat | COA.cat,
                                         at = list(rwide.cat = c("Low", "High"),
                                                   COA.cat = c("Gaussian", "House-of-Cards")))),
                           by = NULL)

emm_kurt_contr_r.m


es_kurt_pc.m <- emmeans(lm_kurt_muts, pairwise ~ pleiocov.cat | COA.cat)
es_kurt_pc.m

emm_kurt_contr_pc.m <- pairs(pairs(emmeans(lm_kurt_muts, ~ pleiocov.cat | COA.cat,
                                          at = list(pleiocov.cat = c("Low", "High"),
                                                    COA.cat = c("Gaussian", "House-of-Cards")))),
                            by = NULL)

emm_kurt_contr_pc.m

#####################################################################################################################

# Mutation counts

lm_count_muts <- lm_robust(count ~ rwide.cat + locisigma.cat + COA.cat + pleiorate.cat + pleiocov.cat +
                             COA.cat*rwide.cat + COA.cat*locisigma.cat + COA.cat*pleiorate.cat + COA.cat*pleiocov.cat,
                           data = d_muts_stats[d_muts_stats$COA.cat != "Other" & d_muts_stats$COA.cat != "Null" & d_muts_stats$Po == "Adapted",],
                           se_type = "HC3") # Not looking at the in-between models for now
summary(lm_count_muts)

lm_count_muts_OLS <- lm(count ~ rwide.cat + locisigma.cat + COA.cat + pleiorate.cat + pleiocov.cat +
                         COA.cat*rwide.cat + COA.cat*locisigma.cat + COA.cat*pleiorate.cat + COA.cat*pleiocov.cat,
                       data = d_muts_stats[d_muts_stats$COA.cat != "Other" & d_muts_stats$COA.cat != "Null" & d_muts_stats$Po == "Adapted",]) # Not looking at the in-between models for now


calc.relimp(lm_count_muts_OLS)



es_countmuts.m <- emmeans(lm_count_muts, pairwise ~  COA.cat)
es_countmuts.m

es_countmuts_l.m <- emmeans(lm_count_muts, pairwise ~ locisigma.cat | COA.cat)
es_countmuts_l.m

emm_countmuts_contr_l.m <- pairs(pairs(emmeans(lm_count_muts, ~ locisigma.cat | COA.cat,
                                          at = list(locisigma.cat = c("Low meaniance", "High meaniance"),
                                                    COA.cat = c("Gaussian", "House-of-Cards")))),
                            by = NULL)

emm_countmuts_contr_l.m

es_countmuts_pr.m <- emmeans(lm_count_muts, pairwise ~ pleiorate.cat | COA.cat)
es_countmuts_pr.m

emm_countmuts_contr_pr.m <- pairs(pairs(emmeans(lm_count_muts, ~ pleiorate.cat | COA.cat,
                                           at = list(pleiorate.cat = c("Low", "High"),
                                                     COA.cat = c("Gaussian", "House-of-Cards")))),
                             by = NULL)

emm_countmuts_contr_pr.m


es_countmuts_r.m <- emmeans(lm_count_muts, pairwise ~ rwide.cat | COA.cat)
es_countmuts_r.m

emm_countmuts_contr_r.m <- pairs(pairs(emmeans(lm_count_muts, ~ rwide.cat | COA.cat,
                                          at = list(rwide.cat = c("Low", "High"),
                                                    COA.cat = c("Gaussian", "House-of-Cards")))),
                            by = NULL)

emm_countmuts_contr_r.m


es_countmuts_pc.m <- emmeans(lm_count_muts, pairwise ~ pleiocov.cat | COA.cat)
es_countmuts_pc.m

emm_countmuts_contr_pc.m <- pairs(pairs(emmeans(lm_count_muts, ~ pleiocov.cat | COA.cat,
                                           at = list(pleiocov.cat = c("Low", "High"),
                                                     COA.cat = c("Gaussian", "House-of-Cards")))),
                             by = NULL)

emm_countmuts_contr_pc.m

##########################################################################################################

# For Table 2: need all effects of each parameter on dist, var, cov

# Distance

dplot_combined_dist_ls <- d_combined[d_combined$atopt  == "Adapted", c(4, 15, 17)] %>%
  group_by(COA.cat, locisigma.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))


dplot_combined_dist_r <- d_combined[d_combined$atopt  == "Adapted", c(4, 16, 17)] %>%
  group_by(COA.cat, rwide.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))



dplot_combined_dist_pr <- d_combined[d_combined$atopt  == "Adapted", c(4, 11, 17)] %>%
  group_by(COA.cat, pleiorate.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))


dplot_combined_dist_pc <- d_combined[d_combined$atopt  == "Adapted", c(4, 12, 17)] %>%
  group_by(COA.cat, pleiocov.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

###################################

# Var  

# Exclude outliers
dplot_combined_var_ls <- d_combined[d_combined$atopt  == "Adapted"  & d_combined$varmean < 50, c(21, 15, 17)] %>%
  group_by(COA.cat, locisigma.cat) %>%
  summarise_all(list(varmean_mean = mean, varmean_se = std.error, n = length))

dplot_combined_var_r <- d_combined[d_combined$atopt  == "Adapted", c(21, 16, 17)] %>%
  group_by(COA.cat, rwide.cat) %>%
  summarise_all(list(varmean_mean = mean, varmean_se = std.error, n = length))


dplot_combined_var_pr <- d_combined[d_combined$atopt  == "Adapted", c(21, 11, 17)] %>%
  group_by(COA.cat, pleiorate.cat) %>%
  summarise_all(list(varmean_mean = mean, varmean_se = std.error, n = length))


dplot_combined_var_pc <- d_combined[d_combined$atopt  == "Adapted", c(21, 12, 17)] %>%
  group_by(COA.cat, pleiocov.cat) %>%
  summarise_all(list(varmean_mean = mean, varmean_se = std.error, n = length))


###############################################


# cov
# Exclude outliers
dplot_combined_cov_ls <- d_combined[d_combined$atopt  == "Adapted" & d_combined$covmean > -1.2, c(22, 15, 17)] %>%
  group_by(COA.cat, locisigma.cat) %>%
  summarise_all(list(covmean_mean = mean, covmean_se = std.error, n = length))

dplot_combined_cov_r <- d_combined[d_combined$atopt  == "Adapted", c(22, 16, 17)] %>%
  group_by(COA.cat, rwide.cat) %>%
  summarise_all(list(covmean_mean = mean, covmean_se = std.error, n = length))

dplot_combined_cov_pr <- d_combined[d_combined$atopt  == "Adapted", c(22, 11, 17)] %>%
  group_by(COA.cat, pleiorate.cat) %>%
  summarise_all(list(covmean_mean = mean, covmean_se = std.error, n = length))

dplot_combined_cov_pc <- d_combined[d_combined$atopt  == "Adapted", c(22, 12, 17)] %>%
  group_by(COA.cat, pleiocov.cat) %>%
  summarise_all(list(covmean_mean = mean, covmean_se = std.error, n = length))


#############################################################################################################################################

# For Table 3: need a table of effects of each parameter on var, kurt, count

dmeans_muts_ls <- d_muts_stats[d_muts_stats$Po  == "Adapted", c(15, 13, 18:20)] %>%
  group_by(COA.cat, locisigma.cat) %>%
  summarise_all(list(mean = mean, se = std.error, n = length))

dmeans_muts_ls

dmeans_muts_r <- d_muts_stats[d_muts_stats$Po  == "Adapted", c(15, 10, 18:20)] %>%
  group_by(COA.cat, rwide.cat) %>%
  summarise_all(list(mean = mean, se = std.error, n = length))

dmeans_muts_r

dmeans_muts_pr <- d_muts_stats[d_muts_stats$Po  == "Adapted", c(15, 11, 18:20)] %>%
  group_by(COA.cat, pleiorate.cat) %>%
  summarise_all(list(mean = mean, se = std.error, n = length))

dmeans_muts_pr


dmeans_muts_pc <- d_muts_stats[d_muts_stats$Po  == "Adapted", c(15, 12, 18:20)] %>%
  group_by(COA.cat, pleiocov.cat) %>%
  summarise_all(list(mean = mean, se = std.error, n = length))

dmeans_muts_pr

