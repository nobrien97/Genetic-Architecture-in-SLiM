
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
d_muts_stats <- read.csv("d_muts_stats.csv")

########################

# Fix up some factor labels d_eucdist_c and d_raw_c

levels(d_eucdist_c$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")
levels(d_raw_c$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")
levels(d_eucdist_c$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_raw_c$tau.cat) <- c("Null", "Strong", "Medium", "Weak")

d_eucdist_c$locisigma.cat <- factor(d_eucdist_c$locisigma.cat, levels = c("Low variance", "Medium variance", "High variance"))
d_raw_c$locisigma.cat <- factor(d_raw_c$locisigma.cat, levels = c("Low variance", "Medium variance", "High variance"))
d_combined$locisigma.cat <- factor(d_combined$locisigma.cat, levels = c("Low variance", "Medium variance", "High variance"))
d_muts_stats$locisigma.cat <- factor(d_muts_stats$locisigma.cat, levels = c("Low variance", "Medium variance", "High variance"))

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



# Table 2: Contingency table/Chi sq

library(chisq.posthoc.test)

conttab_dist <- table(d_combined[d_combined$COA.cat != "Other",]$atopt, d_combined[d_combined$COA.cat != "Other",]$COA.cat)
chisq.test(conttab_dist)

chisq.posthoc.test(conttab_dist)


# Variance linear model

d_combined_stat <- d_combined
d_combined_stat$rwide <- scale(d_combined_stat$rwide)
d_combined_stat$locisigma <- scale(d_combined_stat$locisigma)
d_combined_stat$pleiorate <- scale(d_combined_stat$pleiorate)
d_combined_stat$pleiocov <- scale(d_combined_stat$pleiocov)



d_combined_stat$COA.cat <- factor(d_combined_stat$COA.cat, levels = c("Null", "Gaussian", "House-of-Cards", "Other"))


lm_var_end <- lm(varmean ~ rwide + locisigma + COA.cat + pleiorate + pleiocov +
                   COA.cat*rwide + COA.cat*locisigma + COA.cat*pleiorate + COA.cat*pleiocov,
                 data = d_combined_stat[d_combined_stat$COA.cat != "Other" & d_combined$atopt == "Adapted",]) # Not looking at the in-between models for now

summary(lm_var_end)


lm_cov_end <- lm(covmean ~ rwide + locisigma + COA.cat + pleiorate + pleiocov +
                   COA.cat*rwide + COA.cat*locisigma + COA.cat*pleiorate + COA.cat*pleiocov,
                 data = d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",]) # Not looking at the in-between models for now

summary(lm_cov_end)

lm_dist_end <- lm(distance ~ rwide + locisigma + COA.cat + pleiorate + pleiocov +
                    COA.cat*rwide + COA.cat*locisigma + COA.cat*pleiorate + COA.cat*pleiocov,
                  data = d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",]) # Not looking at the in-between models for now

summary(lm_dist_end)



# Use lm_robust to adjust for non-uniform errors, non-normality accounted for by sample size
library(estimatr)
library(xtable)

lm_cov_end <- lm_robust(covmean ~ (rwide + locisigma + COA.cat + pleiorate)^2,
                        data = d_raw_end[d_raw_end$COA.cat != "Other",])


