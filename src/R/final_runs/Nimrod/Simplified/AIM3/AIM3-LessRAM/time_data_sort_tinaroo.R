
# Over time: need to get pop means, mean var and mean cov so file shouldn't be too big
# Will need to do with selection model

# For null, will have to calculate the optimum after the fact
# phenomeans at generation 50,000
# optimum = phenomeans+(phenomeans*(3.1746) # mu*nloci*50,000 = 6.3e-6*100*50,000 = 31.5
# 

# Test here first with a small section, do big null file on supercomputer, and selection model
set.seed(873662137)
setwd("/90days/s4395747")
source("src_G_mat.R")
library(plyr)
library(tidyverse)

# Sel model first, it's smaller


# Selection model: don't need to do anything to optima

# Linux version

d_sel <- data.table::fread("/90days/s4395747/d_sel_time.csv", header = T, integer64="character")
d_sel <- d_sel[, -1]
d_sel$seed <- as.numeric(d_sel$seed)

d_sel <- arrange(d_sel, gen, modelindex, seed)



# sel distances
d_sel_opt <- data.table::fread("/90days/s4395747/out_8T_stabsel_opt_c.csv", header = F, integer64="character")
names(d_sel_opt) <- c("seed", "modelindex", paste0("opt", 0:7))
d_sel_opt$seed <- as.numeric(d_sel_opt$seed)
d_sel_opt <- arrange(d_sel_opt, modelindex, seed)

# Run on supercomputer: extract means, and optima, do distance

ls_sel_dist <- MCeuc_dist_new(d_sel, d_sel_opt, 24)

d_sel_eucdist <- data.frame(
  gen = rep(unique(d_sel$gen), each = length(unique(d_sel$seed))*length(unique(d_sel$modelindex))),
  seed = sort(rep(unique(d_sel$seed), each = length(unique(d_sel$modelindex)))),
  modelindex = sort(unique(d_sel$modelindex)),
  distance = unlist(ls_sel_dist)
)


# Add parameters

ls_combos_sel <- read.csv("/90days/s4395747/lscombos_sel.csv")

d_sel_eucdist$delmu <- rep(ls_combos_sel$delmu, times = 20100)
d_sel_eucdist$rwide <- rep(ls_combos_sel$rwide, times = 20100)
d_sel_eucdist$locisigma <- rep(ls_combos_sel$locisigma, times = 20100)
d_sel_eucdist$pleiorate <- rep(ls_combos_sel$pleiorate, times = 20100)
d_sel_eucdist$pleiocov <- rep(ls_combos_sel$pleiocov, times = 20100)
d_sel_eucdist$tau <- rep(ls_combos_sel$tau, times = 20100)

saveRDS(d_sel_eucdist, "d_sel_eucdist.RDS")



d_null <- data.table::fread("/90days/s4395747/out_8T_null_means_c_2.csv", header = F, integer64="character")


names(d_null)[1:6] <- c("gen", "seed", "modelindex", "rsd", "rwide", "delmu")
# d_null$seed <- as.factor(d_null$seed)

names(d_null)[7:34] <-  c(paste0("pleiocov_0", 1:7), paste0("pleiocov_1", 2:7), paste0("pleiocov_2", 3:7), paste0("pleiocov_3", 4:7), paste0("pleiocov_4", 5:7), paste0("pleiocov_5", 6:7), paste0("pleiocov_6", 7))

names(d_null)[35:42] <- paste0("mean", 0:7)

names(d_null)[43:50] <- paste0("var", 0:7)

names(d_null)[51:78] <- c(paste0("phenocov_0", 1:7), paste0("phenocov_1", 2:7), paste0("phenocov_2", 3:7), paste0("phenocov_3", 4:7), paste0("phenocov_4", 5:7), paste0("phenocov_5", 6:7), paste0("phenocov_6", 7))

names(d_null)[79:106] <- c(paste0("phenocor_0", 1:7), paste0("phenocor_1", 2:7), paste0("phenocor_2", 3:7), paste0("phenocor_3", 4:7), paste0("phenocor_4", 5:7), paste0("phenocor_5", 6:7), paste0("phenocor_6", 7))

names(d_null)[107] <- "H"

d_null$seed <- as.numeric(d_null$seed)


d_null_opt <- d_null[d_null$gen == 50000]
d_null_opt <- cbind(d_null_opt[, c(2:3)], (d_null_opt[, c(35:42)] + (d_null_opt[, c(35:42)]*(100/31.5))))
names(d_null_opt) <- c("seed", "modelindex", paste0("opt", 0:7))


# For testing, double check results make sense with rowMeans and that we are taking the right columns 
d_null_test <- d_null[d_null$seed == 1140555014 & d_null$modelindex == 1 & d_null$gen == 150000,]

d_null$varmean <- rowMeans(d_null[, c(43:50)])
d_null$covmean <- rowMeans(d_null[, c(51:78)])

d_null <- d_null[, -c(7:34, 43:107)]
write.csv(d_null_opt, "d_null_opt.csv", row.names = F)
write.csv(d_null, "d_null_time.csv", row.names = F)
saveRDS(d_null, "d_null_time.csv")


d_null <- data.table::fread("/90days/s4395747/d_null_time.csv", header = T, integer64="character")
d_null_opt <- data.table::fread("/90days/s4395747/d_null_opt.csv", header = T, integer64="character")

# Add an extra variable now so the function takes the correct values (the means)

d_null$tau <- 0.0
d_null <- d_null[,c(1:6, 17, 7:16)]

nrow(d_null) 
# Split the file up to calculate distances separately - otherwise too much RAM used
# 20582601 rows total
# Could arrange by modelindex and seed etc. then create two lists of 512 models
# Then do the same for opt, keeping 512 models


d_null$seed <- as.numeric(d_null$seed)
d_null_opt$seed <- as.numeric(d_null_opt$seed)

d_null <- arrange(d_null, gen, modelindex, seed)
d_null_opt <- arrange(d_null_opt, modelindex, seed)

d_null_pt1 <- d_null[d_null$modelindex < 129,]

d_null_pt8 <- d_null[d_null$modelindex > 128 & d_null$modelindex < 257,]

d_null_pt2 <- d_null[d_null$modelindex > 256 & d_null$modelindex < 385,]

d_null_pt3 <- d_null[d_null$modelindex > 384 & d_null$modelindex < 513,]

d_null_pt4 <- d_null[d_null$modelindex > 512 & d_null$modelindex < 641,]

d_null_pt5 <- d_null[d_null$modelindex > 640 & d_null$modelindex < 769,]

d_null_pt6 <- d_null[d_null$modelindex > 768 & d_null$modelindex < 897,]

d_null_pt7 <- d_null[d_null$modelindex > 896 & d_null$modelindex < 1025,]



d_null_opt_pt1 <- d_null_opt[d_null_opt$modelindex < 129,]

d_null_opt_pt8 <- d_null_opt[d_null_opt$modelindex > 128 & d_null_opt$modelindex < 257,]

d_null_opt_pt2 <- d_null_opt[d_null_opt$modelindex > 256 & d_null_opt$modelindex < 385,]

d_null_opt_pt3 <- d_null_opt[d_null_opt$modelindex > 384 & d_null_opt$modelindex < 513,]

d_null_opt_pt4 <- d_null_opt[d_null_opt$modelindex > 512 & d_null_opt$modelindex < 641,]

d_null_opt_pt5 <- d_null_opt[d_null_opt$modelindex > 640 & d_null_opt$modelindex < 769,]

d_null_opt_pt6 <- d_null_opt[d_null_opt$modelindex > 768 & d_null_opt$modelindex < 897,]

d_null_opt_pt7 <- d_null_opt[d_null_opt$modelindex > 896 & d_null_opt$modelindex < 1025,]


write.csv(d_null_pt1, "d_null_pt1.csv", row.names = F)
write.csv(d_null_pt8, "d_null_pt8.csv", row.names = F)
write.csv(d_null_pt2, "d_null_pt2.csv", row.names = F)
write.csv(d_null_pt3, "d_null_pt3.csv", row.names = F)
write.csv(d_null_pt4, "d_null_pt4.csv", row.names = F)
write.csv(d_null_pt5, "d_null_pt5.csv", row.names = F)
write.csv(d_null_pt6, "d_null_pt6.csv", row.names = F)
write.csv(d_null_pt7, "d_null_pt7.csv", row.names = F)

write.csv(d_null_opt_pt1, "d_null_opt_pt1.csv", row.names = F)
write.csv(d_null_opt_pt8, "d_null_opt_pt8.csv", row.names = F)
write.csv(d_null_opt_pt2, "d_null_opt_pt2.csv", row.names = F)
write.csv(d_null_opt_pt3, "d_null_opt_pt3.csv", row.names = F)
write.csv(d_null_opt_pt4, "d_null_opt_pt4.csv", row.names = F)
write.csv(d_null_opt_pt5, "d_null_opt_pt5.csv", row.names = F)
write.csv(d_null_opt_pt6, "d_null_opt_pt6.csv", row.names = F)
write.csv(d_null_opt_pt7, "d_null_opt_pt7.csv", row.names = F)

rm(list=ls())
source("src_G_mat.R")


d_null_pt1 <- data.table::fread("/90days/s4395747/d_null_pt1.csv", header = T, integer64="character")
d_null_opt_pt1 <- data.table::fread("/90days/s4395747/d_null_opt_pt1.csv", header = T, integer64="character")
d_null_pt1$seed <- as.numeric(d_null_pt1$seed)
d_null_opt_pt1$seed <- as.numeric(d_null_opt_pt1$seed)

ls_null_dist_1 <- MCeuc_dist_new(d_null_pt1, d_null_opt_pt1, 24)
d_null_eucdist_pt1 <- data.frame(
  gen = rep(unique(d_null_pt1$gen), each = length(unique(d_null_pt1$seed))*length(unique(d_null_pt1$modelindex))),
  seed = rep(unique(d_null_pt1$seed), each = length(unique(d_null_pt1$modelindex))),
  modelindex = sort(unique(d_null_pt1$modelindex)),
  distance = unlist(ls_null_dist_1)
)
saveRDS(d_null_eucdist_pt1, "d_null_eucdist_pt1.RDS")


d_null_pt8 <- data.table::fread("/90days/s4395747/d_null_pt8.csv", header = T, integer64="character")
d_null_opt_pt8 <- data.table::fread("/90days/s4395747/d_null_opt_pt8.csv", header = T, integer64="character")
d_null_pt8$seed <- as.numeric(d_null_pt8$seed)
d_null_opt_pt8$seed <- as.numeric(d_null_opt_pt8$seed)

ls_null_dist_8 <- MCeuc_dist_new(d_null_pt8, d_null_opt_pt8, 24)
d_null_eucdist_pt8 <- data.frame(
  gen = rep(unique(d_null_pt8$gen), each = length(unique(d_null_pt8$seed))*length(unique(d_null_pt8$modelindex))),
  seed = rep(unique(d_null_pt8$seed), each = length(unique(d_null_pt8$modelindex))),
  modelindex = sort(unique(d_null_pt8$modelindex)),
  distance = unlist(ls_null_dist_8)
)
saveRDS(d_null_eucdist_pt8, "d_null_eucdist_pt8.RDS")



d_null_pt2 <- data.table::fread("/90days/s4395747/d_null_pt2.csv", header = T, integer64="character")
d_null_opt_pt2 <- data.table::fread("/90days/s4395747/d_null_opt_pt2.csv", header = T, integer64="character")
d_null_pt2$seed <- as.numeric(d_null_pt2$seed)
d_null_opt_pt2$seed <- as.numeric(d_null_opt_pt2$seed)

ls_null_dist_2 <- MCeuc_dist_new(d_null_pt2, d_null_opt_pt2, 24)
d_null_eucdist_pt2 <- data.frame(
  gen = rep(unique(d_null_pt2$gen), each = length(unique(d_null_pt2$seed))*length(unique(d_null_pt2$modelindex))),
  seed = rep(unique(d_null_pt2$seed), each = length(unique(d_null_pt2$modelindex))),
  modelindex = sort(unique(d_null_pt2$modelindex)),
  distance = unlist(ls_null_dist_2)
)
saveRDS(d_null_eucdist_pt2, "d_null_eucdist_pt2.RDS")


d_null_pt3 <- data.table::fread("/90days/s4395747/d_null_pt3.csv", header = T, integer64="character")
d_null_dist_3 <- d_null_dist_3 %>% distinct()
d_null_opt_pt3 <- data.table::fread("/90days/s4395747/d_null_opt_pt3.csv", header = T, integer64="character")
d_null_pt3$seed <- as.numeric(d_null_pt3$seed)
d_null_opt_pt3$seed <- as.numeric(d_null_opt_pt3$seed)

ls_null_dist_3 <- MCeuc_dist_new(d_null_pt3, d_null_opt_pt3, 24)
d_null_eucdist_pt3 <- data.frame(
  gen = rep(unique(d_null_pt3$gen), each = length(unique(d_null_pt3$seed))*length(unique(d_null_pt3$modelindex))),
  seed = rep(unique(d_null_pt3$seed), each = length(unique(d_null_pt3$modelindex))),
  modelindex = sort(unique(d_null_pt3$modelindex)),
  distance = unlist(ls_null_dist_3)
)
saveRDS(d_null_eucdist_pt3, "d_null_eucdist_pt3.RDS")


d_null_pt4 <- data.table::fread("/90days/s4395747/d_null_pt4.csv", header = T, integer64="character")
d_null_opt_pt4 <- data.table::fread("/90days/s4395747/d_null_opt_pt4.csv", header = T, integer64="character")
d_null_pt4$seed <- as.numeric(d_null_pt4$seed)
d_null_opt_pt4$seed <- as.numeric(d_null_opt_pt4$seed)

ls_null_dist_4 <- MCeuc_dist_new(d_null_pt4, d_null_opt_pt4, 24)
d_null_eucdist_pt4 <- data.frame(
  gen = rep(unique(d_null_pt4$gen), each = length(unique(d_null_pt4$seed))*length(unique(d_null_pt4$modelindex))),
  seed = rep(unique(d_null_pt4$seed), each = length(unique(d_null_pt4$modelindex))),
  modelindex = sort(unique(d_null_pt4$modelindex)),
  distance = unlist(ls_null_dist_4)
)
saveRDS(d_null_eucdist_pt4, "d_null_eucdist_pt4.RDS")


d_null_pt5 <- data.table::fread("/90days/s4395747/d_null_pt5.csv", header = T, integer64="character")
d_null_opt_pt5 <- data.table::fread("/90days/s4395747/d_null_opt_pt5.csv", header = T, integer64="character")
d_null_pt5$seed <- as.numeric(d_null_pt5$seed)
d_null_opt_pt5$seed <- as.numeric(d_null_opt_pt5$seed)

ls_null_dist_5 <- MCeuc_dist_new(d_null_pt5, d_null_opt_pt5, 24)
d_null_eucdist_pt5 <- data.frame(
  gen = rep(unique(d_null_pt5$gen), each = length(unique(d_null_pt5$seed))*length(unique(d_null_pt5$modelindex))),
  seed = rep(unique(d_null_pt5$seed), each = length(unique(d_null_pt5$modelindex))),
  modelindex = sort(unique(d_null_pt5$modelindex)),
  distance = unlist(ls_null_dist_5)
)
saveRDS(d_null_eucdist_pt5, "d_null_eucdist_pt5.RDS")


d_null_pt6 <- data.table::fread("/90days/s4395747/d_null_pt6.csv", header = T, integer64="character")
d_null_opt_pt6 <- data.table::fread("/90days/s4395747/d_null_opt_pt6.csv", header = T, integer64="character")
d_null_pt6$seed <- as.numeric(d_null_pt6$seed)
d_null_opt_pt6$seed <- as.numeric(d_null_opt_pt6$seed)

ls_null_dist_6 <- MCeuc_dist_new(d_null_pt6, d_null_opt_pt6, 24)
d_null_eucdist_pt6 <- data.frame(
  gen = rep(unique(d_null_pt6$gen), each = length(unique(d_null_pt6$seed))*length(unique(d_null_pt6$modelindex))),
  seed = rep(unique(d_null_pt6$seed), each = length(unique(d_null_pt6$modelindex))),
  modelindex = sort(unique(d_null_pt6$modelindex)),
  distance = unlist(ls_null_dist_6)
)
saveRDS(d_null_eucdist_pt6, "d_null_eucdist_pt6.RDS")



d_null_pt7 <- data.table::fread("/90days/s4395747/d_null_pt7.csv", header = T, integer64="character")
d_null_opt_pt7 <- data.table::fread("/90days/s4395747/d_null_opt_pt7.csv", header = T, integer64="character")
d_null_pt7$seed <- as.numeric(d_null_pt7$seed)
d_null_opt_pt7$seed <- as.numeric(d_null_opt_pt7$seed)

ls_null_dist_7 <- MCeuc_dist_new(d_null_pt7, d_null_opt_pt7, 24)
d_null_eucdist_pt7 <- data.frame(
  gen = rep(unique(d_null_pt7$gen), each = length(unique(d_null_pt7$seed))*length(unique(d_null_pt7$modelindex))),
  seed = rep(unique(d_null_pt7$seed), each = length(unique(d_null_pt7$modelindex))),
  modelindex = sort(unique(d_null_pt7$modelindex)),
  distance = unlist(ls_null_dist_7)
)
saveRDS(d_null_eucdist_pt7, "d_null_eucdist_pt7.RDS")


d_null_eucdist_pt1 <- readRDS("d_null_eucdist_pt1.RDS")
d_null_eucdist_pt2 <- readRDS("d_null_eucdist_pt8.RDS") # number 8 is actually number 2, due to RAM problems and needing to split stuff up even more
d_null_eucdist_pt3 <- readRDS("d_null_eucdist_pt2.RDS")
d_null_eucdist_pt4 <- readRDS("d_null_eucdist_pt3.RDS")
d_null_eucdist_pt5 <- readRDS("d_null_eucdist_pt4.RDS")
d_null_eucdist_pt6 <- readRDS("d_null_eucdist_pt5.RDS")
d_null_eucdist_pt7 <- readRDS("d_null_eucdist_pt6.RDS")
d_null_eucdist_pt8 <- readRDS("d_null_eucdist_pt7.RDS")


d_null_eucdist <- data.table::rbindlist(list(d_null_eucdist_pt1, d_null_eucdist_pt2, d_null_eucdist_pt3, d_null_eucdist_pt4, d_null_eucdist_pt5, d_null_eucdist_pt6, d_null_eucdist_pt7, d_null_eucdist_pt8))

# Add parameters 

d_null_eucdist <- arrange(d_null_eucdist, gen, modelindex, seed)

ls_combos_null <- read.csv("/90days/s4395747/lscombos_null.csv")

d_null_eucdist$delmu <- rep(ls_combos_null$delmu, each = 100, times = 201)
d_null_eucdist$rwide <- rep(ls_combos_null$rwide, each = 100, times = 201)
d_null_eucdist$pleiorate <- rep(ls_combos_null$pleiorate, each = 100, times = 201)
d_null_eucdist$pleiocov <- rep(ls_combos_null$pleiocov, each = 100, times = 201)
d_null_eucdist$locisigma <- rep(ls_combos_null$locisigma, each = 100, times = 201)
d_null_eucdist$tau <- 0.0

saveRDS(d_null_eucdist, "d_null_eucdist.RDS")

# Combine the two

d_sel_eucdist <- readRDS("d_sel_eucdist.RDS")

d_sel_eucdist <- arrange(d_sel_eucdist, gen, modelindex, seed)

d_sel_eucdist$modelindex <- d_sel_eucdist$modelindex + 1024

d_eucdist_c <- data.table::rbindlist(list(d_null_eucdist, d_sel_eucdist), use.names=T)
d_eucdist_c <- arrange(d_eucdist_c, gen, modelindex, seed)
d_eucdist_c$delmu.cat <- cut(d_eucdist_c$delmu, breaks = 3, labels = c("Low", "Medium", "High"))
d_eucdist_c$locisigma.cat <- cut(d_eucdist_c$locisigma, breaks = 3, labels = c("Low", "Medium", "High")) 
d_eucdist_c$pleiorate.cat <- cut(d_eucdist_c$pleiorate, breaks = 3, labels = c("Low", "Medium", "High")) 
d_eucdist_c$pleiocov.cat <- cut(d_eucdist_c$pleiocov, breaks = 3, labels = c("Low", "Medium", "High")) 
d_eucdist_c$rwide.cat <- cut(d_eucdist_c$rwide, breaks = 3, labels = c("Low", "Medium", "High")) 

# Custom breaks for Tau so we can differentiate from null models and selection models
tau_bp <- c(-Inf, 0.1, 333.3, 666.6, Inf)

d_eucdist_c$tau.cat <- cut(d_eucdist_c$tau, breaks = tau_bp, labels = c("Null", "High", "Medium", "Low")) 

saveRDS(d_eucdist_c, "d_eucdist_c.RDS")

d_eucdist_c <- readRDS("d_eucdist_c.RDS")

# Delmu * rwide * ls * tau
d_mean_eucdist <- d_eucdist_c[, -c(2:3, 5:10)] %>%
  group_by(gen, delmu.cat, rwide.cat, locisigma.cat, pleiorate.cat, pleiocov.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error, dist_var = var))

d_mean_eucdist_var <- d_eucdist_c[, -c(1:3, 5:10)] %>%
  group_by(delmu.cat, rwide.cat, locisigma.cat, pleiorate.cat, pleiocov.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error, dist_var = var))


# d_mean_eucdist_notau for null and sel models
d_eucdist_c_null <- d_eucdist_c[d_eucdist_c$tau.cat == "Null",]
d_eucdist_c_sel <- d_eucdist_c[d_eucdist_c$tau.cat != "Null",]

# Delmu * rwide * ls
d_mean_eucdist_notau_null <- d_eucdist_c_null[, -c(2:3, 5:10, 13)] %>%
  group_by(gen, delmu.cat, rwide.cat, locisigma.cat, pleiorate.cat, pleiocov.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error, dist_var = var))

d_mean_eucdist_notau_sel <- d_eucdist_c_sel[, -c(2:3, 5:10, 13)] %>%
  group_by(gen, delmu.cat, rwide.cat, locisigma.cat, pleiorate.cat, pleiocov.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error, dist_var = var))

d_mean_eucdist_null_var <- d_eucdist_c_null[, -c(1:3, 5:10, 13)] %>%
  group_by(delmu.cat, rwide.cat, locisigma.cat, pleiorate.cat, pleiocov.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error, dist_var = var))

d_mean_eucdist_sel_var <- d_eucdist_c_sel[, -c(1:3, 5:10, 13)] %>%
  group_by(delmu.cat, rwide.cat, locisigma.cat, pleiorate.cat, pleiocov.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error, dist_var = var))



saveRDS(d_mean_eucdist, "d_mean_eucdist.RDS")
saveRDS(d_mean_eucdist_notau_null, "d_mean_eucdist_notau_null.RDS")
saveRDS(d_mean_eucdist_notau_sel, "d_mean_eucdist_notau_sel.RDS")

saveRDS(d_mean_eucdist_var, "d_mean_eucdist_var.RDS")
saveRDS(d_mean_eucdist_null_var, "d_mean_eucdist_null_var.RDS")
saveRDS(d_mean_eucdist_sel_var, "d_mean_eucdist_sel_var.RDS")



################################################################################################
################################################################################################
# Variance/covariance: import the null_time thing from before, combine with sel and export RDS #
################################################################################################
################################################################################################

# Remove previous stuff
rm(list = ls())

set.seed(873662137)
setwd("/90days/s4395747")
source("src_G_mat.R")
library(plyr)
library(tidyverse)


d_sel <- data.table::fread("/90days/s4395747/d_sel_time.csv", header = T, integer64="character")
d_sel <- d_sel[, -1]
d_sel$seed <- as.numeric(d_sel$seed)

d_sel <- arrange(d_sel, gen, modelindex, seed)
d_sel$modelindex <- d_sel$modelindex + 1024



d_null <- data.table::fread("/90days/s4395747/d_null_time.csv", header = T, integer64="character")
d_null <- d_null %>% distinct() # Get rid of any duplicate rows
d_null$seed <- as.numeric(d_null$seed)
d_null <- arrange(d_null, gen, modelindex, seed)

d_null[d_null$gen == 50000 & d_null$modelindex == 2 & d_null$seed == 829606697,]


ls_combos_null <- read.csv("/90days/s4395747/lscombos_null.csv")

d_null$delmu <- rep(ls_combos_null$delmu, each = 100, times = 201)
d_null$rwide <- rep(ls_combos_null$rwide, each = 100, times = 201)
d_null$pleiorate <- rep(ls_combos_null$pleiorate, each = 100, times = 201)
d_null$pleiocov <- rep(ls_combos_null$pleiocov, each = 100, times = 201)
d_null$locisigma <- rep(ls_combos_null$locisigma, each = 100, times = 201)
d_null$tau <- 0.0



ls_combos_sel <- read.csv("/90days/s4395747/lscombos_sel.csv")

d_sel$delmu <- rep(ls_combos_sel$delmu, each = 100, times = 201)
d_sel$rwide <- rep(ls_combos_sel$rwide, each = 100, times = 201)
d_sel$locisigma <- rep(ls_combos_sel$locisigma, each = 100, times = 201)
d_sel$pleiorate <- rep(ls_combos_sel$pleiorate, each = 100, times = 201)
d_sel$pleiocov <- rep(ls_combos_sel$pleiocov, each = 100, times = 201)
d_sel$tau <- rep(ls_combos_sel$tau, each = 100, times = 201)


# Combine
d_raw_c <- data.table::rbindlist(list(d_null, d_sel), use.names=T)


# Add categories, save RDS

d_raw_c <- arrange(d_raw_c, gen, modelindex, seed)
d_raw_c$delmu.cat <- cut(d_raw_c$delmu, breaks = 3, labels = c("Low", "Medium", "High"))
d_raw_c$locisigma.cat <- cut(d_raw_c$locisigma, breaks = 3, labels = c("Low", "Medium", "High")) 
d_raw_c$pleiorate.cat <- cut(d_raw_c$pleiorate, breaks = 3, labels = c("Low", "Medium", "High")) 
d_raw_c$pleiocov.cat <- cut(d_raw_c$pleiocov, breaks = 3, labels = c("Low", "Medium", "High")) 
d_raw_c$rwide.cat <- cut(d_raw_c$rwide, breaks = 3, labels = c("Low", "Medium", "High")) 

# Custom breaks for Tau so we can differentiate from null models and selection models
tau_bp <- c(-Inf, 0.1, 333.3, 666.6, Inf)

d_raw_c$tau.cat <- cut(d_raw_c$tau, breaks = tau_bp, labels = c("Null", "High", "Medium", "Low")) 

saveRDS(d_raw_c, "d_raw_c.RDS")

write.csv(d_raw_c, "d_raw_c.csv")

library(estimatr)
library(emmeans)

d_raw_c_fingen <- subset(d_raw_c, gen == 150000)

d_raw_c_fingen <- d_raw_c[d_raw_c$gen == 150000,]

var_fingen_lm <- lm_robust(varmean ~ delmu.cat * rwide.cat * locisigma.cat + 
			tau.cat * delmu.cat + tau.cat * locisigma.cat + tau.cat * rwide.cat, data = d_raw_c_fingen)

summary(var_fingen_lm)



# Delmu * rwide * ls * tau
d_mean_var <- d_raw_c[, c(1, 15:16, 21:26)] %>%
  group_by(gen, delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error, var = var))

d_var_null <- d_raw_c[d_raw_c$tau.cat == "Null"]
d_var_sel <- d_raw_c[d_raw_c$tau.cat != "Null"]

d_mean_var_null <- d_var_null[, c(1, 15:16, 21:25)] %>%
  group_by(gen, delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error, var = var))

d_mean_var_sel <- d_var_sel[, c(1, 15:16, 21:25)] %>%
  group_by(gen, delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error, var = var))

saveRDS(d_mean_var, "d_mean_var.RDS")
saveRDS(d_mean_var_null, "d_mean_var_null.RDS")
saveRDS(d_mean_var_sel, "d_mean_var_sel.RDS")

d_raw_end <- d_raw_c[d_raw_c$gen == 150000,]

saveRDS(d_raw_end, "d_raw_end.RDS")

# Variance of variance and covariance over time: cycling around the optimum

d_mean_var_var <- d_raw_c[, c(15:16, 21:26)] %>%
  group_by(delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error, var = var))

d_mean_var_var_null <- d_var_null[, c(15:16, 21:25)] %>%
  group_by(delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error, var = var))

d_mean_var_var_sel <- d_var_sel[, c(15:16, 21:25)] %>%
  group_by(delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error, var = var))

saveRDS(d_mean_var_var, "d_mean_var_var.RDS")
saveRDS(d_mean_var_var_null, "d_mean_var_var_null.RDS")
saveRDS(d_mean_var_var_sel, "d_mean_var_var_sel.RDS")



