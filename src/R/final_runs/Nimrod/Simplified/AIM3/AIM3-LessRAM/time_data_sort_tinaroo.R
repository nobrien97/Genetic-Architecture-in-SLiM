
# Over time: need to get pop means, mean var and mean cov so file shouldn't be too big
# Will need to do with selection model

# For null, will have to calculate the optimum after the fact
# phenomeans at generation 50,000
# optimum = phenomeans+(phenomeans*(3.1746) # mu*nloci*50,000 = 6.3e-6*100*50,000 = 31.5
# 

# Test here first with a small section, do big null file on supercomputer, and selection model

setwd("/90days/s4395747")
source("/90days/src_G_mat.R")
library(plyr)
library(tidyverse)

# Sel model first, it's smaller


# Selection model: don't need to do anything to optima

# Linux version

d_sel <- data.table::fread("/90days/d_sel_time.csv", header = F, integer64="character")
d_sel <- d_sel[, -1]
d_sel$seed <- as.numeric(d_sel$seed)



# sel distances
d_sel_opt <- data.table::fread("/90days/out_8T_stabsel_opt_c.csv", header = F, integer64="character")
names(d_sel_opt) <- c("seed", "modelindex", paste0("opt", 0:7))
d_sel_opt$seed <- as.numeric(d_sel_opt$seed)

# Run on supercomputer: extract means, and optima, do distance

ls_sel_means <- MCeuc_dist_new(d_sel, d_sel_opt, 24)

d_sel_eucdist <- data.frame(
  gen = rep(unique(d_sel$gen), each = length(unique(d_sel$seed))*length(unique(d_sel$modelindex))),
  seed = sort(rep(unique(d_sel$seed), each = length(unique(d_sel$modelindex)))),
  modelindex = sort(unique(d_sel$modelindex)),
  distance = unlist(ls_sel_means)
)

saveRDS(d_sel_eucdist, "d_sel_eucdist.RDS")

# Add parameters

ls_combos_sel <- read.csv("/90days/lscombos_sel.csv")

d_sel_eucdist$delmu <- rep(ls_combos_sel$delmu, times = 20100)
d_sel_eucdist$rwide <- rep(ls_combos_sel$rwide, times = 20100)
d_sel_eucdist$locisigma <- rep(ls_combos_sel$locisigma, times = 20100)
d_sel_eucdist$tau <- rep(ls_combos_sel$tau, times = 20100)



d_null <- data.table::fread("/90days/s4395747/out_8T_null_means_c2.csv", header = F, integer64="character")


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
d_null_opt <- cbind(d_null_opt[, c(2:3)], d_null_opt[, c(35:42)] + (d_null_opt[, c(35:42)]*(100/31.5)))
names(d_null_opt) <- c("seed", "modelindex", paste0("opt", 0:7))

d_null$varmean <- rowMeans(d_null[, c(43:50)])
d_null$covmean <- rowMeans(d_null[, c(51:78)])

d_null <- d_null[, -c(43:107)]
write.csv(d_null_opt, "d_null_opt.csv", row.names = F)
write.csv(d_null, "d_null_time.csv", row.names = F)


ls_null_means <- MCeuc_dist_new(d_null, d_null_opt, 24)

d_null_eucdist <- data.frame(
  gen = rep(unique(d_null$gen), each = length(unique(d_null$seed))*length(unique(d_null$modelindex))),
  seed = sort(rep(unique(d_null$seed), each = length(unique(d_null$modelindex)))),
  modelindex = sort(unique(d_null$modelindex)),
  distance = unlist(ls_null_means)
)

saveRDS(d_null_eucdist, "d_null_eucdist.RDS")

# Add parameters 

ls_combos_null <- read.csv("/90days/lscombos_sel.csv")

d_null_eucdist$delmu <- rep(ls_combos_null$delmu, times = 20100)
d_null_eucdist$rwide <- rep(ls_combos_null$rwide, times = 20100)
d_null_eucdist$locisigma <- rep(ls_combos_null$locisigma, times = 20100)
d_null_eucdist$tau <- 0.0


# Combine the two


d_sel_eucdist$modelindex <- d_sel$modelindex + 1024

d_eucdist_c <- bind_rows(d_nul_eucdist, d_sel_eucdist)
d_eucdist_c$delmu.cat <- cut(d_eucdist_c$delmu, breaks = 3, labels = c("Low", "Medium", "High"))
d_eucdist_c$locisigma.cat <- cut(d_eucdist_c$locisigma, breaks = 3, labels = c("Low", "Medium", "High")) 
d_eucdist_c$rwide.cat <- cut(d_eucdist_c$rwide, breaks = 3, labels = c("Low", "Medium", "High")) 

# Custom breaks for Tau so we can differentiate from null models and selection models
tau_bp <- c(-Inf, 0.1, 333.3, 666.6, Inf)

d_eucdist_c$tau.cat <- cut(d_eucdist_c$tau, breaks = tau_bp, labels = c("Null", "High", "Medium", "Low")) 

saveRDS(d_eucdist_c, "d_eucdist_c.RDS")


# Delmu * rwide * ls
d_meaneucdist <- d_raw_mat[, -c(2:3)] %>%
  group_by(gen, delmu.cat, rwide.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error, dist_var = var))
