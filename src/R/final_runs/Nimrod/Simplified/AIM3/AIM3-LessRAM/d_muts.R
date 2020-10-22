set.seed(873662137)
setwd("/90days/s4395747")
source("src_G_mat.R")
library(plyr)
library(tidyverse)

d_null_muts <- read.csv("/30days/s4395747/null_muts_fixed.csv", header = F)

d_sel_muts <- read.csv("/30days/s4395747/sel_muts_fixed.csv", header = F)

names(d_sel_muts) <- c("seed", "modelindex", "type", "pos", "t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7", "freq")

names(d_null_muts) <- c("seed", "modelindex", "type", "pos", "t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7", "freq")


# Remove m2 deleterious mutations
d_sel_muts <- d_sel_muts[d_sel_muts$type != 2,]

d_null_muts <- d_null_muts[d_null_muts$type != 2,]

# Get rid of NA rows
d_sel_muts <- na.omit(d_sel_muts)

d_null_muts <- na.omit(d_null_muts)

# Arrange by seed and muts
d_sel_muts <- arrange(d_sel_muts, seed, modelindex)
d_null_muts <- arrange(d_null_muts, seed, modelindex)

saveRDS(d_sel_muts, "d_sel_muts.RDS")
saveRDS(d_null_muts, "d_null_muts.RDS")

d_sel_muts <- readRDS("d_sel_muts.RDS")
d_null_muts <- readRDS("d_null_muts.RDS")

# R decided to make this a factor for some reason?
d_null_muts$t2 <- as.numeric(d_null_muts$t2)
d_null_muts$t7 <- as.numeric(d_null_muts$t7)

# Takes 90.8GB in R as a data.frame!
# Reduce size by only looking at 50 seeds

# We have the seed set so we know we will sample the same numbers each time

# temporarily convert seed to character so we don't get any floating point funny business
d_null_muts$seed <- as.character(d_null_muts$seed)
seed_sample <- sample(unique(d_null_muts$seed), 50)
d_null_muts <- d_null_muts[d_null_muts$seed %in% seed_sample,]
# Convert back
d_null_muts$seed <- as.numeric(d_null_muts$seed)

saveRDS(d_null_muts, "d_null_muts_sbst.RDS")

d_sel_muts <- d_sel_muts[d_sel_muts$seed == seed_sample,]



# Size effect vs freq graph: combine columns

d_sel_long <- d_sel_muts %>% 
  pivot_longer(
    cols = c(paste0("t", 0:7)),
    names_to = "Trait",
    values_to = "Effect"
  )

d_null_long <- d_null_muts %>% 
  pivot_longer(
    cols = c(paste0("t", 0:7)),
    names_to = "Trait",
    values_to = "Effect"
  )

d_sel_long <- d_sel_long[d_sel_long$Effect != 0.0,]

d_null_long <- d_null_long[d_null_long$Effect != 0.0,]
# Could also use Effect < 0.00000000001 & Effect > -0.00000000001, but R does this anyway with tolerances

ls_combos_sel <- read.csv("/90days/s4395747/lscombos_sel.csv")

d_sel_long$delmu <- ls_combos_sel$delmu[d_sel_long$modelindex]
d_sel_long$rwide <- ls_combos_sel$rwide[d_sel_long$modelindex]
d_sel_long$pleiorate <- ls_combos_sel$pleiorate[d_sel_long$modelindex]
d_sel_long$pleiocov <- ls_combos_sel$pleiocov[d_sel_long$modelindex]
d_sel_long$locisigma <- ls_combos_sel$locisigma[d_sel_long$modelindex]
d_sel_long$tau <- ls_combos_sel$tau[d_sel_long$modelindex]

saveRDS(d_sel_long, "d_sel_long.RDS")


ls_combos_null <- read.csv("/90days/s4395747/lscombos_null.csv")

d_null_long$delmu <- ls_combos_null$delmu[d_null_long$modelindex]
d_null_long$rwide <- ls_combos_null$rwide[d_null_long$modelindex]
d_null_long$pleiorate <- ls_combos_null$pleiorate[d_null_long$modelindex]
d_null_long$pleiocov <- ls_combos_null$pleiocov[d_null_long$modelindex]
d_null_long$locisigma <- ls_combos_null$locisigma[d_null_long$modelindex]
d_null_long$tau <- 0.0

saveRDS(d_null_long, "d_null_long.RDS")

d_null_long <- readRDS("d_null_long.RDS")

d_null_long <- d_null_long[d_null


d_sel_long <- readRDS("d_sel_long.RDS")


d_sel_long$modelindex <- d_sel_long$modelindex + 1024


# Combine the data sets:
d_muts_c <- data.table::rbindlist(list(d_null_long, d_sel_long), use.names=T)







