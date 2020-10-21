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

d_null_muts <- d_sel_muts[d_null_muts$type != 2,]


d_sel_muts <- arrange(d_sel_muts, seed, modelindex)
d_null_muts <- arrange(d_null_muts, seed, modelindex)

saveRDS(d_sel_muts, "d_sel_muts.RDS")                        
saveRDS(d_null_muts, "d_null_muts.RDS")                        

d_sel_muts <- readRDS("d_sel_muts.RDS")
d_null_muts <- readRDS("d_null_muts.RDS")


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


ls_combos_sel <- read.csv("/90days/s4395747/lscombos_sel.csv")

d_sel_long$delmu <- ls_combos_sel$delmu[d_sel_long$modelindex]
d_sel_long$rwide <- ls_combos_sel$rwide[d_sel_long$modelindex]
d_sel_long$pleiorate <- ls_combos_sel$pleiorate[d_sel_long$modelindex]
d_sel_long$pleiocov <- ls_combos_sel$pleiocov[d_sel_long$modelindex]
d_sel_long$locisigma <- ls_combos_sel$locisigma[d_sel_long$modelindex]
d_sel_long$tau <- ls_combos_sel$tau[d_sel_long$modelindex]


plot_sel_freq <- ggplot()





saveRDS(d_sel_muts, "d_sel_muts.RDS")                        
