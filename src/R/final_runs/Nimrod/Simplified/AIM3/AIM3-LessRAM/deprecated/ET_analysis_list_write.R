
setwd("/90days/s4395747")

source("src_G_mat.R")
library(tidyverse)
library(evolqg)
library(parallel)
library(rrapply)
library(tidyr)

d_raw_mat <- readRDS("d_raw_mat.RDS")

# Sample 1000 times to get a nice big matrix of eigenanalysis of eigentensor 1
# Will take around 20 minutes to do this step on Tinaroo
ls_ET <- eigentensor_G(d_raw_mat, 1000, 24)

# Organise into a more reasonable output for transforming to data frame

ls_ET_org <- MCOrg_ET(ls_ET, 24)

ETGvecs <- c(paste0("Gmax.vec", 1:8), paste0("G2.vec", 1:8))

# Generate data frame of eigenvalues and vectors
d_ETG <- ListToDF_ET(ls_ET_org, ETGvecs)

# Remove unnecessary third column (contains the eigentensor: we could do multiple eigentensors, and then we would keep this)
d_ETG <- d_ETG[,-3]

names(d_ETG)[1:2] <- c("Replicate", "delmu.rwide")

d_ETG$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number


# These values are stored as list objects, make them a regular numeric vector
d_ETG$Gmax.val <- as.numeric(d_ETG$Gmax.val)
d_ETG$G2.val <- as.numeric(d_ETG$G2.val)

d_ETG <- separate(data = d_ETG,
                  col = delmu.rwide,
                  into = c("delmu", "rwide"))

# Reorder the factor levels to low - medium - high
d_ETG$delmu <- factor(d_ETG$delmu, levels = c("Low", "Medium", "High"))
d_ETG$rwide <- factor(d_ETG$rwide, levels = c("Low", "Medium", "High"))

write.csv(d_ETG, "d_ET.csv", row.names = F)
