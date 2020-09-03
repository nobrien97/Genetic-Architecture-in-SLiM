# Analysis of Model 1 data: Null 

# Big files: data.table::use fread()


d_null <- data.table::fread("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM1/Output/out_8T_null_means.csv", header = F, integer64="character")

names(d_null)[1:6] <- c("gen", "seed", "modelindex", "rsd", "rwide", "delmu")
# d_null$seed <- as.factor(d_null$seed)

names(d_null)[7:34] <-  c(paste0("pleiocov_0", 1:7), paste0("pleiocov_1", 2:7), paste0("pleiocov_2", 3:7), paste0("pleiocov_3", 4:7), paste0("pleiocov_4", 5:7), paste0("pleiocov_5", 6:7), paste0("pleiocov_6", 7))

names(d_null)[35:42] <- paste0("mean", 0:7)

names(d_null)[43:50] <- paste0("var", 0:7)

names(d_null)[51:78] <- c(paste0("phenocov_0", 1:7), paste0("phenocov_1", 2:7), paste0("phenocov_2", 3:7), paste0("phenocov_3", 4:7), paste0("phenocov_4", 5:7), paste0("phenocov_5", 6:7), paste0("phenocov_6", 7))

names(d_null)[79:106] <- c(paste0("phenocor_0", 1:7), paste0("phenocor_1", 2:7), paste0("phenocor_2", 3:7), paste0("phenocor_3", 4:7), paste0("phenocor_4", 5:7), paste0("phenocor_5", 6:7), paste0("phenocor_6", 7))

names(d_null)[107] <- "H"


nrow(unique(d_null[d_null$gen == 150000,][,c('seed', 'modelindex')]))


d_null_sbst <- d_null[1:500000,]
d_null_sbst$seed <- as.factor(d_null_sbst$seed)

d_null_burn <- read.csv("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM1/Output/out_8T_null_burnin.csv", header = F)
names(d_null_burn) <- c("gen", "seed", "modelindex", "H")
