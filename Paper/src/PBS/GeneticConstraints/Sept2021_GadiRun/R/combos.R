# Generate combinations of variables and save as a csv

library(DoE.wrapper)


# Generate hypercubes until our maximum pairwise correlation is less than 0.1
repeat{
    lhc_seed <- sample(1:2147483647, 1)
    lhc_samples <- lhs.design(
        nruns=64,
        nfactors=6,
        seed=lhc_seed,
        type="maximin",
        factor.names=list(
            rwide=c(0.0,0.5),
            locisigma=c(0.1,10),
            Ne=c(100,10000),
            nloci=c(1,100),
            width=c(0.0002,0.5),
            optShift=c(1,100))
  )
  if (max(abs(cor(lhc_samples)[cor(lhc_samples) != 1])) < 0.1)
    break
}


# Check correlations and space between samples

plot(lhc_samples)
cor(lhc_samples)
max(abs(cor(lhc_samples)[cor(lhc_samples) != 1]))
hist(lhc_samples$locisigma)

# Using GGally to plot correlation matrix
library(tidyverse)
library(GGally)
ggpairs(lhc_samples)


# Add on column for constraint - will need to repeat the
# sampling 3 times for K = Low, Medium, High

K_levels <- c('"c(1.0, 0.0, 0.0)"', '"c(0.0, 1.0, 0.0)"', '"c(0.0, 0.0, 1.0)"')
lhc_samples <-  lhc_samples[rep(seq_len(nrow(lhc_samples)), 3), ]
lhc_samples$K <- rep(K_levels, each = nrow(lhc_samples) / 3)

write.csv(lhc_samples, "./lhc_gaditest.csv")

# Generate cmds.txt

singleRunBashName <- "./constraints_testSR.sh"
modelIndices <- seq_len(nrow(lhc_samples))
# Test run - go for 10 hours, maximum of 2976 cores, 2880 is 15*192 models, so that will do
seeds <- 1:15#length(c(read.csv("seeds.csv"))$Seed)

dfCMDs <- data.frame(name = singleRunBashName,
                     modelindex = rep(modelIndices, times = length(seeds)),
                     seed = rep(seeds, each = length(modelIndices))
                     )

write.table(dfCMDs, "./cmds.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)




# TODO: Have more of a look into this - I don't think it's necessary for purely reading files
# Might be worth combining model outputs at the end though instead of writing to one file
# Split data frame into separate files for each row - 
# to avoid filesystem problems from reading the same file on multiple threads
if (!dir.exists("./combos"))
    dir.create("./combos")

for (row in seq.int(1, nrow(lhc_samples))) {
    write.csv(lhc_samples[row,], paste0("./combos/lhc_sample", row, ".csv"), row.names = F)
}

# Do the same for seeds:

seeds <- read.csv("./May2021-GadiTest/R/seeds_gaditest.csv")

if (!dir.exists("./seeds"))
    dir.create("./seeds")

for (row in seq.int(1, nrow(seeds))) {
    write.table(seeds[row,], paste0("./seeds/seed", row, ".csv"), row.names = F, col.names = F, sep = ",")
}

