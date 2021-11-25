# Generate combinations of variables and save as a csv

library(DoE.wrapper)


# Generate hypercubes until our maximum pairwise correlation is less than 0.1
repeat{
    lhc_seed <- sample(1:2147483647, 1)
    lhc_samples <- lhs.design(
        nruns=32,
        nfactors=2,
        seed=lhc_seed,
        type="maximin",
        factor.names=list(
            rwide=c(0.0,0.5),
            Ne=c(100,10000))
  )
  if (max(abs(cor(lhc_samples)[cor(lhc_samples) != 1])) < 0.05)
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

# Add on column for tree or control: 1 = tree, 0 = control

tree_levels <- c(1, 0)
lhc_samples <-  lhc_samples[rep(seq_len(nrow(lhc_samples)), 2), ]
lhc_samples$tree <- rep(tree_levels, each = nrow(lhc_samples) / 2)

write.csv(lhc_samples, "./lhc_coalburntest.csv")


# Generate cmds.txt
lhc_samples <- read.csv("./lhc_coalburntest.csv")
singleRunBashName <- "./CoalBurnTest2_testSR.sh"
modelIndices <- seq_len(32)#nrow(lhc_samples))
seeds <- 1:10

dfCMDs <- data.frame(name = singleRunBashName,
                     modelindex = rep(modelIndices, times = length(seeds)),
                     seed = rep(seeds, each = length(modelIndices))
                     )

write.table(dfCMDs, "./cmds2.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
