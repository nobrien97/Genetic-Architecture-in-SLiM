# Generate combinations of variables and save as a csv

library(DoE.wrapper)


# Generate hypercubes until our maximum pairwise correlation is less than 0.1
repeat{
  lhc_seed <- sample(1:2147483647, 1)
  lhc_samples <- lhs.design(
    nruns=512,
    nfactors=8,
    seed=lhc_seed,
    type="maximin",
    factor.names=list(
      Xstart = c(0.1, 5), 
      Xstop = c(5.1, 9), 
      Aalpha = c(0.00001, 10), 
      Abeta = c(0.00001, 10), 
      Balpha = c(0.00001, 10), 
      Bbeta = c(0.00001, 10),
      Hilln = c(10, 200), 
      Bthreshold = c(0.00001, 10))
  )
  if (max(abs(cor(lhc_samples)[cor(lhc_samples) != 1])) < 0.05)
    break
}


# Check correlations and space between samples

plot(lhc_samples)
cor(lhc_samples)
max(abs(cor(lhc_samples)[cor(lhc_samples) != 1]))
hist(lhc_samples$Xstart)

# Using GGally to plot correlation matrix
library(tidyverse)
library(GGally)
ggpairs(lhc_samples)

# Add on column for CPP or R: 1 = CPP, 0 = R

cpp_levels <- c(1, 0)
lhc_samples <-  lhc_samples[rep(seq_len(nrow(lhc_samples)), 2), ]
lhc_samples$cpp <- rep(cpp_levels, each = nrow(lhc_samples) / 2)

write.csv(lhc_samples, "./lhc_cppAUC.csv")
