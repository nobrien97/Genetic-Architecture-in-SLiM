# Script for generating the LHC and seeds used in my pilot run
# Generating 3 seeds and 5 combinations for a small data set to set my analysis up with


# LHC samples:

library(DoE.wrapper)

lscombos_null_pilot <- lhs.design(
  nruns = 8, # this drastically increases the amount of runs: nruns * nseeds is the total. Also generating the LHC takes forever with more of these
  nfactors = 8,
  seed = 1370976025, #sampled from 1:2147483647
  type = "optimum",
  factor.names = list(
    Ne = c(10, 5000), # Smaller population size for pilot run, just so we can get results faster
    rwide = c(0.0, 1.241e-4), 
    nloci = c(10, 500),
    locisigma = c(0.1, 10),
    pleiorate = c(0.0, 1.0),
    delmu = c(0.0, 1.0),
    pleiocov = c(-0.5, 0.5), # More even sampling, going from -0.5 to 0.5 rather than 0 to 0.5
    delchr = c(0, 10))
)

write.csv(lscombos_null_pilot, "lscombos_null_pilot.csv")


lscombos_sel_pilot <- lhs.design(
  nruns = 9, # this drastically increases the amount of runs: nruns * nseeds is the total. Also generating the LHC takes forever with more of these
  nfactors = 9,
  seed = 1370976025, #sampled from 1:2147483647
  type = "optimum",
  factor.names = list(
    Ne = c(10, 5000), # Smaller population size so we can get pilot results faster
    rwide = c(0.0, 1.241e-4), 
    nloci = c(10, 500),
    locisigma = c(0.1, 10),
    pleiorate = c(0.0, 1.0),
    delmu = c(0.0, 1.0),
    pleiocov = c(-0.5, 0.5), # More even sampling, going from -0.5 to 0.5 rather than 0 to 0.5
    delchr = c(0, 10),
    wsd = c(0.1, 10))
)

write.csv(lscombos_sel_pilot, "lscombos_sel_pilot.csv")


# Seeds 

rsample <- as.character(runif(3, 1, (2^62 - 1)))
write.table(c("Seed", rsample), file = "seeds_pilot.csv", row.names = FALSE, col.names = FALSE, sep=",")
