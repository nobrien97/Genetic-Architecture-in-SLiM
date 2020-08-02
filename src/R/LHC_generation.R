#################################################################
#               Latin Hypercube Generation Script               #
#################################################################

# This script generates the LHC samples for my final runs, and checks their quality.
# I also justify the ranges I'm using for each variable


# I have 9 factors:

# Population size (Ne): Ranges from 100 to 10,000: very small populations to intermediate sizes. 
# Has the largest performance impact. 

# Genome-wide recombination rate (rwide): Ranges from 0.0 to 1.241e-4: none to high (in plants) recombination
# Average recombination is 1.346e-5
# Nachman et al. 2002 - Drosophila may go above 5cM
# Stapley et al. 2017 - max in plants was 9.22 cM, leading to this max figure - 9.22e-8*1346 bp/locus = 1.241e-4
# 
# This is a mean value used in an rnorm() function that gets recombination rates for each chromosome

# Genome-wide recombination rate standard deviation (rsd): Ranges from Standard deviation used in the rnorm() 
# function to get recombination rates for chromosomes
# Fix this at 0.1 * rwide: can't go much higher on the max recombination, it'll get too slow, and recombination rates 
# vary too much to effectively sample with everything else to do in a timely manner

# Number of loci (nloci): Ranges from 10 to 500: very few loci affecting a trait to many 
# (from experimental results from QTL studies - McKown et al. 2014; Ehrenreich et al. 2010; Lu et al. 2010)

# Distribution of allelic effects (locisigma): Ranges from 0.1 to 10: very small effects to large effects
# Optimum will always be a fixed distance away from the burnt-in population, small and large is a relative term
# to how much closer they are getting to the optimum with a single effect

# Rate of pleiotropy (pleiorate): Ranges from 0.0 to 0.5: 
# No pleiotropy to a relatively realistic amount (Chesmore et al. 2017 found 0.44 in human disease traits)
# Stern and Orgogozo 2008: mentioned that in Ohya et al. 2005, 35% of deletion mutations affected two or more
# characters. Also, Dudley et al. 2005 found 58% reduced growth in two or more conditions

# Rate of deleterious mutation (delmu): Ranges from 0.0 to 1.0: 
# No deleterious mutations to equal chance to a neutral mutation - as expected with neutral theory

# Pleiotropic mutational correlation (pleio_cov) - this defines the magnitude and direction of 
# correlations between traits in pleiotropic mutations, used in the multivariate normal pull to get 
# effects for all traits: Ranges from -0.5 to 0.5: 
# Ranges from strong negative to strong positive mutational correlations, 0 represents pure pleiotropy

# Number of chromosomes to have deleterious mutations (delchr): Ranges from 0 to 10, from none of them to all of them

# Fitness function sd multiplier (wsd): ranges from 0.1 to 10 (very strong to very weak selection)
# fitness sd always relative to the optimum so it is comparable across seeds/models. 0.1 means that the sd 
# will be much smaller relative to the optimum, 10 means much larger, so the curve is flatter (weaker selection)

# Generate latin hypercube sample for all the factors

library(DoE.wrapper)

lscombos_imp <- lhs.design(
  nruns = 100, # this drastically increases the amount of runs: nruns * nseeds is the total. Also generating the LHC takes forever with more of these
  nfactors = 9,
  seed = 1370976025, #sampled from 1:2147483647
  type = "improved",
  factor.names = list(
    Ne = c(10, 10000),
    rwide = c(0.0, 1.241e-4), 
    nloci = c(10, 500),
    locisigma = c(0.1, 10),
    pleiorate = c(0.0, 1.0),
    delmu = c(0.0, 1.0),
    pleiocov = c(-0.5, 0.5), # More even sampling, going from -0.5 to 0.5 rather than 0 to 0.5
    delchr = c(0, 10),
    wsd = c(0.1, 10))
)


lscombos_maximin <- lhs.design(
  nruns = 100, # this drastically increases the amount of runs: nruns * nseeds is the total. Also generating the LHC takes forever with more of these
  nfactors = 9,
  seed = 1370976025, #sampled from 1:2147483647
  type = "maximin",
  factor.names = list(
    Ne = c(10, 10000),
    rwide = c(0.0, 1.241e-4), 
    nloci = c(10, 500),
    locisigma = c(0.1, 10),
    pleiorate = c(0.0, 1.0),
    delmu = c(0.0, 1.0),
    pleiocov = c(-0.5, 0.5), # More even sampling, going from -0.5 to 0.5 rather than 0 to 0.5
    delchr = c(0, 10),
    wsd = c(0.1, 10))
)

lscombos_opt <- lhs.design(
  nruns = 100, # this drastically increases the amount of runs: nruns * nseeds is the total. Also generating the LHC takes forever with more of these
  nfactors = 9,
  seed = 1370976025, #sampled from 1:2147483647
  type = "optimum",
  factor.names = list(
    Ne = c(10, 10000),
    rwide = c(0.0, 1.241e-4), 
    nloci = c(10, 500),
    locisigma = c(0.1, 10),
    pleiorate = c(0.0, 1.0),
    delmu = c(0.0, 1.0),
    pleiocov = c(-0.5, 0.5), # More even sampling, going from -0.5 to 0.5 rather than 0 to 0.5
    delchr = c(0, 10),
    wsd = c(0.1, 10))
)

lscombos_rand <- lhs.design(
  nruns = 100, # this drastically increases the amount of runs: nruns * nseeds is the total. Also generating the LHC takes forever with more of these
  nfactors = 9,
  seed = 1370976025, #sampled from 1:2147483647
  type = "random",
  factor.names = list(
    Ne = c(10, 10000),
    rwide = c(0.0, 1.241e-4), 
    nloci = c(10, 500),
    locisigma = c(0.1, 10),
    pleiorate = c(0.0, 1.0),
    delmu = c(0.0, 1.0),
    pleiocov = c(-0.5, 0.5), # More even sampling, going from -0.5 to 0.5 rather than 0 to 0.5
    delchr = c(0, 10),
    wsd = c(0.1, 10))
)

lscombos_ga <- lhs.design(
  nruns = 100, # this drastically increases the amount of runs: nruns * nseeds is the total. Also generating the LHC takes forever with more of these
  nfactors = 9,
  seed = 1370976025, #sampled from 1:2147483647
  type = "genetic",
  factor.names = list(
    Ne = c(10, 10000),
    rwide = c(0.0, 1.241e-4), 
    nloci = c(10, 500),
    locisigma = c(0.1, 10),
    pleiorate = c(0.0, 1.0),
    delmu = c(0.0, 1.0),
    pleiocov = c(-0.5, 0.5), # More even sampling, going from -0.5 to 0.5 rather than 0 to 0.5
    delchr = c(0, 10),
    wsd = c(0.1, 10))
)


plot(lscombos_ga)
cor(lscombos_ga)
max(abs(cor(lscombos_ga)[cor(lscombos_ga) != 1]))

plot(lscombos_imp)
cor(lscombos_imp)
max(abs(cor(lscombos_imp)[cor(lscombos_imp) != 1]))


plot(lscombos_maximin)
cor(lscombos_maximin)
max(abs(cor(lscombos_maximin)[cor(lscombos_maximin) != 1]))


plot(lscombos_opt)
cor(lscombos_opt)
max(abs(cor(lscombos_opt)[cor(lscombos_opt) != 1]))


plot(lscombos_rand)
cor(lscombos_rand)
max(abs(cor(lscombos_rand)[cor(lscombos_rand) != 1]))


cor(lscombos)
