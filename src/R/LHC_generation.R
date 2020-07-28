#################################################################
#               Latin Hypercube Generation Script               #
#################################################################

# This script generates the LHC samples for my final runs, and checks their quality.
# I also justify the ranges I'm using for each variable


# I have 9 factors:

# Population size (Ne): Ranges from 100 to 10,000: very small populations to intermediate sizes. 
# Has the largest performance impact. 

# Genome-wide recombination rate (rwide): Ranges from 0.0 to 1.907e-3: none to high (in plants) recombination
# Average recombination is 1.346e-5
# Nachman et al. 2002 - Drosophila may go above 5cM
# Stapley et al. 2017 - max in plants was 9.22 cM, leading to this max figure
# 
# This is a mean value used in an rnorm() function that gets recombination rates for each chromosome

# Genome-wide recombination rate standard deviation (rsd): Ranges from Standard deviation used in the rnorm() 
# function to get recombination rates for chromosomes

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
# No deleterious mutations to equal chance to a neutral mutation

# Pleiotropic mutational correlation (pleio_cov) - this defines the magnitude and direction of 
# correlations between traits in pleiotropic mutations, used in the multivariate normal pull to get 
# effects for all traits: Ranges from -0.5 to 0.5: 
# Ranges from strong negative to strong positive mutational correlations, 0 represents pure pleiotropy

# Number of chromosomes to have deleterious mutations: Ranges from 0 to 10, from none of them to all of them


# Generate latin hypercube sample for all the factors

library(DoE.wrapper)

lscombos <- lhs.design(
  nruns = 100, # 100 is too few, but this drastically increases the amount of runs: nruns * nseeds is the total. Also generating the LHC takes forever with more of these
  nfactors = 9,
  seed = 1370976025, #sampled from 1:2147483647
  #  type = "improved",
  factor.names = list(
    rregions = c(0.0, 0.5),
    rwide = c(0.0, 1.346e-5),  # Try it, don't know how long it will run!
    Ne = c(10, 10000),
    pleiocov = c(-0.5, 0.5), # More even sampling, going from -0.5 to 0.5 rather than 0 to 0.5
    pleiorate = c(0.0, 0.5),
    delmu = c(0.0, 1.0))
  #    npercent = c(0.00001, 0.8)) # Based on genome size of 100000 - Sinnott-Armstrong et al. 2020 estimate ~80,000 sites affecting human height, so we sample from 1 to 80000 QTLs
  
)
