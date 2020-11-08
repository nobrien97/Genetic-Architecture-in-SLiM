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

# Strength of selection (tau): ranges from 10 to 1000 (very strong to very weak selection)
# fitness always relative to the optimum so it is comparable across seeds/models. 10 means that the sd 
# will be much smaller relative to the optimum, 1000 means much larger, so the curve is flatter (weaker selection)

# Generate latin hypercube sample for all the factors

library(DoE.wrapper)


# Null and recom model LHS Generation: Generate 512 samples in case we can do more, fewer factors 
# since we now fix nloci, Ne, and remove delchr
# Using maximin algorithm: not as good as optimum in minimising correlations, but much much faster 
# (won't be possible to use optimum for 512 samples)

lscombos_nul_maxi <- lhs.design(
  nruns = 1024, # this drastically increases the amount of runs: nruns * nseeds is the total. Also generating the LHC takes forever with more of these
  nfactors = 5,
  seed = 1370976025, #sampled from 1:2147483647
  type = "maximin",
  factor.names = list(
    rwide = c(0.0, 1.241e-4), 
    locisigma = c(0.1, 10),
    pleiorate = c(0.0, 0.5),
    delmu = c(0.0, 1.0),
    pleiocov = c(0.0, 0.5)) # going from 0 to 0.5 to avoid non-positive-definite errors
)

plot(lscombos_nul_maxi)
cor(lscombos_nul_maxi)
max(abs(cor(lscombos_nul_maxi)[cor(lscombos_nul_maxi) != 1]))
hist(lscombos_nul_maxi$locisigma)

write.csv(lscombos_nul_maxi, "lscombos_null.csv")

# Selection model LHS Generation

lscombos_sel_maxi <- lhs.design(
  nruns = 128, # this drastically increases the amount of runs: nruns * nseeds is the total. Also generating the LHC takes forever with more of these
  nfactors = 6,
  seed = 1782079225, #sampled from 1:2147483647
  type = "maximin",
  factor.names = list(
    rwide = c(0.0, 1.241e-4), 
    locisigma = c(0.1, 10),
    pleiorate = c(0.0, 0.5),
    delmu = c(0.0, 1.0),
    pleiocov = c(0.0, 0.5), # More even sampling, going from -0.5 to 0.5 rather than 0 to 0.5
    tau = c(10, 1000)) # From virtually no selection (very weak, weaker than drift) to very strong
)
plot(lscombos_sel_maxi)
hist(lscombos_sel_maxi$tau)

write.csv(lscombos_sel_maxi, "lscombos_sel.csv")

# Seed generation

rsample <- as.character(round(runif(100, 1, (2^32 - 1)))) # Pull from 32 bit integer range for now to reduce problems, below
# Seems to be a problem with the seeds themselves: 
# they are treated as the same if too large, even though they are below 64 bit limit
# They are above the 32 bit limit though, which is the strange thing: may be a bug?
# Mutation information and G matrices are still different among seeds, may be a problem with sample()
write.table(c("Seed", rsample), file = "seeds.csv", row.names = FALSE, col.names = FALSE, sep=",")

plot_rsample <- hist(as.numeric(rsample))

# Augmenting the selection one to 256 nruns since we have time

ls_combos_sel_aug <- lhs.augment(lscombos_sel_maxi, m = 128, type = "optAugment", seed = 2086426644) # seed sampled from sample(1:2147483647, 1)
write.csv(ls_combos_sel_aug, "lscombos_sel.csv")
