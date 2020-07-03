##############################################################################################################
#  Simple test script for including two traits, with arguments coming from latin square, for use on cluster
##############################################################################################################


#  Parallel script modified from SLiM-Extras example R script, info at
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.


# Change this!
setwd("D:/uni/Honours/Programming/SLiM/2Traits_LatinSquare")

# Parallelisation libraries 

library(foreach)
library(doParallel)
library(future)


# Seed generation - 100 seeds from the total range of possible values in a 32 bit signed integer

rsample <- as.character(runif(100, 1, (2^62 - 1)))

write.table(rbind(rsample), file = "seeds.csv", row.names = FALSE, col.names = FALSE, sep=",")

# Generate latin hypercube sample for all the factors

library(DoE.wrapper)

lscombos <- lhs.design(
  nruns = 100, # 100 is too few, but this drastically increases the amount of runs: nruns * nseeds is the total. Also generating the LHC takes forever with more of these
  nfactors = 6,
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

# Time test

 system.time(lhs.design(
  nruns = 150, # 100 is too few, but this drastically increases the amount of runs: nruns * nseeds is the total. Also generating the LHC takes forever with more of these
  nfactors = 7,
  seed = 1370976025, #sampled from 1:2147483647
  #  type = "improved",
  factor.names = list(
    rregions = c(0.0, 0.5),
    rwide = c(0.0, 1e-6), # Try it, don't know how long it will run!
    Ne = c(10, 10000),
    pleiocov = c(-0.5, 0.5), # More even sampling, going from -0.5 to 0.5 rather than 0 to 0.5
    pleiorate = c(0.0, 1.0),
    delmu = c(0.0, 1.0),
    npercent = c(0.00001, 0.8))) 
  
)

 

# Write to disk

write.csv(lscombos, file = "lscombos.csv", row.names = FALSE)



# Diagnostic plot for sampling

hist(lscombos$rregions)

plot(lscombos)
cor(lscombos)

# Parameter setup - import Latin Square data frame
lscombos <- read.csv("lscombos150.csv", header = T)


# - run on all available cores, treat them as nodes

cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

#Run SLiM through the WSL, defining parameters in command line

foreach(i=rsample) %:%
  foreach(j=1:nrow(lscombos)) %dopar% {
    # Use string manipulation functions to configure the command line args, feeding from a data frame of latin square combinations
    # then run SLiM with system(),
    slim_out <- system(sprintf("wsl slim -d seed=%s -d rregion=%s -d Ne=%i -d QTL_cov=%f -d pleiorate=%f -d delmu=%f -d rwide=%s -d -d modelindex=%i neutral2T10GRjune_restruct.slim", 
                               as.character(i), as.character(lscombos$rregions[j]), as.integer(round(lscombos$Ne[j])), lscombos$pleiocov[j], lscombos$pleiorate[j], lscombos$delmu[j], as.character(lscombos$rwide[j]), j, intern=T))
  }
stopCluster(cl)

# Time test

system.time(system("wsl slim -d seed=123 -d rregion=0.5 -d Ne=10000 -d QTL_cov=0.5 -d pleiorate=1.0 -d delmu=0.0 -d rwide=1.346e-5 neutral2T10GRjune_restruct2_opt_timetest.slim"))
# Run gc() after to find the maximum memory usage
gc()


# 100 loci
# With matrix output, wsl slim -d seed=123 -d rregion=0.5 -d Ne=10000 -d QTL_cov=0.5 -d pleiorate=1.0 -d delmu=0.0 -d rwide=1.346e-5
#user  system elapsed 
#0.02    0.00  492.15

# Without matrix output, wsl slim -d seed=123 -d rregion=0.5 -d Ne=10000 -d QTL_cov=0.5 -d pleiorate=1.0 -d delmu=0.0 -d rwide=1.346e-5
#user  system elapsed 
#0.08    0.00  133.95 


#133.95*100 = 13395   + (492.15 - 133.95)  = 13753.2/run  * 10000 = 137532000/60/60 = 38203.33 hrs; 




# This time with only 10% as many linkage groups - 10% of the loci are selected as recombination end points

#With matrix output
#user  system elapsed 
#0.03    0.00  315.72 

#Without matrix output
#user  system elapsed 
#0.08    0.05   47.66 

# 4766s for running the sim, + (315.72 - 47.66) = 268.06 for matrix output calc, 5034.06
# 5034.06*10000 runs = 50340600s or 13983.5 hours / 192 =72.83 hours
# add 10% for run just in case: 80.113hrs
# /576 = 26.7 hrs
# /384 = 40.06 hrs



# Estimate for Job array: 24 cores, 1 node, doing 500 operations (5 lscombos, 100 seeds)
# 13983.5 / 20 = 699.175 cpu hrs for 500 runs
# 699.175/24 = 29.132 hrs on 24 cores, add 20% to be safe: ~35 hrs




gens = c(1000, 2500, 5000, 7500, 10000)
runtime = c(24.87, 57.47, 351.77, 402.86, 703.15)
genstime = data.frame(gens, runtime)

abline(mC <- lm(runtime ~ gens, data = genstime))
plot(gens, runtime)


# Generating matrices is SLOW, as shown in above plot
# Scales linearly with pop size
# It is super fast if there are no mutations or substitutions though (it doesn't have to go through each genome)




#    user  system elapsed 
# 24.347   0.081  ~26 
# Result on Awoonga for 1000 generations at rregion 0.5. Ne 10000 QTL_cov 0.5 pleiorate 1.0 delmu 1.0 rwide 1e-06
# For a full run multiply this worst case scenario to get to total generations e.g. 24.625*50 to get to 50000 generations, 1231.25s

# For 100,000 generations, this is 2600s/run, ~43 minutes
# *100 seeds and 150 runs, gets as to 10833 hours cputime
# For 96 cores (4 nodes), 112.85 hours
# For 192 cores (8 nodes), 56.42 hours
# Add on top of this the time to generate matrices (above)



# 100000 generations next time
# More seeds, need to treat seeds as a nested random effect, as it is correlated with population size

# Everything seems to scale pretty linearly: rregion definitely has the largest impact, considering how 
# high the recombination rate can be (0.5): this interacts with population size as well, so some runs will be
# very slow, but this will be linearly offset by the runs which are super fast, which from the histograms of
# lscombos are equally frequent, so the average is a pretty good indicator of the total time it will take I 
# think.


# Data import

slim_out <- read.csv("out2T10L_lhstest.csv", header = F)

names(slim_out) <- c("gen", "seed", "rregions", "rwide", "delmu", "pleiocov", "Ne", "pheno0mean", "pheno0var", "pheno1mean", "pheno1var", "phenocov", "phenocor")

slim_out$seed <- as.factor(slim_out$seed)

# Think about which interactions to think about - rather than a full factorial of all effects + all interactions
# Genome recombination rate: maybe 
# need another column, model - for each parameter combination, add another column with model number (e.g. one combo might have model 1, next will have model 2 and so on - pretty much opposite of seed)

# There is more evolution in larger populations because the allelic effects linger in the population for longer without being lost by drift: heterozygosity stays for longer
# seeds represent a repeat of evolution - different initial conditions (maybe a different ecotype, or environment etc.)
# Space that the total data is sampling given by trait value vs consecutive number (row number)
#   Use this to determine the rare runs that deviate from the null (0) and can increase the sampling on these runs to do stats!
# Population of models is the number of seeds - enables us 
# Unlocking data sets by calibrating populations of models to data density (our models are populations of models, data density is data density)
# experimental data cloud must be one model from our population of models!

# These models summarise a population of models that can be compared to experiment
# computational methods in biology





meanstest <- read.csv("out_2T100L_means.csv", header = F)

names(meanstest) <- c("gen", "seed", "modelindex", "rregions", "rwide", "delmu", "pleiocov", "Ne", "ploci", "pheno0mean", "pheno0var", "pheno1mean", "pheno1var", "phenocov", "phenocor")

library(ggplot2)

seedchoice <- sample(unique(as.factor(meanstest$seed)), 1)

ggplot(meanstest,
       aes(x = gen, y = phenocov, colour = as.factor(seed))) +
  geom_line() +
#  gghighlight(seed == as.factor(seedchoice[1]) | seed == as.factor(seedchoice[2])) +
  scale_colour_manual(values=c("#CC6666", "#66CC99")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Generation", y = "pheno0 mean")



# PBS array test run - make sure my ranges are right for 20 nodes, 5 lscombos per node

for (i in ((PBSIND*5) - 4):(PBSIND*5)) {
  print(i)
}



# Tinaroo speed test

# Gens = 1000   system.time(system("/home/$USER/SLiM/slim -d seed=123 -d rregion=0.5 -d Ne=10000 -d QTL_cov=0.5 -d pleiorate=0.5 -d delmu=0.0 -d rwide=1.346e-5 /home/$USER/SLiM/Scripts/null2T100L_june_restruct2_opt_nomat.slim")) 
# user  system elapsed 
# 56.104   0.549  57.871 

# With matrix:
# user  system elapsed 
# 116.102   5.917 123.196 

# 57.871*99 + 123.196 = 5852.425s per run
# 16256.73611 hrs for 10000 runs
# /20 for cpu time per node = 812.8368
# /24 for actual time per node = 33.8682hrs
# Add 20% just in case: 40.64 hrs



