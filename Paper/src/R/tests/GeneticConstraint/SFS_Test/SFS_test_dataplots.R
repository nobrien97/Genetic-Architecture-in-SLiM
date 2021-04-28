# In this script we load data from SFS_test and try to plot some SFSs, mutation frequencies, and allelelic
# effects, as well as trait effects

# For this run, we had Ne = 8000, rwide = 0.5, nloci = 30, genomelength = 40, 
# locisigma = c(1, 10)
# All other settings can be seen in the Paper/src/PBS/slim/polygen_maint.slim file

setwd("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/Paper/data/tests/GeneticConstraint/SFS_Test/")

source("../../../../src/R/includes/plot_function.R")


d_burnmeans <- read.csv("out_stabsel_burnin.csv", header = F)

names(d_burnmeans) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean")

d_means <- read.csv("out_stabsel_means.csv", header = F)

names(d_means) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", "dist", "mean_w", "delta")

d_means$sigma <- as.factor(d_means$modelindex)

levels(d_means$sigma) <- c("\u03c3 = 1", "\u03c3 = 10")

d_muts <- read.csv("out_stabsel_muts.csv", header = F) # Ignore extra empty column (fixed in next version of script)

names(d_muts) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", "constraint", "originGen", "value", "chi", "Freq", "fixGen")

# Why is one line missing a generation value? Because it stuck sim.generation with the previous fixGen
# Easy fix, but will have to think about the code to stop this happening in future
d_muts[9835,] <- c(78300, unlist(d_muts[9835, 1:11], use.names = F))
d_muts[9834,]$fixGen <- 76559

d_muts$constraint <- as.factor(d_muts$constraint)

levels(d_muts$constraint) <- c("Low", "Medium", "High")

d_muts$Tfix <- NA

d_muts[d_muts$fixGen != "N/A",]$Tfix <- as.integer(d_muts[d_muts$fixGen != "N/A",]$fixGen) - as.integer(d_muts[d_muts$fixGen != "N/A",]$originGen)  # Time to fixation

d_muts$sigma <- as.factor(d_muts$modelindex)

levels(d_muts$sigma) <- c("\u03c3 = 1", "\u03c3 = 10")

write.csv(d_muts, "d_muts.csv", row.names = F)

write.csv(d_means, "d_phenomeans.csv", row.names = F)


d_SFS <- read.csv("out_stabsel_sfs.csv", header = F)
names(d_SFS) <- c("gen", "seed", "modelindex", "mutID", "count") 
d_SFS$sigma <- as.factor(d_SFS$modelindex)

levels(d_SFS$sigma) <- c("\u03c3 = 1", "\u03c3 = 10")


write.csv(d_SFS, "d_SFS.csv", row.names = F)


# Plot means over time

plot_maker(d_means[d_means$gen > 84500,], type = "l", x = "gen", y = "phenomean", xlab = 
             "Time (generations)", ylab = "Mean phenotype", group = "seed", 
           facet = c("sigma", "h"), facet.lab = "Additive effect size distribution", leg.enabled = F,
           pal = "highlight")




plot_maker(d_SFS[d_SFS$gen > 84500 & d_SFS$gen < 87000,], "h", x = "count", xlab = "Allele count",
           facet = c("sigma", "h"), facet.lab = "Additive effect size distribution")

