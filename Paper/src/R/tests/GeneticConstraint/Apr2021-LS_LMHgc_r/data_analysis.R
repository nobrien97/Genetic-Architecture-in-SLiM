# This data is from April 2021's test run of adaptation, with 18 models and 48 seeds
# locisigma = (1, 10]
# genetic constraint = ((1, 0, 0), (0, 1, 0), (0, 0, 1)] - all low, all medium, all highly constrained
# rwide = 0.0, 0.2, 0.5
setwd("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/Paper/data/tests/GeneticConstraint/Apr2021-LS_LMHgc_r")


d_burnmeans <- read.csv("out_stabsel_burnin.csv", header = F)

names(d_burnmeans) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean")

d_means <- read.csv("out_stabsel_means.csv", header = F)

names(d_means) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", "dist", "mean_w", "deltaPheno")

d_means$sigma <- 0.0

d_means[d_means$modelindex %% 2 != 0,]$sigma <- 1
d_means[d_means$modelindex %% 2 == 0,]$sigma <- 10

d_means$sigma <- as.factor(d_means$sigma)

levels(d_means$sigma) <- c("\u03c3 = 1", "\u03c3 = 10")

d_means$gencon <- 0

d_means[d_means$modelindex %in% c(1, 2, 7, 8, 13, 14),]$gencon <- "Low"
d_means[d_means$modelindex %in% c(3, 4, 9, 10, 15, 16),]$gencon <- "Medium"
d_means[d_means$modelindex %in% c(5, 6, 11, 12, 17, 18),]$gencon <- "High"


d_means$gencon <- as.factor(d_means$gencon)


d_means$recomRate <- 0.0

d_means[d_means$modelindex %in% c(7:12),]$recomRate <- 0.2
d_means[d_means$modelindex %in% c(13:18),]$recomRate <- 0.5


# Plot the delta fitness over time

loopseq <- seq(1, length(unique(d_means$gen)))

d_means$delta_w <- 0
d_means <- d_means %>% arrange(seed, modelindex, gen)

for (seed in unique(d_means$seed)) {
  for (index in unique(d_means$modelindex)) {
    w_vec <- c(0.0, d_means[d_means$modelindex == index & d_means$seed == seed,]$mean_w[-length(d_means[d_means$modelindex == index & d_means$seed == seed,]$mean_w)])
        d_means[d_means$modelindex == index & d_means$seed == seed,]$delta_w <-
          d_means[d_means$modelindex == index & d_means$seed == seed,]$mean_w - 
          w_vec
  }
}

# Always a good idea to == floating points, what are you talking about?
d_means[d_means$delta_w == 1.000,]$delta_w <- NA

d_means <- d_means[,c(1:3, 10:12, 4:9, 13)]

write.table(d_means, "d_phenomeans.csv", sep = ",", row.names = F)



d_muts <- read.csv("out_stabsel_muts.csv", header = F) # Ignore extra empty column (fixed in next version of script)

names(d_muts) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", "constraint", "originGen", "value", "chi", "Freq", "fixGen")

d_muts$constraint <- as.factor(d_muts$constraint)

levels(d_muts$constraint) <- c("Low", "Medium", "High")

d_muts$Tfix <- NA

d_muts[d_muts$fixGen != "N/A",]$Tfix <- as.integer(d_muts[d_muts$fixGen != "N/A",]$fixGen) - d_muts[d_muts$fixGen != "N/A",]$originGen  # Time to fixation

d_muts$sigma <- as.factor(d_muts$modelindex)

levels(d_muts$sigma) <- c("\u03c3 = 1", "\u03c3 = 10")


source("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/Paper/src/R/includes/plot_function.R")

