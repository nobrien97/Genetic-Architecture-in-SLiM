# Script to plot heterozygosity over time


library(foreach)
library(doParallel)
library(future)

rsample <- as.character(runif(20, 1, (2^62 - 1)))
write.table(rbind(rsample), file = "seeds.csv", row.names = FALSE, col.names = FALSE, sep=",")

cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

#Run SLiM through the WSL, defining parameters in command line

setwd("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/SLiM/")

foreach(i=rsample) %dopar% {
    # Use string manipulation functions to configure the command line args, feeding from a data frame of latin square combinations
    # then run SLiM with system(),
    slim_out <- system(sprintf("wsl slim -s %s heterozygosity_burnin_test.slim", 
                               as.character(i), intern=T))
  }
stopCluster(cl)

system.time(system("wsl slim -s 123 heterozygosity_burnin_test.slim", 
                           ))

library(ggplot2)
library(gghighlight)

# Simple function to calculate standard error

std.error <-  function(x) {
  n <- length(x)
  sd <- sd(x)
  sd/sqrt(n)
}

dat_het <-read.csv("heterozygosity.csv", header = F)

names(dat_het) <- c("gens", "seed", "Ne", "rwide", "delmu", "mean", "var")
dat_het$delmu <- factor(dat_het$delmu, labels = c("No~deleterious~mutations", "Equal~chance~deleterious/neutral")) 
dat_het$rwide <- factor(dat_het$rwide, labels = c("r~'='~0", "r~'='~1.346~'x'~10^{-5}")) 
dat_het$seed <- as.factor(dat_het$seed)


# Calculate group means and standard errors to plot

library(dplyr)

dat_het_means <- dat_het[c(1, 3:6)] %>%
                   group_by(Ne, delmu, rwide, gens) %>%
                   summarise_all(list(groupmean = mean, se = std.error))
  

library(coda)

mcmc_het <- mcmc(dat_het$mean[dat_het$Ne == 10000 & dat_het$rwide > 0.0 & dat_het$delmu == 0.0])
autocorr.plot(mcmc_het, lag.max = 100)
autocorr.plot(dat_het$mean, lag.max = 50)


seedchoice <- sample(unique(dat_het$seed), 2)


 
 plot_het <- ggplot(dat_het_means,
                         aes(x = gens, y = groupmean, color = as.factor(Ne))) +
    geom_ribbon(aes(ymin = (groupmean - se), ymax = (groupmean + se)), alpha=0.3, show.legend = F, linetype=0) +
    geom_line() +
    facet_grid(delmu ~ ., scales = "free", labeller = label_parsed) +
    ggtitle("Mean heterozygosity across QTL loci over time") +
#    theme_classic() +
#    theme(legend.position = "right") +
    labs(x = "Generation", y = "Mean heterozygosity", color = "Population size")
  
