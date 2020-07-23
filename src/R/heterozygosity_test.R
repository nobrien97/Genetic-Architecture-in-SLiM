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

dat_het <-read.csv("heterozygosity.csv", header = F)

names(dat_het) <- c("gens", "seed", "mean", "var")
#dat_het$gens = seq(500, (length(dat_het$mean) * 500), 500)


library(coda)

mcmc_het <- mcmc(dat_het$mean)
autocorr.plot(mcmc_het, lag.max = 50)
autocorr.plot(dat_het$mean, lag.max = 50)

EPibounds <- as.data.frame(c(0.0252+(0.0252*0.05), 0.0252-(0.0252*0.05)))
names(EPibounds) <- "pibounds"

plot_het_mean <- ggplot(dat_het,
                        aes(x = gens, y = mean, color = as.factor(seed))) +
  geom_line() +
  geom_hline(
    data = EPibounds,
    aes(yintercept = pibounds), linetype="dashed"
  ) +
  ggtitle("Mean heterozygosity across QTL loci over time") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Generation", y = "Mean heterozygosity")
  
plot_het_var <- ggplot(dat_het,
                        aes(x = gens, y = var, color = as.factor(seed))) +
  geom_line() +
  ggtitle("Variance in heterozygosity across QTL loci among individuals over time") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Generation", y = "Heterozygosity variance")
