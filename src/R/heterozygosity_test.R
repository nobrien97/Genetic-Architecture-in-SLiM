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
library(forcats)

# Simple function to calculate standard error

std.error <-  function(x) {
  n <- length(x)
  sd <- sd(x)
  sd/sqrt(n)
}

dat_het <-read.csv("heterozygosity.csv", header = F)

dat_het1000 <- read.csv("heterozygosity1000.csv", header = F)

names(dat_het1000) <- c("gens", "seed", "mean", "var")
dat_het1000$seed <- as.factor(dat_het1000$seed)

names(dat_het) <- c("gens", "seed", "Ne", "rwide", "delmu", "mean", "var")
dat_het$delmu <- factor(dat_het$delmu, levels = c('No~deleterious~mutations', 'Equal~chance~deleterious/neutral')) 
dat_het$rwide <- factor(dat_het$rwide, levels = c("r~'='~0", "r~'='~1.346~'x'~10^{-5}")) 
dat_het$seed <- as.factor(dat_het$seed)


# Calculate group means and standard errors to plot

library(dplyr)

dat_het_means <- dat_het[c(1, 3:6)] %>%
                   group_by(Ne, delmu, rwide, gens) %>%
                   summarise_all(list(groupmean = mean, se = std.error))

dat_het1000_means <- dat_het1000[c(1, 3)] %>%
                      group_by(gens) %>%
                      summarise_all(list(groupmean = mean, se = std.error))  

library(coda)

mcmc_het <- mcmc(dat_het$mean[dat_het$Ne == 10000 & dat_het$rwide > 0.0 & dat_het$delmu == 0.0])
autocorr.plot(mcmc_het, lag.max = 100)
autocorr.plot(dat_het$mean, lag.max = 50)

EPi1000bounds <- as.data.frame(c(0.0252+(0.0252*0.05), 0.0252-(0.0252*0.05)))
names(EPibounds) <- "pibounds"

plot_het1000 <- ggplot(dat_het1000_means,
                   aes(x = gens, y = groupmean)) +
  geom_ribbon(aes(ymin = (groupmean - se), ymax = (groupmean + se)), alpha=0.3, show.legend = F, linetype=0) +
  geom_line() +
  geom_hline(
    data = EPi1000bounds,
    aes(yintercept = pibounds), linetype="dashed"
  ) +
  ggtitle("Mean heterozygosity across QTL loci over time (Null model, Ne = 1000, No deleterious mutations)") +
  #    theme_classic() +
  #    theme(legend.position = "right") +
  labs(x = "Generation", y = "Mean heterozygosity", color = "Population size")


# Labels for expected equilibrium points for each population size (4Nu)

EPibounds <- as.data.frame(c(4*500*6.3e-6+((4*500*6.3e-6)*0.05), 4*500*6.3e-6-((4*500*6.3e-6)*0.05),  
                             4*2500*6.3e-6+((4*2500*6.3e-6)*0.05), 4*2500*6.3e-6-((4*2500*6.3e-6)*0.05),
                             4*10000*6.3e-6+((4*10000*6.3e-6)*0.05), 4*10000*6.3e-6-((4*10000*6.3e-6)*0.05)) )
EPibounds$delmu <- as.factor("No~deleterious~mutations")
names(EPibounds) <- "pibounds"
 
plot_het <- ggplot(dat_het_means,
                         aes(x = gens, y = groupmean, color = as.factor(Ne))) +
    geom_ribbon(aes(ymin = (groupmean - se), ymax = (groupmean + se)), alpha=0.3, show.legend = F, linetype=0) +
    geom_line() +
    geom_hline(
    data = EPibounds,
    aes(yintercept = pibounds), linetype="dashed"
  ) +
    facet_grid(delmu ~ ., scales = "free", labeller = label_parsed) +
    ggtitle("Mean heterozygosity across QTL loci over time (AIM 1 model)") +
#    theme_classic() +
#    theme(legend.position = "right") +
    labs(x = "Generation", y = "Mean heterozygosity", color = "Population size")

plot_het_500 <- ggplot(subset(dat_het_means, Ne == 500),
                   aes(x = gens, y = groupmean)) +
  geom_ribbon(aes(ymin = (groupmean - se), ymax = (groupmean + se)), alpha=0.3, show.legend = F, linetype=0) +
  geom_line() +
 
  facet_grid(delmu ~ ., scales = "free", labeller = label_parsed) +
  ggtitle("Mean heterozygosity across QTL loci over time (AIM 1 model), Ne = 500") +
  #    theme_classic() +
  #    theme(legend.position = "right") +
  labs(x = "Generation", y = "Mean heterozygosity")




 
 
dat_hetrecom <-read.csv("heterozygosity_recom.csv", header = F)
 
names(dat_hetrecom) <- c("gens", "seed", "Ne", "rwide", "delmu", "mean", "var")
dat_hetrecom$delmu <- factor(dat_hetrecom$delmu, labels = c('No~deleterious~mutations', 'Equal~chance~deleterious/neutral'), ordered = T) 
dat_hetrecom$seed <- as.factor(dat_hetrecom$seed)  

dat_hetrecom_means <- dat_hetrecom[c(1, 3:6)] %>%
  group_by(Ne, delmu, rwide, gens) %>%
  summarise_all(list(groupmean = mean, se = std.error))


plot_hetrecom <- ggplot(dat_hetrecom_means,
                   aes(x = gens, y = groupmean, color = as.factor(Ne))) +
  geom_ribbon(aes(ymin = (groupmean - se), ymax = (groupmean + se)), alpha=0.3, show.legend = F, linetype=0) +
  geom_line() +
  facet_grid(delmu ~ ., scales = "free_y", labeller = label_parsed) +
  geom_hline(
    data = EPibounds,
    aes(yintercept = pibounds), linetype="dashed"
  ) +
  ggtitle("Mean heterozygosity across QTL loci over time (AIM 2 model)") +
  #    theme_classic() +
  #    theme(legend.position = "right") +
  labs(x = "Generation", y = "Mean heterozygosity", color = "Population size")

plot_hetrecom_500 <- ggplot(subset(dat_hetrecom_means, Ne == 500),
                       aes(x = gens, y = groupmean)) +
  geom_ribbon(aes(ymin = (groupmean - se), ymax = (groupmean + se)), alpha=0.3, show.legend = F, linetype=0) +
  geom_line() +
  
  facet_grid(delmu ~ ., scales = "free", labeller = label_parsed) +
  ggtitle("Mean heterozygosity across QTL loci over time (AIM 2 model), Ne = 500") +
  #    theme_classic() +
  #    theme(legend.position = "right") +
  labs(x = "Generation", y = "Mean heterozygosity")


# Heterozygosity test on new data set, where H is generated from all mutations, not just QTLs
# Ne = 1000

setwd("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/Matrix_test_July/Output")
dat_matex <- read.csv("out_8T100L_null_means.csv", header = F)

names(dat_matex)[1:6] <- c("gens", "seed", "modelindex", "rsd", "rwide", "delmu") 
names(dat_matex)[35] <- "Ne"
names(dat_matex)[108:117] <- c("meanH1", "meanH2", "meanH3", "meanH4", "meanH5", "meanH6", "meanH7", "meanH8", "meanH9", "meanH10")
dat_matex <- dat_matex[,c(1:6, 35, 108:117)]
dat_matex$seed <- as.factor(dat_matex$seed)
dat_matex$meanH <- (dat_matex$meanH1 + dat_matex$meanH2 + dat_matex$meanH3 + dat_matex$meanH4 + dat_matex$meanH5 + dat_matex$meanH6 + dat_matex$meanH7 + dat_matex$meanH8 + dat_matex$meanH9 + dat_matex$meanH10)/10

dat_matex_means <- dat_matex[c(1, 6, 8:18)] %>%
  group_by(delmu, gens) %>%
  summarise_all(list(groupmean = mean, se = std.error))

plot_matex_H <- ggplot(dat_matex_means,
                            aes(x = gens, y = meanH_groupmean)) +
  geom_ribbon(aes(ymin = (meanH_groupmean - meanH_se), ymax = (meanH_groupmean + meanH_se)), alpha=0.3, show.legend = F, linetype=0) +
  geom_line() +
  
  facet_grid(delmu ~ ., scales = "free", labeller = label_parsed) +
  ggtitle("Mean heterozygosity across the genome over time (AIM 1 model), Ne = 1000") +
  #    theme_classic() +
  #    theme(legend.position = "right") +
  labs(x = "Generation", y = "Mean heterozygosity")



# Null Time test 3

plot_nulltime3 <- ggplot(null_time3,
                         aes(x=gen, y=time)) +
  geom_line()


# Recom time test 3

plot_recomtime3 <- ggplot(recom_time3,
                         aes(x=gen, y=time)) +
  geom_line()


# Stab sel time test 2: old model format

plot_stabseltime2 <- ggplot(slim_stabsel8T_timetest,
                          aes(x=gen, y=time)) +
  geom_line()
