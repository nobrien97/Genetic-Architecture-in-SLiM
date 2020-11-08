# Script to plot heterozygosity over time


library(foreach)
library(doParallel)
library(future)

rsample <- as.character(runif(20, 1, (2^62 - 1)))
write.table(cbind(rsample), file = "seeds.csv", row.names = FALSE, col.names = FALSE, sep=",")

rsample <- read.csv("seeds.csv", header = F)
rsample <- as.character(as.vector(t(rsample)))


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

dat_het <-read.csv("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/src/SLiM/Heterozygosity/heterozygosity.csv", header = F)


names(dat_het) <- c("gens", "seed", "Ne", "rwide", "delmu", "mean", "var")
dat_het$delmu <- as.factor(dat_het$delmu)
dat_het$rwide <- factor(dat_het$rwide, levels = c("r~'='~0", "r~'='~1.346~'x'~10^{-5}")) 
dat_het$seed <- as.factor(dat_het$seed)


# Calculate group means and standard errors to plot

library(tidyverse)

dat_het_means <- dat_het[c(1, 3:6)] %>%
                   group_by(Ne, delmu, rwide, gens) %>%
                   summarise_all(list(groupmean = mean, se = std.error))

# Labels for expected equilibrium points for each population size (4Nu)

EPibounds <- as.data.frame(c(4*500*6.3e-6+((4*500*6.3e-6)*0.05), 4*500*6.3e-6-((4*500*6.3e-6)*0.05),  
                             4*2500*6.3e-6+((4*2500*6.3e-6)*0.05), 4*2500*6.3e-6-((4*2500*6.3e-6)*0.05),
                             4*10000*6.3e-6+((4*10000*6.3e-6)*0.05), 4*10000*6.3e-6-((4*10000*6.3e-6)*0.05)) )
EPibounds$delmu <- as.factor("No~deleterious~mutations")
names(EPibounds) <- "pibounds"
 
plot_het <- ggplot(dat_het_means[dat_het_means$delmu == 0,],
                         aes(x = gens, y = groupmean, color = as.factor(Ne))) +
    geom_ribbon(aes(ymin = (groupmean - se), ymax = (groupmean + se)), alpha=0.3, show.legend = F, linetype=0) +
    geom_line() +
    geom_hline(
    data = EPibounds,
    aes(yintercept = pibounds), linetype="dashed"
  ) +
    theme_classic() +
    scale_color_manual(values = c("red", "purple", "blue")) +
#    theme_classic() +
#    theme(legend.position = "right") +
    labs(x = "Generation", y = "Mean heterozygosity", color = "Population size") +
    theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

ggsave("figs1_het.png", plot_het, height = 8, width = 8, dpi = 800)



