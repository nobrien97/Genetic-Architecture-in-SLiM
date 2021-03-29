
# Import data

d_burnmeans <- read.csv("out_stabsel_burnin.csv", header = F)

names(d_burnmeans) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean")

d_means <- read.csv("out_stabsel_means.csv", header = F)

names(d_means) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", "dist", "mean_w")

d_muts <- read.csv("out_stabsel_muts.csv", header = F)[,1:10] # Ignore extra empty column (fixed in next version of script)

names(d_muts) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", "originGen", "value", "Freq", "fixGen")

d_pos <- read.csv("out_stabsel_pos.csv", header = F)

names(d_pos) <- c("modelindex", "seed", "nloci", "QTL1", "QTL2", "QTL3")

d_dict <- read.csv("out_stabsel_dict.csv", header = F)

names(d_dict) <- c("modelindex", "seed", "pos1", "pos2", "pos3", "pos4", "sep", "Type1", "Type2", "Type3", "Type4")


# need to add new column to muts: constraint - do so by pivot_longer d_dict, sort by pos, then use QTL position 
# of mut to get the right mut

d_dict_long <- d_dict %>% pivot_longer(cols = c("pos1", "pos2", "pos3", "pos4"), 
                            values_to = "pos")

d_dict_long$name <- rep(c(0, 1, 2, 2), 24)

d_dict_long <- d_dict_long[,-c(3:7)]


d_dict_long <- d_dict_long[order(d_dict_long$seed, d_dict_long$pos),]


d_muts$constraints <- -1


for (seed in unique(d_dict_long$seed)) {
    mini_frame <- d_dict_long[d_dict_long$seed == seed,]
    
      for (pos in unique(d_muts[d_muts$seed == seed,]$position)) {
        d_muts[d_muts$seed == seed & d_muts$position == pos,]$constraints <- mini_frame[mini_frame$pos == pos,]$name
      }
    
}

library(tidyverse)


# In this example, d_dict and d_pos show that this seed has QTLs at pos 0, 2, 3; and those are constrained according to 
# medium, low, high constraint respectively

ggplot(d_muts[d_muts$fixGen != "N/A",], 
       aes(value)) +
  facet_grid(~constraints) +
  geom_freqpoly() +
  theme_classic() +
  labs(x = "Trait Value", y = "Density")
