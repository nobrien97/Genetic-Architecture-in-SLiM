
# Import data

setwd("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/Paper/data/tests/GeneticConstraint/noLD/")

d_burnmeans <- read.csv("out_stabsel_burnin.csv", header = F)

names(d_burnmeans) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean")

d_means <- read.csv("out_stabsel_means.csv", header = F)

names(d_means) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", "dist", "mean_w")

d_muts <- read.csv("out_stabsel_muts.csv", header = F) # Ignore extra empty column (fixed in next version of script)

names(d_muts) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", "constraint", "originGen", "value", "Freq", "fixGen")

d_muts$constraint <- as.factor(d_muts$constraint)

levels(d_muts$constraint) <- c("Low", "Medium", "High")


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
library(colorspace)

# In this example, d_dict and d_pos show that this seed has QTLs at pos 0, 2, 3; and those are constrained according to 
# medium, low, high constraint respectively

plot_fixations <- ggplot(d_muts[d_muts$fixGen != "N/A" & d_muts$mutType == 3 & d_muts$gen == 150000,], 
                          aes(x = value, colour = constraint)) +
                          scale_color_discrete_sequential() +
                          geom_density() +
                          theme_classic() +
                          ggtitle("Distribution of effect sizes among fixations") +
                          labs(x = "Trait Value", y = "Density", colour = "Genetic Constraint") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        legend.key.width = unit(1.6, "cm"),
        legend.title.align = 0.5,
        panel.spacing.y = unit(1, "lines"))

plot_fixations

ggsave("plot_fixations.png", plot = plot_fixations, height = 15, width = 15)


plot_seg <- ggplot(d_muts[d_muts$fixGen == "N/A" & d_muts$mutType == 3 & d_muts$gen == 150000,], 
              aes(x = value, colour = constraint)) + 
              scale_color_discrete_sequential() +
              geom_density() +
              theme_classic() +
              ggtitle("Distribution of effect sizes among segregating mutations") +
              labs(x = "Trait Value", y = "Density", colour = "Genetic Constraint") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        legend.key.width = unit(1.6, "cm"),
        legend.title.align = 0.5,
        panel.spacing.y = unit(1, "lines"))

plot_seg
ggsave("plot_seg.png", plot = plot_seg, height = 15, width = 15)




plot_del_fixations <- ggplot(d_muts[d_muts$fixGen != "N/A" & d_muts$mutType == 2,], 
                            aes(x = value, colour = constraint)) +
                            scale_color_discrete_sequential() +
                            geom_density() +
                            theme_classic() +
                            ggtitle("Distribution of deleterious fitness effects among fixations") +
                            labs(x = "Selection coefficient", y = "Density", colour = "Genetic Constraint") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        legend.key.width = unit(1.6, "cm"),
        legend.title.align = 0.5,
        panel.spacing.y = unit(1, "lines"))


plot_del_fixations
ggsave("plot_del_fixations.png", plot = plot_del_fixations, height = 15, width = 15)



plot_del_seg <- ggplot(d_muts[d_muts$fixGen == "N/A" & d_muts$mutType == 2,], 
                   aes(x = value, colour = constraint)) + 
  scale_color_discrete_sequential() +
  geom_density() +
  theme_classic() +
  ggtitle("Distribution of deleterious fitness effects among segregating mutations") +
  labs(x = "Selection coefficient", y = "Density", colour = "Genetic Constraint") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        legend.key.width = unit(1.6, "cm"),
        legend.title.align = 0.5,
        panel.spacing.y = unit(1, "lines"))

plot_del_seg
ggsave("plot_del_seg.png", plot = plot_del_seg, height = 15, width = 15)



plot_phenoburn <- ggplot(d_burnmeans, 
                         aes(x = gen, y = phenomean, group = seed)) +
  scale_color_discrete_sequential() +
  geom_line() +
  theme_classic() +
  ggtitle("Drift of phenotypic means over time") +
  labs(x = "Time (generations)", y = "Phenotypic mean") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        legend.key.width = unit(1.6, "cm"),
        legend.title.align = 0.5,
        panel.spacing.y = unit(1, "lines"))

plot_phenoburn


plot_pheno <- ggplot(d_means, 
                         aes(x = gen, y = phenomean, group = seed)) +
  scale_color_discrete_sequential() +
  geom_line() +
  theme_classic() +
  ggtitle("Movement of phenotypic means over time with stabilising selection") +
  labs(x = "Time (generations)", y = "Phenotypic mean") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        legend.key.width = unit(1.6, "cm"),
        legend.title.align = 0.5,
        panel.spacing.y = unit(1, "lines"))

plot_pheno



# Extreme values: Here, high constraint means neutral mutations are 50x and deleterious 100x more likely than QTL mutations
# Medium constraint: neutral are 5x more likely and deleterious 10x more likely than QTL
# low constraint: each type are equally likely

setwd("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/Paper/data/tests/GeneticConstraint/extremeconstraints/")

d_burnmeans_ex <- read.csv("out_stabsel_burnin.csv", header = F)

names(d_burnmeans_ex) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean")

d_means_ex <- read.csv("out_stabsel_means.csv", header = F)

names(d_means_ex) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", "dist", "mean_w")

d_muts_ex <- read.csv("out_stabsel_muts.csv", header = F) # Ignore extra empty column (fixed in next version of script)

names(d_muts_ex) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", "constraint", "originGen", "value", "Freq", "fixGen")

d_muts_ex$constraint <- as.factor(d_muts_ex$constraint)

levels(d_muts_ex$constraint) <- c("Low", "Medium", "High")

d_muts_ex$Tfix <- NA

d_muts_ex[d_muts_ex$fixGen != "N/A",]$Tfix <- as.integer(d_muts_ex[d_muts_ex$fixGen != "N/A",]$fixGen) - d_muts_ex[d_muts_ex$fixGen != "N/A",]$originGen  # Time to fixation


plot_ex_fixations <- ggplot(d_muts_ex[d_muts_ex$fixGen != "N/A" & d_muts_ex$mutType == 3 & d_muts_ex$gen == 150000,], 
                         aes(x = value, colour = constraint)) +
  scale_color_discrete_sequential() +
  geom_density() +
  theme_classic() +
  ggtitle("Distribution of effect sizes among fixations") +
  labs(x = "Trait Value", y = "Density", colour = "Genetic Constraint") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        legend.key.width = unit(1.6, "cm"),
        legend.title.align = 0.5,
        panel.spacing.y = unit(1, "lines"))

plot_ex_fixations


plot_ex_seg <- ggplot(d_muts_ex[d_muts_ex$fixGen == "N/A" & d_muts_ex$mutType == 3,], 
                   aes(x = value, colour = constraint)) + 
  scale_color_discrete_sequential() +
  geom_density() +
  theme_classic() +
  ggtitle("Distribution of effect sizes among segregating mutations") +
  labs(x = "Trait Value", y = "Density", colour = "Genetic Constraint") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        legend.key.width = unit(1.6, "cm"),
        legend.title.align = 0.5,
        panel.spacing.y = unit(1, "lines"))

plot_ex_seg



plot_ex_Tfix <- ggplot(d_muts_ex[!is.na(d_muts_ex$Tfix) & d_muts_ex$mutType == 3,], 
                       aes(x = Tfix, colour = constraint)) + 
  scale_color_discrete_sequential() +
  geom_density() +
  theme_classic() +
  ggtitle("Time to fixation among QTL mutations") +
  labs(x = "Time (generations)", y = "Density", colour = "Genetic Constraint") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        legend.key.width = unit(1.6, "cm"),
        legend.title.align = 0.5,
        panel.spacing.y = unit(1, "lines"))

plot_ex_Tfix


source("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/Paper/src/R/includes/plot_function.R")

plot_maker(d_muts_ex[!is.na(d_muts_ex$Tfix) & d_muts_ex$mutType == 3,], type = "d", x="Tfix", xlab = 
             "Time to fixation (generations)", colour.lab = "Genetic Constraint", colour="constraint")


# Extreme values w/ large initial distance to the optimum 

setwd("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/Paper/data/tests/GeneticConstraint/faropt/")

d_burnmeans_far <- read.csv("out_stabsel_burnin.csv", header = F)

names(d_burnmeans_far) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean")

d_means_far <- read.csv("out_stabsel_means.csv", header = F)

names(d_means_far) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", "dist", "mean_w")

d_muts_far <- read.csv("out_stabsel_muts.csv", header = F) # Ignore extra empty column (fixed in next version of script)

names(d_muts_far) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", "constraint", "originGen", "value", "Freq", "fixGen")

d_muts_far$constraint <- as.factor(d_muts_far$constraint)

levels(d_muts_far$constraint) <- c("Low", "Medium", "High")

d_muts_far$Tfix <- NA

d_muts_far[d_muts_far$fixGen != "N/A",]$Tfix <- as.integer(d_muts_far[d_muts_far$fixGen != "N/A",]$fixGen) - d_muts_far[d_muts_far$fixGen != "N/A",]$originGen  # Time to fixation




