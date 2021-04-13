# This file looks at data from locisigma = 10.0 and locisigma 1.0, with sampling increased to once every 
# 20 and 100 generations for population means and mutations, respectively 
# modelindex: 1 = locisigma = 1.0
#             2 = locisigma = 10.0

setwd("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/Paper/data/tests/GeneticConstraint/NewShapeIncreasedSamples//")

d_burnmeans <- read.csv("out_stabsel_burnin.csv", header = F)

names(d_burnmeans) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean")

d_means <- read.csv("out_stabsel_means.csv", header = F)

names(d_means) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", "dist", "mean_w")

d_muts <- read.csv("out_stabsel_muts.csv", header = F) # Ignore extra empty column (fixed in next version of script)

names(d_muts) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", "constraint", "originGen", "value", "chi", "Freq", "fixGen")

d_muts$constraint <- as.factor(d_muts$constraint)

levels(d_muts$constraint) <- c("Low", "Medium", "High")

d_muts$Tfix <- NA

d_muts[d_muts$fixGen != "N/A",]$Tfix <- as.integer(d_muts[d_muts$fixGen != "N/A",]$fixGen) - d_muts[d_muts$fixGen != "N/A",]$originGen  # Time to fixation


source("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/Paper/src/R/includes/plot_function.R")


# End of the run: after the adaptive walk, and some maintenance of variation
plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 85000,], type = "d", x="value", xlab = 
             "Additive effect size", colour.lab = "Genetic Constraint", colour="constraint", pal = "Magenta")


# Took about 1000 generations for all replicates to adapt

test_plot <- plot_maker(d_means, type = "l", x = "gen", y = "phenomean", xlab = 
             "Time (generations)", ylab = "Mean phenotype", group = "seed", colour = "as.factor(modelindex)", leg.enabled = F)



# Plot effect size at generation 75000 (pre-selection)

plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 75000,], type = "d", x="value", xlab = 
             "Additive effect size", facet = c("modelindex", "h"), colour.lab = "Genetic Constraint", colour="constraint", pal = "Magenta")

# Plot effect size following adaptation (generation 80000)

d_muts$sigma <- as.factor(d_muts$modelindex)
levels(d_muts$sigma) <- c("\u03B1 = 1", "\u03B1 = 10")

plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 80000,], type = "d", x="value", 
           facet = c("sigma", "v"), facet.lab = "Additive effect size distribution", 
           xlab = "Additive effect size", colour.lab = "Genetic Constraint",
           leg.pos = "bottom",
           colour="constraint", pal = "Magenta",
           savename = "gencon_fixations_dens_sigma.png")

plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 80000,], type = "d", x="chi", 
           facet = c("sigma", "v"), facet.lab = "Additive effect size distribution", 
           xlab = "Standardised effect size", colour.lab = "Genetic Constraint",
           leg.pos = "bottom",
           colour="constraint", pal = "Magenta",
           savename = "gencon_fixations_dens_stdsigma.png")

plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 80000,], type = "d", x="chi", 
           facet = c("sigma", "v"), facet.lab = "Additive effect size distribution", 
           xlab = "Standardised effect size", colour.lab = "Genetic Constraint",
           leg.pos = "bottom", dens.count = T,
           colour="constraint", pal = "Magenta",
           savename = "gencon_fixations_denscount_stdsigma.png")

plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 80000,], type = "d", x="abs(chi)", 
           facet = c("sigma", "v"), facet.lab = "Additive effect size distribution", 
           xlab = "Absolute standardised effect size", colour.lab = "Genetic Constraint",
           leg.pos = "bottom",
           colour="constraint", pal = "Magenta",
           savename = "gencon_fixations_dens_absstdsigma.png")



plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 80000,], type = "h", x="chi", xlab = 
             "Standardised effect size", colour.lab = "Genetic Constraint", colour="constraint", pal = "Magenta", 
           savename = "gencon_fixations_hist.png")


