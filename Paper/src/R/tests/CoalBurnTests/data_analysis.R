# This data is from coalescent burn in tests - interested in heterozygosity over time, patterns of diversity and allele frequencies
# locisigma = (1, 10]
# genetic constraint = ((1, 0, 0), (0, 1, 0), (0, 0, 1)] - all low, all medium, all highly constrained
# rwide = 0.0, 0.2, 0.5
setwd("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/Paper/data/tests/CoalBurnTests")


d_burnmeans <- read.csv("out_stabsel_burnin.csv", header = F)

names(d_burnmeans) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean")

combos <- read.csv("../../../src/PBS/GeneticConstraints/CoalBurnTests/R/lhc_coalburntest.csv")


# Add predictor variable columns
d_burnmeans[names(combos)[2:4]] <- NA

# Can this be vectorised? If not RCPP
# Add the right LHC samples to the dataframe
for (index in unique(d_burnmeans$modelindex)) {
  d_burnmeans[d_burnmeans$modelindex == index, names(combos)[2:4]] <- combos[index, 2:4]
}


d_means <- read.csv("out_stabsel_means.csv", header = F)

names(d_means) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", "dist", "mean_w", "deltaPheno")

# Clean data (didn't remove any rows but fixed error with the below for loop ??????)
d_means <- d_means[!is.na(d_means$modelindex),]

# Add predictor variable columns
d_means[names(combos)[2:4]] <- NA

# Can this be vectorised? 
# Add the right LHC samples to the dataframe
for (index in unique(d_means$modelindex)) {
  d_means[d_means$modelindex == index, names(combos)[2:4]] <- combos[index, 2:4]
}



d_means$sigma <- 1

d_means$constraint <- "Low"

d_means$constraint <- as.factor(d_means$constraint)




d_muts <- read.csv("out_stabsel_muts.csv", header = F) # Ignore extra empty column (fixed in next version of script)

names(d_muts) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", "constraint", "originGen", "value", "chi", "Freq", "fixGen")


# Combine the two dataframes with dplyr
library(tidyverse)
df_combined <- inner_join(d_muts, d_means, by = c("gen", "modelindex", "seed"))

write.csv(df_combined, "df_combined.csv", row.names = FALSE)
write.csv(d_burnmeans, "df_burnin.csv", row.names = FALSE)


d_muts$constraint <- as.factor(d_muts$constraint)

levels(d_muts$constraint) <- c("Low", "Medium", "High")

d_muts$Tfix <- NA

d_muts[d_muts$fixGen != "N/A",]$Tfix <- as.integer(d_muts[d_muts$fixGen != "N/A",]$fixGen) - d_muts[d_muts$fixGen != "N/A",]$originGen  # Time to fixation

d_muts$sigma <- as.factor(d_muts$modelindex)

levels(d_muts$sigma) <- c("\u03c3 = 1", "\u03c3 = 10")


write.table(d_means, "d_phenomeans.csv", sep = ",", row.names = F)