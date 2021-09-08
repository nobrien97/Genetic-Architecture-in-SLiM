df_means <- read.csv("out_stabsel_means.csv", header = F)
names(df_means) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", "dist", "mean_w", "deltaPheno")

df_means <- na.omit(df_means)

# Load combos
combos <- read.csv("../../../../src/PBS/GeneticConstraints/Sept2021_GadiRun/R/lhc_gaditest.csv")

# Add predictor variable columns
df_means[names(combos)[2:8]] <- NA

# Can this be vectorised? 
# Add the right LHC samples to the dataframe
for (index in unique(df_means$modelindex)) {
  df_means[df_means$modelindex == index, names(combos)[2:8]] <- combos[index, 2:8]
}

# Rename constraints
df_means[df_means$K == "\"c(1.0, 0.0, 0.0)\"", "K"] <- "Low"
df_means[df_means$K == "\"c(0.0, 1.0, 0.0)\"", "K"] <- "Medium"
df_means[df_means$K == "\"c(0.0, 0.0, 1.0)\"", "K"] <- "High"


# Burn-in data

df_burnin <- read.csv("out_stabsel_burnin.csv", header = F)
names(df_burnin) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean")


# Add predictor variable columns
df_burnin[names(combos)[2:8]] <- NA

# Can this be vectorised? If not RCPP
# Add the right LHC samples to the dataframe
for (index in unique(df_burnin$modelindex)) {
  df_burnin[df_burnin$modelindex == index, names(combos)[2:8]] <- combos[index, 2:8]
}

# Rename constraints
df_burnin[df_burnin$K == "\"c(1.0, 0.0, 0.0)\"", "K"] <- "Low"
df_burnin[df_burnin$K == "\"c(0.0, 1.0, 0.0)\"", "K"] <- "Medium"
df_burnin[df_burnin$K == "\"c(0.0, 0.0, 1.0)\"", "K"] <- "High"



# Import mutation output

df_muts <- read.csv("out_stabsel_muts.csv", header = F)
names(df_muts) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", "constraint", "originGen", "value", "chi", "Freq", "Count", "fixGen")



# Combine the two dataframes with dplyr
library(tidyverse)
df_combined <- inner_join(df_muts, df_means, by = c("gen", "modelindex", "seed"))


# For the muts df alone, we add the rest of the predictors
# Add predictor variable columns
df_muts[names(combos)[2:8]] <- NA

# Can this be vectorised? 
# Add the right LHC samples to the dataframe
for (index in unique(df_muts$modelindex)) {
  df_muts[df_muts$modelindex == index, names(combos)[2:8]] <- combos[index, 2:8]
}

# Rename constraints
df_muts[df_muts$K == "\"c(1.0, 0.0, 0.0)\"", "K"] <- "Low"
df_muts[df_muts$K == "\"c(0.0, 1.0, 0.0)\"", "K"] <- "Medium"
df_muts[df_muts$K == "\"c(0.0, 0.0, 1.0)\"", "K"] <- "High"






# Write files
write.csv(df_burnin, "df_burnin.csv", row.names = FALSE)
write.csv(df_means, "df_means.csv", row.names = FALSE)
write.csv(df_muts, "df_muts.csv", row.names = FALSE)
write.csv(df_combined, "df_combined.csv", row.names = FALSE)
