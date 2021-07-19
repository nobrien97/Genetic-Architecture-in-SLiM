# Generate combinations of variables and save as a csv

library(DoE.wrapper)

lhc_seed = sample(1:2147483647, 1)

lhc_samples <- lhs.design(
    nruns = 1024,
    nfactors = 6,
    seed = lhc_seed,
    type = "maximin",
    factor.names = list(
        rwide = c(0.0, 0.5),
        locisigma = c(0.1, 10),
        Ne = c(100, 20000),
        nloci = c(1, 100),
        width = c(0.0002, 0.5),
        opt = c(1, 100))
)

# Check correlations and space between samples

plot(lhc_samples)
cor(lhc_samples)
max(abs(cor(lhc_samples)[cor(lhc_samples) != 1]))
hist(lhc_samples$locisigma)

# Add on column for constraint - will need to repeat the sampling 3 times for K = Low, Medium, High

K_levels <- c('"c(1.0, 0.0, 0.0)"', '"c(0.0, 1.0, 0.0)"', '"c(0.0, 0.0, 1.0)"')
lhc_samples <-  lhc_samples[rep(seq_len(nrow(lhc_samples)), 3), ]
lhc_samples$K <- rep(K_levels, each = nrow(lhc_samples)/3)

write.csv(lhc_samples, "./lhc_gaditest.csv", header = T)



# TODO: Have more of a look into this - I don't think it's necessary for purely reading files
# Might be worth combining model outputs at the end though instead of writing to one file
# Split data frame into separate files for each row - 
# to avoid filesystem problems from reading the same file on multiple threads
if (!dir.exists("./combos"))
    dir.create("./combos")

for (row in seq.int(1, nrow(lhc_samples))) {
    write.csv(lhc_samples[row,], paste0("./combos/lhc_sample", row, ".csv"), row.names = F)
}

# Do the same for seeds:

seeds <- read.csv("./May2021-GadiTest/R/seeds_gaditest.csv")

if (!dir.exists("./seeds"))
    dir.create("./seeds")

for (row in seq.int(1, nrow(seeds))) {
    write.table(seeds[row,], paste0("./seeds/seed", row, ".csv"), row.names = F, col.names = F, sep = ",")
}




