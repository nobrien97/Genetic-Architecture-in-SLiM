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

write.csv(lhc_samples, "./R/lhc_samples.csv", row.names = F)

