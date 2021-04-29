# Generate combinations of variables and save as a csv

#parameters
p <- list()
p$locisigma <- c(1.0, 10.0)
p$gencon <- c('"c(1.0, 0.0, 0.0)"', '"c(0.0, 1.0, 0.0)"', '"c(0.0, 0.0, 1.0)"')
p$rwide <- c(0.0, 0.2, 0.5)

write.csv(expand.grid(p), "combos.csv")
