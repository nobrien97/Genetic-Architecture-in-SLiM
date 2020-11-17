library(Cairo)
CairoWin()
LHCscatter <- scatter3D(lscombos$rwide, lscombos$locisigma, lscombos$delmu, colvar = NULL, ticktype = "detailed", xlab = "Recombination rate ", ylab = "Additive effect size variance", zlab = "QTL mutation rate", phi = 2)
CairoPNG("LHCscatter.png", LHCscatter)