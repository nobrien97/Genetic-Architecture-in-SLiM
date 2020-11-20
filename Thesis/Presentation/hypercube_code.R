library(Cairo)
library(plot3D)
setwd("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM1")

lscombos_null <- read.csv("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM1/lscombos_null.csv")
lscombos_sel <- read.csv("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM3/lscombos_sel.csv")

lscombos_sel$X <- lscombos_sel$X + 1024
lscombos_null$tau <- 0.0

lscombos <- data.table::rbindlist(list(lscombos_null, lscombos_sel))

CairoWin()

size <- par(cex = 1)
LHCscatter <- scatter3D(lscombos$rwide, lscombos$locisigma, lscombos$delmu, colvar = NULL, xlab = "", 
                        ylab = "", zlab = "", phi = 2
                        )

CairoPNG("LHCscatter.png", LHCscatter)

