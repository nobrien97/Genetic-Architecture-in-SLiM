# Graph LD - triangular matrix thing
library(tidyverse)

plotLD <- function(mat){
    LDTbl <- as.data.frame(mat) %>%
        mutate(Locus1=row.names(.)) %>%
        pivot_longer(-Locus1,names_to="Locus2",values_to="r2")

    plt <- ggplot(LDTbl, mapping = aes(x = Locus1, y = Locus2, fill = r2)) +
    geom_bin_2d() +
    labs(x = "Locus X", y = "Locus Y", fill = expression(r^2))
    return(plt)
}

setwd("/mnt/c/GitHub/Genetic-Architecture-in-SLiM/Paper/data/tests/CoalBurnTest2")
mateg <- read.table("out_stabsel_ld_burnin4002686600_55.tsv", header = F)
mateg <- as.matrix(mateg)
colnames(mateg) <- 1:100
rownames(mateg) <- 1:100

plotLD(mateg)
