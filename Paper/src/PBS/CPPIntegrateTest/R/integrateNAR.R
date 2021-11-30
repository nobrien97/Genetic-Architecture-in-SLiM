# Testing the numerical solutions to a NAR ODE from deSolve and SLiM
# to make sure everything works as expected.

setwd("/mnt/c/GitHub/Genetic-Architecture-in-SLiM/Paper/src/PBS/CPPIntegrateTest/")
# Load in combos 
lhc_samples <- read.csv("R/lhc_cppAUC.csv")
for (midx in 1:nrow(lhc_samples)) {
  
  combos <- lhc_samples[midx,]
  combos$Xstart <- round(combos$Xstart, 1) # Round to nearest .1
  combos$Xstop <- round(combos$Xstop, 1)
  
  if (combos$cpp) {
    system(sprintf("slim -d modelindex=%i -d Xstart=%f -d Xstop=%f -d Aalpha=%f -d Abeta=%f -d Balpha=%f -d Bbeta=%f -d Hilln=%f -d Bthreshold=%f ./slim/integrateCPP.slim",
                   midx, combos$Xstart, 
                   combos$Xstop, combos$Aalpha, 
                   combos$Abeta, combos$Balpha, combos$Bbeta,
                   combos$Hilln, combos$Bthreshold,
                   intern=F))
    
  } else {
    # Mostly taken from Stella Knief's honours project: 
    # https://github.com/sknief/honours/blob/master/3_HPC/OutbackRuns/ODE/ODE_SCRIPT.r
    library(tidyverse)
    library(deSolve)
    library(DescTools)
    Freya <-function(t, state, parameters) {
      with(as.list(c(state,parameters)), {
        dA <- Abeta * (t > Xstart && t <= Xstop) * 1/(1 + A^Hilln) - Aalpha*A
        dB <- Bbeta * A^Hilln/(Bthreshold^Hilln + A^Hilln) - Balpha*B
        list(c(dA, dB))
      })
    }
    
    #Values for the ODE
    state <- c(A = 0, B = 0)
    times <- seq(0, 10, by = 0.1)
    
    params <- c(Xstart = combos$Xstart, Xstop = combos$Xstop, Aalpha = combos$Aalpha, 
                Abeta = combos$Abeta, Balpha = combos$Balpha, Bbeta = combos$Bbeta,
                Hilln = combos$Hilln, Bthreshold = combos$Bthreshold)
    
    solution <- ode(state, times, Freya, params) %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(X = ifelse(time > params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
      select(time, X, A, B)
    
    out <- cbind(midx, AUC(solution$time, solution$A, absolutearea = T),
                 AUC(solution$time, solution$B, absolutearea = T))
    colnames(out) = c("modelindex", "A", "B")
    
    write.table(out, file = "out_AUC_R.tsv", append = T, sep = "\t", row.names = F,
                col.names = F)
    
  }
}

# Combine data frames
out_slim <- read.table("out_AUC_slim.tsv", sep = "\t", header = F)
out_R <- read.table("out_AUC_R.tsv", sep = "\t", header = F)

names(out_slim) <- c("modelindex", "A", "B")
names(out_R) <- c("modelindex", "A", "B")

out_slim[names(lhc_samples)[2:4]] <- NA
for (index in unique(out_slim$modelindex)) {
  out_slim[out_slim$modelindex == index, names(lhc_samples)] <- lhc_samples[index, ]
}

out_R[names(lhc_samples)[2:4]] <- NA
for (index in unique(out_R$modelindex)) {
  out_R[out_R$modelindex == index, names(lhc_samples)] <- lhc_samples[index, ]
}

df_combined <- rbind(out_slim, out_R)
df_combined$modelindex[df_combined$modelindex > 512] <- 1:512


# Scale data

df_combined$AScaled <- scale(df_combined$A)
df_combined$BScaled <- scale(df_combined$B)

# Plot distances?
df_combined$distA <- NA
df_combined$distB <- NA

for (idx in unique(df_combined$modelindex)) {
  df_combined[idx,]$distA <- c(dist(df_combined$A[df_combined$modelindex == idx]))
  df_combined[idx,]$distB <- c(dist(df_combined$B[df_combined$modelindex == idx]))
  df_combined[idx,]$scaleddistA <- c(dist(df_combined$AScaled[df_combined$modelindex == idx]))
  df_combined[idx,]$scaleddistB <- c(dist(df_combined$BScaled[df_combined$modelindex == idx]))
}

# Calculate maximum distances

max(df_combined$distA[!is.na(df_combined$distA)])
max(df_combined$distB[!is.na(df_combined$distB)])

write.table(df_combined, file = "out_AUC.tsv", sep = "\t", row.names = F)
  