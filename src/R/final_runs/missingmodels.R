################################################################################
#    For some reason, Nimrod didn't do 25 of the runs on AIM1_final_512.nsh    #
#      Now I have to find what the runs were so I can rerun them manually      #
################################################################################

# Import data
null_chr512 <- read.csv("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM1/Output/out_8T_null_chr_512.csv", header = F)

names(null_chr512)[1:2] <- c("modelindex", "seed")

length(unique(null_chr512$modelindex)) # Yep, we're missing stuff: 25 stuffs to be precise
length(unique(null_chr512$seed))

# Chop the dataframe to something more manageable, we only need seeds and models
# Then identify which models and seeds have stuff missing
null_chr512 <- null_chr512[c(1:2)]
library(plyr)
models <- count(null_chr512$modelindex)
models <- models$x[models$freq < 100] # Which models have been repeated less than 100 times (there are 100 seeds)?
seeds <- count(null_chr512$seed)
seeds <- seeds[seeds$freq < 256,]$x # which seeds have been repeated less than 256 times (there are 256 models)?


# Set up new frame to store our 25 rogue combinations in
id <- data.frame (modelindex = rep(0, 25), seed = rep(0, 25))

# For each combination of model and seed, figure out which combinations don't have an entry in null_chr512, then
# put them in id
for (i in seq_along(models)) {
  for (j in seq_along(seeds)) {
    if (nrow(subset(null_chr512, modelindex == models[i] & seed == seeds[j])) != 1) {
      id[i,]$modelindex <- models[i]
      id[i,]$seed <- seeds[j]
    }
  }
}

#> id
#   modelindex       seed
#1         263  283359780
#2         275  283359780
#3         297  283359780
#4         319  283359780
#5         330  283359780
#6         367  283359780
#7         380  283359780
#8         383  283359780
#9         400  283359780
#10        407  283359780
#11        409  283359780
#12        410 2414432728
#13        420 2414432728
#14        431 2414432728
#15        436  283359780
#16        440 2414432728
#17        448  283359780
#18        449  283359780
#19        459 2414432728
#20        466  283359780
#21        470 2414432728
#22        477  283359780
#23        489  283359780
#24        496  283359780
#25        498  283359780

# Will have to run these model indexes and seeds again because why not, sounds like good fun

write.csv(id, "AIM1_supp_combos.csv")


################################################################################
#          It also didn't do 24 of the runs on AIM3_final_112.nsh              #            
#      Now I have to find what the runs were so I can rerun them manually      #
################################################################################

# Import data
sel_chr112 <- read.csv("F:/Uni/AIM3/OUTPUT/out_8T_stabsel_chr_112.csv", header = F)

names(sel_chr112)[1:2] <- c("modelindex", "seed")

length(unique(sel_chr112$modelindex)) # Yep, we're missing stuff: 25 stuffs to be precise
length(unique(sel_chr112$seed))

# Chop the dataframe to something more manageable, we only need seeds and models
# Then identify which models and seeds have stuff missing
sel_chr112 <- sel_chr112[c(1:2)]
library(tidyverse)
models <- count(sel_chr112$modelindex)
models <- models$x[models$freq < 100] # Which models have been repeated less than 100 times (there are 100 seeds)?
seeds <- count(sel_chr112$seed)
seeds <- seeds[seeds$freq < 16,]$x # which seeds have been repeated less than 16 times (there are 16 models)?
modseeds <- crossing(models, seeds)

# Set up new frame to store our 24 rogue combinations in
id <- data.frame(model = 1, seed = 1)

# For each combination of model and seed, figure out which combinations don't have an entry in sel_chr112, then
# put them in id
for (i in seq_along(modseeds$models)) {
    if (nrow(subset(sel_chr112, modelindex == modseeds$models[i] & seed == modseeds$seeds[i])) != 1) {
      missing_new <- data.frame(model = modseeds$models[i], seed = modseeds$seeds[i])
      id <- rbind(id, missing_new)
    }
}

#> id
#model       seed
#1     97 1931309650
#2     99  749138265
#3     99 1111628040
#4     99 3522021942
#5     99 3925196109
#6    100  749138265
#7    100 1763963450
#8    101 2987488969
#9    102  749138265
#10   103  749138265
#11   103  884298832
#12   103 2160468435
#13   103 3056777668
#14   103 3735160950
#15   105 1535973786
#16   105 2787103171
#17   106 1931309650
#18   107 1931309650
#19   108  553958512
#20   108  884298832
#21   108 1111628040
#22   108 3165192143
#23   109 2987488969
#24   111 3165192143

# Will have to run these model indexes and seeds again because why not, sounds like good fun

write.csv(id, "AIM3_supp_combos.csv")
