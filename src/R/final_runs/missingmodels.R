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
