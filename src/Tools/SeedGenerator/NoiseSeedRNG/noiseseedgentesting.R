# Testing noise-based seed generator

df <- read.csv("seeds.csv")

# Test uniformity
hist(df$Seed)
range(df$Seed)

system.time(system("/mnt/c/GitHub/Genetic-Architecture-in-SLiM/src/Tools/SeedGenerator/NoiseSeedRNG/noiseseedgenerator -l -n 10000000"))
system.time(system("/mnt/c/GitHub/Genetic-Architecture-in-SLiM/src/Tools/SeedGenerator/Linux/seedgenerator -l -n 10000000"))
