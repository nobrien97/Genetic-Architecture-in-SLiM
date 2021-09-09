# Testing noise-based seed generator

df <- read.csv("seeds.csv")

# Test uniformity
hist(df$Seed)
range(df$Seed)
