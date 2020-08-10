## Results for all worst-case runs of each model.
## Due to a bug/git hiccup, stab sel only ran for 100,000 generations, but the figure is linear, so I am 
## extrapolating to the final time.


library(ggplot2)

# Except where otherwise mentioned, the input parameters are as follows:

# -s 123 Ne=10000 nloci=500 locisigma=10.0 pleiorate=1.0 delmu=1.0 rwide=1.346e-5 pleio_cov=0.5 delchr=10 
# and for selection models: wsd=10.0


setwd("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/Speedtests/time_test_diag3")

# Import data

null_time3 <- read.table("null_time3.txt", header=T)
null_time4 <- read.table("null_time3_lowerrecom.txt", header=T)
recom_time3 <- read.table("recom_time3.txt", header=T)
stabsel_time3 <- read.table("stabsel_time3.txt", header=T)




# Null Time test 3 - recombination rate 1.241e-4

plot_nulltime3 <- ggplot(null_time3,
                         aes(x=gen, y=time)) +
  geom_line()

# Null Time test 4 - recombination rate 1.346e-5

plot_nulltime4 <- ggplot(null_time3_lowerrecom,
                         aes(x=gen, y=time)) +
  geom_line()


# Recom time test 3

plot_recomtime3 <- ggplot(recom_time3,
                          aes(x=gen, y=time)) +
  geom_line()


# Stab sel time test 3 (Ne = 10,000)

plot_stabseltime3 <- ggplot(stabsel_time3,
                            aes(x=gen, y=time)) +
       geom_line()


# Only ran for 100,000 generations, distinct difference between burn-in time per gen and the selection model,
# so we make a linear model based off only the selection model run data (without the burn in)

sel_sbst <- subset(stabsel_time3, gen > 50000)

lm(time ~ gen, data = sel_sbst)

# 3.305x + -123000
# for 150,000 gens = 103.542 hours - assumes first 50k gens were at the same rate, they weren't

# So lm(selsbst(150,000)) - lm(selsbst(50,000)) + burn-in time for more accurate estimate

# (372750 - 42250) + 42562.7 = 373062.7s = 103.629 hrs = 4.32 days


# Stab sel time test Ne = 8000

plot_stabseltimeN8000 <- ggplot(stabsel_N8000,
                                aes(x=gen, y=time)) +
  geom_line()


# Only ran for 100,000 generations, distinct difference between burn-in time per gen and the selection model,
# so we make a linear model based off only the selection model run data (without the burn in)

sel8000 <- subset(stabsel_N8000, gen > 50000)

lm(time ~ gen, data = sel8000)

# 1.954x + -71957.555
# for 150,000 gens = 61.428 hours - assumes first 50,000 gens were this slow, they aren't

sel_burnin8000 <- subset(stabsel_N8000, gen <= 50000)

lm(time ~ gen, data = sel_burnin8000)

# 0.5371x + -1892.6575
# for 50,000 gens = 24962.3425 seconds, actual time 24959.9

# for total estimated time: (lm(sel8000(150,000) - lm(sel8000(50,000) ) + burnin time
# (221142.445 - 25742.445) + 24959.9 = 220349.9 = 61.21 hours







