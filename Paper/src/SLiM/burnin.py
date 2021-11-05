# Python script to handle neutral burn-in through coalescent modelling. 
# A SLiM model constructs the initial .tree file, which is then populated
# with random mutations. 
# This script expects several arguments in a particular order:
#   The .tree filepath
#   The modelindex
#   The seed

import msprime, pyslim
from sys import argv
from pandas import read_csv

# SLiM manual 3.6 pp 415 section 17.2
ts = pyslim.load(argv[0])
ts = ts.simplify()

combos = read_csv("./lhc_combos.csv")
rwide = combos.rwide.get(argv[1] - 1)

# Mutation rate is the same as in the SLiM model: we need to adjust for the
# different types of mutations that can occur by dividing the overall rate
# 8.045e-6 into the proper amount, given by K
# K = Low, each type is (8.045e-6)/3
# K = Med, each type is 31.25%    62.5%     6.25% of the mutation rate
#                       2.514e-6  5.028e-6  5.028e-7
# K = High, each type is 33.11%    66.23%    0.66%
#                        2.664e-6  5.328e-6  5.310e-8


# We use a dictionary key to choose the right value depending on the modelindex
constraint = combos.K.get(argv[1] - 1)

def kLow():
    return (8.045e-6/3, 8.045e-6/3, 8.045e-6/3)

def kMed():
    return(2.514e-6, 5.028e-6, 5.028e-7)

def kHigh():
    return(2.664e-6, 5.328e-6, 5.310e-8)

mutRates = {"c(1.0, 0.0, 0.0)": kLow(),
            "c(0.0, 1.0, 0.0)": kMed(),
            "c(0.0, 0.0, 1.0)": kHigh()
}

mutRates = mutRates[constraint]

# Should we change the seed? Will have to see how random it is with different rates
# Will probably be the same with kLow(), so might have to change this

mutated = msprime.sim_mutations(ts, rate = mutRates[0], random_seed = argv[2],
                                model = msprime.SLiMMutationModel(type = 1))

mutated = msprime.sim_mutations(mutated, rate = mutRates[1], random_seed = argv[2],
                                model = msprime.SLiMMutationModel(type = 2))

mutated = msprime.sim_mutations(ts, rate = mutRates[2], random_seed = argv[2],
                                model = msprime.SLiMMutationModel(type = 3))
mutated.dump(argv[0])