# Python script to handle neutral burn-in through coalescent modelling. 
# We simulate an initial ancestry and fill in mutations between
# This script expects several arguments in a particular order:
#   The modelindex
#   The seed

import msprime, pyslim
from sys import argv, path
from os.path import abspath
from pandas import read_csv
from math import ceil
# noiseseedgen.py; found in Tools/NoiseSeedRNG
#path.append(abspath('../SeedGenerator/NoiseSeedRNG'))
#import noiseseedgen

# Noise RNG function, rather than importing
def randNoise32(position, seed):
    #Uses bit manipulation to generate 32-bit pseudo-random numbers.

    #Parameters:
    #       position (int): Where in the noise function should be sampled?
    #       seed (int): Input seed for bit mashing.

    #Returns:
    #       mangled (int32): 32-bit bit-mangled value.
    PRIME1 = 0x75BD0FB
    PRIME2 = 0x75BD12D
    mangled = position
    mangled *= PRIME1
    mangled += seed
    mangled ^= (mangled >> 8)
    mangled *= PRIME2
    mangled ^= (mangled << 6)
    mangled *= mangled
    mangled ^= (mangled >> 5)
    # Cast into 32 bit value
    return mangled & 0xffffffff



# Read in parameters
combos = read_csv("$HOME/tests/CoalBurnTests/R/lhc_coalburntest.csv")
Ne = ceil(combos.Ne.get(int(argv[1]) - 1))
rwide = combos.rwide.get(int(argv[1]) - 1)

# Change the seed so we take into account both modelindex and seed
seed = randNoise32(int(argv[1]), int(argv[2]))

pop = msprime.Demography()
pop.add_population(name = "A", initial_size = Ne)

print("Constructing tree sequence...")
ts = msprime.sim_ancestry(samples = Ne, demography = pop, sequence_length = 100, discrete_genome = True, 
                          recombination_rate = rwide, random_seed = seed, model = "dtwf")
print("Tree sequence done.")

# Mutation rate is the same as in the SLiM model: we need to adjust for the
# different types of mutations that can occur by dividing the overall rate
# 8.045e-6 into the proper amount, given by K
# K = Low, each type is (8.045e-6)/3
# K = Med, each type is 31.25%    62.5%     6.25% of the mutation rate
#                       2.514e-6  5.028e-6  5.028e-7
# K = High, each type is 33.11%    66.23%    0.66%
#                        2.664e-6  5.328e-6  5.310e-8

# We use a dictionary key to choose the right value depending on the modelindex
'''
constraint = combos.K.get(int(argv[1]) - 1)

def kLow():
    return (8.045e-6/3, 8.045e-6/3, 8.045e-6/3)

def kMed():
    return(2.514e-6, 5.028e-6, 5.028e-7)

def kHigh():
    return(2.664e-6, 5.328e-6, 5.310e-8)

mutRates = {'\"c(1.0, 0.0, 0.0)\"': kLow(),
            '\"c(0.0, 1.0, 0.0)\"': kMed(),
            '\"c(0.0, 0.0, 1.0)\"': kHigh()
}

mutRates = mutRates[constraint]
'''


print("Adding mutations...")
mutated = pyslim.annotate_defaults(ts, model_type="WF", slim_generation=1)
mutated = msprime.sim_mutations(mutated, rate = (8.045e-6)/2, random_seed = seed, 
                                model = msprime.SLiMMutationModel(type = 1))
mutated.simplify()
filepath = "treeburnin_{seed}_{modelindex}.trees".format(seed = str(argv[2]), modelindex = str(argv[1]))
print(f"Dumping tree to ./{filepath}")
mutated.dump(filepath)
quit()