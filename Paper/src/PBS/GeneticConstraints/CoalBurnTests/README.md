# Coalescent burn-in tests

The new versions of `msprime` have support for SLiM mutations, making a coalescent burn-in approach possible. This seems very promising, being orders of magnitude faster than a forward simulation.
The downside is we're only burning in neutral variation, without the context of background selection, however we can do a hybrid approach, and have a few generations of forward burn-in to make sure everything is alright (i.e. heterozygosity is in expected bounds).

In these tests I:
- Verify that the coalescent approach generates levels of diversity similar to that of a forward-approach
- Time the process to see what kind of speed increase there is
- Test a hybrid approach to see if effects on heterozygosity are different when background selection is used in the forward stage of the burn-in

For the coalescent model, it is subject to a further 2500 generations of burn-in to smooth over differences between WF and the backwards time WF burn-in.

The control model is subject to 52500 generations of burn in, 50000 to compare to the coalescent, and 2500 to match the difference smoothing.

There is then 1000 generations of stabilising selection for both models followed by a shift and another 5000 

I vary recombination rate and Ne to verify the models behave over these parameters
640 total models: 4 nodes, 24 cores means each core does ~7 simulations