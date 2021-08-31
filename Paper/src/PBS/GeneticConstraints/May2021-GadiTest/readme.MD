## May 2021 - Gadi Test
This test uses Latin Hypercube Sampling to generate a large combination of parameters to test the Gadi HPC.
Rather than a traditional PBS job array, I'm using an nci-parallel job, similar to Embedded Nimrod at UQ Tinaroo

Here, I try to run on 20,736 cores, running 48 seeds for 1024 LHC samples, for 3 levels of constraint
This totals to $(1024*3)*48 = 147456$ runs 

The LHC includes Ne, rwide, locisigma, nloci, fitness function width (selection strength), and fitness optimum position.

My home PC (i7 4770) runs worst-case scenario combinations (Ne = 20000; rwide = 0.5) in around 8.3 hours, for 100,000 total generations
(50k burn-in, 25k each of initial stabilising and post-environmental shift stabilising selection)

$8.3 * 147456 = 1223884.8$ hrs
/20736 = ~59 hours + remainder of 8.3 hours for all sims to finish, 67.3 hours

This should be faster on Gadi, as it has newer processors with higher IPC and clock speed, but how much faster is uncertain as of now.

Gadi's queues typically have walltimes of under 48 hours, and using 20736 cores reduces this to only 5 hours. 
Hence, we may have to split this job into subjobs, using 2976 cores at a time (as this number allows a maximum wall time of 10 hours)
Alternatively, we could get permissions to run 20736 cores for a little over 5 hours, and split into subjobs that way: we would need to do 8 runs
to complete the job this way, but we could automate that.