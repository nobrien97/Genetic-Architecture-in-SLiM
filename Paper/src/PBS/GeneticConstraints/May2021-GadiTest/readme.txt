This test uses Latin Hypercube Sampling to generate a large combination of parameters to test the Gadi HPC.
Rather than a traditional PBS job array, I'm using an nci-parallel job, similar to Embedded Nimrod at UQ Tinaroo

Here, I try to run on 10,000 cores, running 48 seeds for 1024 LHC samples, for 3 levels of constraint
This totals to (1024*3)*48 = 147456 runs 
