# Genetic Architecture in SLiM
Tracking evolution in forward-modelling software SLiM 3.4: https://messerlab.org/slim/

This repository will contain Eidos and R code tracking my progress through my honours project. I aim to use SLiM to model a variety of genetic architectures and their effect on polygenic adaptation/random evolution (altering number of loci affecting a trait, their size effects, linkage, amount of deleterious mutation, extent of pleiotropy). 

Several branches exist: My original aims are below, but were changed throughout my project. The master branch represents the early models built to test these hypotheses, but were too slow to run and too difficult to sample with the number of parameters.
						Genome-restructure, fitness-function, and stabsel-opt are testing branches related to aspects of the master models
						Simplification represents a toned down version of the models (fewer parameters) and the initial analyses - using G matrices and eigentensors
						LS_D_R represents a shift in focus for analysis: mainly R code based on mean variances and covariances rather than G matrices, and a focus on additive effect size (LS), mutation rate (D), and recombination (R)

Code in the branch LS_D_R represents my final thesis, answering which architectures and mutation/selection regimes constrain adaptation, and which models promote adaptedness vs adaptability.




###############################################################################################################
AIM 1: The sampling distribution of complex neutral genetic architectures.
Using simulations of neutral populations at varying levels of genetic diversity, I will quantify differences in genetic variation and covariation as a result of their genetic architecture. I will then be able to evaluate the random sampling distribution of G for the parameter space detailed in the research approach below. This will also help me test my simulations against known theoretical results in population and evolutionary quantitative genetics.  

AIM 2: The effect of recombination on the sampling distribution of neutral complex genetic architectures.
Using forward simulation of populations and their genomes, I will determine differences in genetic variation and covariation between populations as a result of differing recombination rates. I will then be able to evaluate the effects of linkage on the sampling distribution of G for the parameter space detailed in the research approach below. 

AIM 3: The effect of recombination on the sampling distribution of complex genetic architectures after adaptation
I will introduce beneficial and deleterious mutations to simulations developed under the two null models produced in AIMs 1 and 2 to measure how genetic architecture and recombination rate affect genetic variation and covariation under selection.  I will measure how they affect adaptation walks, and maintenance of variation in natural populations. 

By realising these aims, I will clarify how the genetic architecture of complex traits affects evolution under both genetic drift and selection. Through investigating variation in genetic architectures, we can also learn what type of genetic architectures and levels of recombination rates might put populations at risk of extinction, and which architectures might facilitate persistence over evolutionary time. I will further discuss the significance of this research below, followed by a description of my expectations and methodology.
