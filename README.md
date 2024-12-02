#  A comprehensive representation of selection at loci with multiple alleles that allows complex forms of genotypic fitness

## Project description
In this repository we collected the scripts accompanying our research
project about drift and selection at multiallelic loci.
The scripts can be used by interested readers to reproduce our results, or
they can be modified for new research questions.

## Main scripts
The main scripts of interest are: `det_and_stoch_random_A_regions.m`, `det_and_stoch_additive_regions.m`, `det_and_stoch_nfds_regions.m`, `detection_prob.m`, and `threshold_sim.m`. These scripts were used in the manuscript.

In particular, `det_and_stoch_random_A_regions.m` uses the a matrix of fitness effects calculated by the function `calc_A_fluctuating_regions.m` to calculate the allele frequency trajectories. The script `det_and_stoch_additive_regions.m` uses the a matrix of fitness effects calculated by the function `calc_A_additive.m` to calculate the allele frequency trajectories for a case where allele-specific fitness effects are combined additively to genotype-specific fitness effects. The script `det_and_stoch_nfds.m` uses the a matrix of fitness effects calculated by the function `calc_A_nfds.m` to calculate the allele frequency trajectories for a case where genotype-specific fitness effects are negatively related to the respective genotype frequencies.
These three scripts calculate determinisitc as well as stochastic trajectories, where the stochastic trajectories are reproducible because of the seed used.

The script `detection_prob.m` calculates and plots the probability to sample an allele at least once depending on its frequency and the sequencing depth.
The script `threshold_sim.m` calculates and plots the frequency trajectories of 10 alleles and the effect that an observation threshold has on them.

## Manuscript
The preprint of this manuscript has been deposited on bioRxiv: https://doi.org/10.1101/2024.11.08.622587 

Information about the published article will come as soon as possible.