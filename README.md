#  A comprehensive representation of genotypic selection beyond two alleles, and its estimation from time-series data

## Project description
In this repository we collected the scripts accompanying our research
project about drift and selection at multiallelic loci, and how to estimate model parameters from emprical data.
The scripts can be used by interested readers to reproduce our results, or
they can be modified to answer new research questions.

## Main scripts
We have used Matlab and Python to illustrate the dynimics of the model, for applying the model to empirical data analysis, and for providing some supplementary information. Note that in some script the matrix of fitness effects is still called by its former designation A instead of F.

### Illustration of model dynamics
The script `det_and_stoch_random_A_regions.m` uses the a matrix of fitness effects calculated by the function `calc_A_fluctuating_regions.m` (but now hard-coded for replicatability) to calculate the allele frequency trajectories. The resulting plot is Figure 1 in the manuscript.
The script `det_and_stoch_additive_regions.m` uses the a matrix of fitness effects calculated by the function `calc_A_additive.m` to calculate the allele frequency trajectories for a case where allele-specific fitness effects are combined additively to genotype-specific fitness effects. The result is Figure 2 in the manuscript.
The script `det_and_stoch_nfds_regions.m` uses the a matrix of fitness effects calculated by the function `calc_A_nfds.m` to calculate the allele frequency trajectories for a case where genotype-specific fitness effects are negatively related to the respective genotype frequencies. The resulting plot is Figure 3 in the manuscript.
These three scripts calculate determinisitc as well as stochastic trajectories, where the stochastic trajectories are reproducible because of the seed used.
In principle, it is easy to programm other cases of the selective force, e.g. heterozygote advantage, fluctuating selection etc. If you are interested in this or need help with implementing a specific selection scenario please contact Nikolas Vellnow (nikolas.vellnow@tu-dortmund.de). He is happy to help.

### Applications to empirical data analysis
#### Processing data from yeast experimental evolution
We downloaded the files `DBVPG6044_hap_freqs_10kb.txt`, `DBVPG6765_hap_freqs_10kb.txt`, `Y12_hap_freqs_10kb.txt` and `YPS128_hap_freqs_10kb.txt` with haplotype frequencies from https://doi.org/10.5061/dryad.8gtht76mz. The Python script `haplotypes_yeast_exp_evol.py` can be used to extract the trajectories for the haplotypes at a specific position, e.g. position 457847 on chromosome C16, and to save it as a `.csv`-file `hap_freqs_rep_3_chr_C16_pos_457847.csv`.
The Python script `extract_haplo_trajectories_yeast_exp_evol.py` can be used to extract haplotype frequencies for all positions across the whole yeast genome and save them in a single `.csv`-file, e.g. `haplo_trajectories_rep_3.csv`.

#### Estimation of F
The script `fit_A_to_emp_trajectory.m` uses haplotype frequency trajectory data saved in a `.csv`-file, e.g. `hap_freqs_rep_3_chr_C16_pos_457847.csv`, to estimate an F-matrix. The script calls the functions `SMatVec.m`, `CostSparseFn.m`, `convert_F_to_F_ij.m` and `convert_F_to_w.m`. Also, the "Optimization Toolbox" and the "Global Optimization Toolbox" for Matlab have to be installed to use this script. Figure 4 in the manuscript (as well as Figure S1 in the Supplementary Material) is the result of this script.
The script `scan_yeast_genome_F.m` takes a `.csv`-file with trajectory data across many genomic windows, e.g. `haplo_trajectories_rep_3.csv`, estimates F for each window and saves the information in a `.mat`-file, e.g. `F_results.mat`. The script calls the functions `fit_F_sparse_func.m`, with its dependencies. Also, the "Optimization Toolbox" and the "Global Optimization Toolbox" for Matlab have to be installed to use this script. The script `scan_yeast_genome_analysis.m` then takes the estimated F matrices, e.g. `F_results.mat`, analyses them and outputs a plot like the one in Figure 5 in the manuscript.

#### Extracting trajectories from pooled sequencing data
We downloaded `.bam`-files for different sampling time points from the DEST webpage https://dest.bio/data-files/sync-bam-bed. We had to modify the bam files so that the read groups had unique identifiers, which we did with the shell script `fix_read_groups.sh`.
The shell script `run_extraction.sh` takes these modified `.bam`-files as input and uses freebayes to call haplotypes/alleles and outputs these to a `.vcf`-file. It then filters the sites so that only sites with at least 3 alleles are included in the output `filtered_sites.tsv`. A textfile with a list of the bam files and the `.fasta`-file of the *D. melanogaster* reference genome is necessary for this script.
The Python script `parse_filtered_sites_w_hist.py` takes the textfile `filtered_sites.tsv` as input, calculates the number of sites and the number alleles at those sites, and creates a histogram. The threshold for hard filtering with the MIN_READS constant in line 27 of the code. Figure 6 in the manuscript ( and FIgure S3 in the Supplementary Material) is a plot generated with this script.


## Supplementary material
The Matlab script `MultiAlleleTestTDependenceMethods.m` simulates how the length of allele frequency trajectories affects the estimate of **F** and creates 4 plots that were used to create Figure S2 in the manuscript.
The script `detection_prob.m` calculates and plots the probability to sample an allele at least once depending on its frequency and the sequencing depth. The script `threshold_sim.m` calculates and plots the frequency trajectories of 10 alleles and the effect that an observation threshold has on them. The resulting plots of theses two scripts are Figure S4 (a) and (b), respectively.

## Manuscript
The preprint of this manuscript has been deposited on bioRxiv: https://doi.org/10.1101/2024.11.08.622587 

Information about the published article will come as soon as possible.