# coarse-graining-phylogenetics

prerequisite
1. this script requires 'perl data language' to function, see "http://pdl.perl.org/?page=install" on installation method
2. the perl script pro_CGP_test.pl depends on the accompanying *.pm files
3. the script has 5 arguments: (a) file_strain_pair_names, (b) file_strain_pair_segment_SSPs, (c) segment size, (d) minimum number of simulation steps, (e) number of steps to print tree/parameters once
4. example code:

    perl pro_CGP_test.pl data_segmentSize30.pair_name data_segmentSize30.pair_SNP 30 1000000 1000 > output.tree 2> output.parameters


We reconstructed the tree of 10 strains, which has 45 strain pairs. The 45 lines in file ‘data_segmentSize30.pair_name’ contains the names (NCBI IDs) of the 45 pairs; each of the 45 arrays in ‘data_segmentSize30.pair_SNP’ contains the number of SSPs in every segment. The segment size in the data file provided is 30, we wanted the simulation to last for 1,000,000 steps and print out the results every 1,000 steps; note that the stopping criteria of the algorithm is to have the logarithm of the posterior probability not increase by more than 1 in 200,000 steps, and hence the script will last for 200,000 steps even if we set minimum-number-of-simulation-steps equal 1. 

The first output file ‘output.tree’ contains the the posterior trees; the second output file ‘output.parameters’ contains the posterior parameters; ‘output.parameters’ has 7 columns: (1) step ID, (2) logarithm of posterior probability, (3) μ, (4) ρ, (5) θ, (6) δTE, (7) time-unit of each step.




