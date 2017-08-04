#!/usr/bin/perl -w

use strict;
use warnings;

use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

use POSIX;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Storable qw/dclone/;
use CGP_lib_tree_operation qw/f_extract_clade_distMatrix_hash f_read_file_pairName_pairSNP f_cut_and_graft f_read_newick_file/;
use CGP_lib_SSP_theory qw/get_parameter_probability2 f_new_dist_map/;
use IO::Handle;



# perl script_CGP3_MCMC.pl tmpoutput_pairName.dat tmpoutput_pairSSP.dat 10 1000 outputFileHeader
my ($fileWithName,$fileNSNP,$NSNPperSegment,$NMCstep_interval,$outputFileHeader,$N_step_relaxation,$file_newick_tree) = @ARGV;
$| = 1;
my $maxNSNP = $NSNPperSegment;
my $mcmcTemperature = 10;	# dummy temperature
my $minMCMCstep = 1; # minimum number of MCMC steps
my $maxdN = $NSNPperSegment;

my ($mu) = (0.02);
my ($rho,$deltaTE) = (0.002,3);
my $tmp_resolution = 20;
my $dt = 0.1;
$dt = ($dt>1/$mu/$tmp_resolution) ? 1/$mu/$tmp_resolution : $dt;
if ($rho>0) {
	$dt = ($dt>1/$rho/$tmp_resolution) ? 1/$rho/$tmp_resolution : $dt;
}
my $N = 5000;
my $theta = min(2*$mu*$N,$maxNSNP/2);

$N_step_relaxation = 200000 unless $N_step_relaxation;





# read the SSP distribution
# inputs: 
#		$fileWithName: file with names of strain pairs
#		$fileNSNP: file with number of SSP of the pairs
#		$maxNSNP: maximum number of SSP per segment
# outputs:
#		$arr_strain: array of strain names
#		$arr_strainPair: array of strain_pairs
#		$arr_data_SNP: matrix storing number of SSPs at different segments / genome pairs
#		$arr_data_SNP_dist: array of SSP distributions
#		$hash_strainPair2SNPdist: hash mapping strain pairs to distributions
#		$maxDelta: maximum number of SSP in a segment found in the data
#		$dummy1: dummy variable
#		$max_ddelta: maximum divergence between any pair of genomes
# the file is read 3 times to get an optimal number of array size ("l_s^cutoff")
my ($arr_strain,$arr_strainPair,$arr_data_SNP,$arr_data_SNP_dist,$hash_strainPair2SNPdist,$maxDelta,$dummy1,$max_ddelta) = @{f_read_file_pairName_pairSNP($fileWithName,$fileNSNP,$maxNSNP)};
$maxNSNP = $max_ddelta+5 if $maxNSNP>$max_ddelta+5;
($arr_strain,$arr_strainPair,$arr_data_SNP,$arr_data_SNP_dist,$hash_strainPair2SNPdist,$maxDelta,$dummy1,$max_ddelta) = @{f_read_file_pairName_pairSNP($fileWithName,$fileNSNP,$maxNSNP)};
if ($maxNSNP>100) {
	$maxNSNP = 100;
	($arr_strain,$arr_strainPair,$arr_data_SNP,$arr_data_SNP_dist,$hash_strainPair2SNPdist,$maxDelta,$dummy1,$max_ddelta) = @{f_read_file_pairName_pairSNP($fileWithName,$fileNSNP,$maxNSNP)};
}


# if a tree file is given, read the tree file, and convert it into a distance matrix; the topology of the tree in the mcmc will not change if a initial tree is given
my $mutateTreeTopology = 1;
my $initial_distance_matrix;
if ($file_newick_tree && (-e $file_newick_tree)) {
	my $newick_txt = f_read_newick_file($file_newick_tree);
	$initial_distance_matrix = f_extract_clade_distMatrix_hash($newick_txt);
	$mutateTreeTopology = 0;
}


# perform MCMC
# inputs:
#		$arr_data_SNP_dist: SSP distributions of different strain pairs
#		$minMCMCstep: minimum number of MCMC step, set to 1
#		$NMCstep_interval: number of interval step to print data
#		$endPairDiv: terminate the calculation when a pair of strains has reach divergence $endPairDiv
#		$endPairDiv: array size, same as l_s^cutoff in the model
#		$dt: time elapse per step
#		$maxNSNP: maximum number of SSP in a segment observed in data, used to set the array size
#		$arr_strain: array of strain names
#		$arr_strainPair: array of strain pairs
#		$mcmcTemperature: temperature of Monte Carlo simulation
#		$step0: ID of first simulation step, set to 0
#		$hash_strainPair2SNPdist: hash mapping strain pairs to their SSP distributions
my $endPairDiv = 1.3*$maxDelta;
my ($fout_tree,$fout_score);
open $fout_tree,">".$outputFileHeader.".record_tree" or die "can't write file $outputFileHeader\.record_tree\n";
open $fout_score,">".$outputFileHeader.".record_score" or die "can't write file record_score\n";
evo_mcmc5($arr_data_SNP_dist,$minMCMCstep,$NMCstep_interval,$mu,$rho,$theta,$deltaTE,$endPairDiv,$dt,$maxNSNP,$initial_distance_matrix,$arr_strain,$arr_strainPair,$mcmcTemperature,time(),0,$hash_strainPair2SNPdist,$mutateTreeTopology,$fout_tree,$fout_score);


exit 1;















sub evo_mcmc5 {
	my ($arr_data_SNP_dist,$minMCMCstep,$NMCstep_interval,$mu,$rho,$theta,$deltaTE,$endPairDiv,$dt,$matrixSize,$start_dist_matrix_hash_input,$arr_strains,$arrLine2StrainPairs,$MCtemperature,$curtime,$step0,$hash_strainPair2SNPdist,$mutateTreeTopology,$fout_tree,$fout_score) = @_;
	my ($maxx_score,$maxx_score_step) = (-1e100,0);
	my ($sav_score,$sav_score_step) = (-1e100,0);
	my $step_to_stop = $minMCMCstep;
	my $maxxTree;
	my @maxxPara;
	my $savTree;
	my @savPara;

	# generate the model
	# input of the function: model parameters (mu,rho,theta,deltaTE,arraySize,cutOffTimeStep)
	# output:
	#		$strainpair2Prob_hash: cross entropies for different coalescent time of strain pairs
	#			$strainpair2Prob_hash->{$strain1}{$strain2}[$timeStep] = $negativeCrossEntropy
	#		$arrTime: list of time at different steps
	my $f_genModel = sub {
		my ($arr_data_SNP_dist,$mu,$rho,$theta,$deltaTE,$endPairDiv,$dt,$matrixSize,$max_time_step) = @_;
		my ($strainpair2Prob_array,$arrTime);
		my $tmp = get_parameter_probability2($arr_data_SNP_dist,$mu,$rho,$theta,$deltaTE,$endPairDiv,$dt,$matrixSize,1,$max_time_step);
		($strainpair2Prob_array,$arrTime) = @{$tmp};
		my $strainpair2Prob_hash = {};
		for my $q (0..$#$arrLine2StrainPairs) {
			my ($s1,$s2) = @{$arrLine2StrainPairs->[$q]};
			unless (scalar(@{$strainpair2Prob_array->[$q]})) {
				print STDERR "eerrrror $mu,$rho,$theta,$deltaTE,$endPairDiv,$dt,$matrixSize\t$s1,$s2\n";
			}
			$strainpair2Prob_hash->{$s1}{$s2} = $strainpair2Prob_array->[$q];
			$strainpair2Prob_hash->{$s2}{$s1} = $strainpair2Prob_array->[$q];
		}
		return [$strainpair2Prob_hash,$arrTime];
	};	
	
	# generate the coalescent time matrix of the strain pairs
	# input:
	#		$input_strainpair2Prob_hash: cross entropies for different coalescent time of strain pairs, predicted from certain (mu,rho,theta,deltaTE)
	# output:
	#		$tmp_start_dist_matrix_hash: estimated coalescent time matrix of the strain pairs
	#		$maxTime: maximum coalescent time among the strain pairs
	my $f_initialize_time_matrix = sub {
		my ($input_strainpair2Prob_hash) = @_;
		my $tmp_start_dist_matrix_hash;
		my $maxTime = 0;
		for my $q (0..$#$arrLine2StrainPairs) {
			my ($s1,$s2) = @{$arrLine2StrainPairs->[$q]};
			my @tmpProbTime = (-1e50,-1);
			for my $r (0..$#{$input_strainpair2Prob_hash->{$s1}{$s2}}) {
				if ($input_strainpair2Prob_hash->{$s1}{$s2}[$r] > $tmpProbTime[0]) {
					@tmpProbTime = ($input_strainpair2Prob_hash->{$s1}{$s2}[$r],$r);
				}
			}
			$tmp_start_dist_matrix_hash->{$s1}{$s2} = $tmpProbTime[1];
			$tmp_start_dist_matrix_hash->{$s2}{$s1} = $tmpProbTime[1];
			if ($maxTime<$tmpProbTime[1]) {
				$maxTime = $tmpProbTime[1];
			}
		}
		return [$tmp_start_dist_matrix_hash,$maxTime];
	};	
	
	# generate the matrix (cross-entropy-at-different-coalescent-time vs strain-pairs)
	my ($strainpair2Prob_hash,$arrTime) = @{$f_genModel->($arr_data_SNP_dist,$mu,$rho,$theta,$deltaTE,$endPairDiv,$dt,$matrixSize)};
	
	# resacle my and rho
	# calculate the SSP distributions from the intial mu, rho, theta, deltaTE
	my ($tmp_dist_matrix,$maxTime) = @{$f_initialize_time_matrix->($strainpair2Prob_hash,$arrTime)};
	if ($maxTime>200) {
		my $rratio = $maxTime/200;
		if ($mu*$rratio>0.02) {
			$rratio = 0.02/$mu;
		}
		$mu *= $rratio;
		$rho *= $rratio;
		($strainpair2Prob_hash,$arrTime) = @{$f_genModel->($arr_data_SNP_dist,$mu,$rho,$theta,$deltaTE,$endPairDiv,$dt,$matrixSize)};			
	}
	my $start_dist_matrix_hash;
	unless ($start_dist_matrix_hash_input) {
		# if an initial distance-matrix-of-strain-pairs ($start_dist_matrix_hash) is not given, then fill this matrix up with most-likely-coalescent-time estimated from initial model parameters
		($tmp_dist_matrix,$maxTime) = @{$f_initialize_time_matrix->($strainpair2Prob_hash,$arrTime)};
		$start_dist_matrix_hash = $tmp_dist_matrix;
	} else {
		# if an initial tree is given, then rescale the initial tree
		my $tmpmaxdist = 0;
		for my $q (0..$#$arr_strains-1) {			
			my $sQ = $arr_strains->[$q];
			die "$sQ is missing in the initial tree" unless exists $start_dist_matrix_hash_input->{$sQ};
			for my $r ($q+1..$#$arr_strains) {
				my $sR = $arr_strains->[$r];
				die "$sR is missing in the initial tree" unless exists $start_dist_matrix_hash_input->{$sQ}{$sR};
				$tmpmaxdist = $start_dist_matrix_hash_input->{$sQ}{$sR} if $start_dist_matrix_hash_input->{$sQ}{$sR} > $tmpmaxdist;
			}
		}		
		my $tmpratio = 200/$tmpmaxdist;
		$maxTime = 0;
		for my $q (0..$#$arr_strains-1) {
			my $sQ = $arr_strains->[$q];
			die "$sQ is missing in the initial tree" unless exists $start_dist_matrix_hash_input->{$sQ};
			for my $r ($q+1..$#$arr_strains) {
				my $sR = $arr_strains->[$r];
				die "$sR is missing in the initial tree" unless exists $start_dist_matrix_hash_input->{$sQ}{$sR};
				my $curTime = floor($start_dist_matrix_hash_input->{$sQ}{$sR} * $tmpratio)+1;				
				$start_dist_matrix_hash->{$sQ}{$sR} = $curTime;
				$start_dist_matrix_hash->{$sR}{$sQ} = $curTime;
				$maxTime = $curTime if $maxTime < $curTime;
			}
		}
		($strainpair2Prob_hash,$arrTime) = @{$f_genModel->($arr_data_SNP_dist,$mu,$rho,$theta,$deltaTE,$endPairDiv,$dt,$matrixSize,$maxTime+200)};
	}
	
	# convert the strain-pair-distance-matrix into a phylogenetic-tree-object, with tree topology inferred from single linkaged clustering
	# input:
	#		$start_dist_matrix_hash: strian pair distance matrix
	# output:
	#		$clusterData: object of the ultrametric phylogeneic tree
	#		$distmap: strain-pair-distance-matrix modified from $start_dist_matrix_hash, which is equivalent to a single linkaged clustering tree
	#		$maxTreeHeight: root node age
	my ($clusterData,$distmap,$maxTreeHeight) = @{distmatrix2cluster($start_dist_matrix_hash)};
	if ($maxTreeHeight>$#$arrTime) {
		($strainpair2Prob_hash,$arrTime) = @{$f_genModel->($arr_data_SNP_dist,$mu,$rho,$theta,$deltaTE,$endPairDiv,$dt,$matrixSize,$maxTreeHeight*2)};
		($clusterData,$distmap,$maxTreeHeight) = @{distmatrix2cluster($start_dist_matrix_hash)};
	}
	
	# $absoluteTreeHeight is used as the termination coalescent-time for the calculation of matrix (cross-entropy-at-different-coalescent-time vs strain-pairs)
	my $absoluteTreeHeight = max(int($maxTreeHeight*1.2),$#$arrTime);
	
	# $treeProb is the log-posterior-likelihood calculated from the tree (stored in $distmap and $clusterData) and cross-entropies of different strain pairs (stored in $strainpair2Prob_hash)
	my $treeProb = get_tree_probability2($strainpair2Prob_hash,$distmap);
	$MCtemperature = abs($treeProb/10000);

	# MCMC step counter
	my $curStep = -1;
	
	# arrays storing the parameters of possible moves in an MCMC step
	my @changeProb_IncreaseProb;
	my @changeProb_DecreaseProb;
	my @changeData_IncreaseProb;
	my @changeData_DecreaseProb;
	
	# temporary function that decides which move to execute in an MCMC step
	my $f_update_tree_topology = sub {
		if (scalar(@changeProb_IncreaseProb)>0) {					
			my $ii = choose_index_from_prob_array(@changeProb_IncreaseProb);
			($clusterData,$distmap,$maxTreeHeight,$treeProb) = @{$changeData_IncreaseProb[$ii]};
		} elsif (scalar(@changeProb_DecreaseProb)>0) {
			push @changeProb_DecreaseProb,1;
			my $ii = choose_index_from_prob_array(@changeProb_DecreaseProb);
			if (defined($ii) && ($ii<$#changeProb_DecreaseProb)) {
				($clusterData,$distmap,$maxTreeHeight,$treeProb) = @{$changeData_DecreaseProb[$ii]};
			}
		}

		$#changeProb_IncreaseProb = -1;
		$#changeProb_DecreaseProb = -1;
		$#changeData_IncreaseProb = -1;
		$#changeData_DecreaseProb = -1;
	};
	
	# start mcmc simulation
	while ($curStep<$step_to_stop) {
		# print and write the model parameters and tree
		if ($curStep%$NMCstep_interval==0) {
			my $printstep = $step0 + $curStep;
print "step=$printstep\tlikelihood=$treeProb\tmu=$mu\trho=$rho\ttheta=$theta\tdeltaTE=$deltaTE\tdt=$dt\n";
print gen_newick($clusterData).";\n";
			print $fout_score "$printstep\t$treeProb\t$mu\t$rho\t$theta\t$deltaTE\t$dt\n";
			print $fout_tree gen_newick($clusterData).";\n";
			$fout_score->autoflush;
			$fout_tree->autoflush;
		}	
		
		# the following move decide whether to change the model parameters or not
		if (rand()<1/(scalar @$arrLine2StrainPairs)) {
			$curStep++;
						
			my @tmpParameter;
			my @moves; # store the possible ways to change the model parameters
			$moves[0] = sub {
				my $rb = $rho*(1+0.4*exp(-1/(1-rand())));
				push @tmpParameter,[$rb,$theta,$deltaTE] if $rb<1;
			};
			$moves[1] = sub {
				push @tmpParameter,[$rho/(1+0.4*exp(-1/(1-rand()))),$theta,$deltaTE];
			};
			$moves[2] = sub {
				my $thetab = $theta*(1+0.4*exp(-1/(1-rand())));
				push @tmpParameter,[$rho,$thetab,$deltaTE] if ($thetab<$maxdN);
			};
			$moves[3] = sub {
				push @tmpParameter,[$rho,$theta/(1+0.4*exp(-1/(1-rand()))),$deltaTE];
			};
			$moves[4] = sub {
				my $deltaTEb = $deltaTE*(1+0.4*exp(-1/(1-rand())));
				push @tmpParameter,[$rho,$theta,$deltaTEb] if $deltaTEb<$maxdN;
			};
			$moves[5] = sub {
				push @tmpParameter,[$rho,$theta,$deltaTE/(1+0.4*exp(-1/(1-rand())))];
			};
			
			#choose one move and change the parameters
			@moves = shuffle @moves;
			for my $q (0..$#moves) {
				$moves[$q]->();
				last if @tmpParameter>0;
			}
			
			my @changeDecreaseProb;
			for ($tmpParameter[0]) {	##### @tmpParameter has 1 element
				next if $_->[0]>1.1 || $_->[1]>1000;
				my ($tmpStrainpair2Prob_hash,$tmpArrTime) = @{$f_genModel->($arr_data_SNP_dist,$mu,$_->[0],$_->[1],$_->[2],$endPairDiv,$dt,$matrixSize,int($maxTreeHeight*1.5))};
				my $newProb = get_tree_probability2($tmpStrainpair2Prob_hash,$distmap);
				if ($newProb>$treeProb) {
					($strainpair2Prob_hash,$arrTime,$treeProb) = ($tmpStrainpair2Prob_hash,$tmpArrTime,$newProb);
					($rho,$theta,$deltaTE) = @$_;
					last;
				} else {
					my @pprob = map {$changeDecreaseProb[$_]->[5]} 0..$#changeDecreaseProb;
					push @pprob,1;
					my $ii = choose_index_from_prob_array(@pprob);
					if (defined($ii) && ($ii<$#pprob)) {
						my $dummy;
						($rho,$theta,$deltaTE,$strainpair2Prob_hash,$arrTime,$dummy,$treeProb) = @{$changeDecreaseProb[$ii]}[0..6];
					}
				}
			}
		}
		
		# the following move decide whether to cut a branch and graft it to somewhere
		if ((rand()<1/(scalar @$arrLine2StrainPairs)) && $mutateTreeTopology) {
			$curStep++;
			my @clusterToMove = shuffle keys %$clusterData;
			foreach my $curcluster (@clusterToMove) {	#####
				my $loopCompleted = 0;
				my @tmpDistMap = @{f_cut_and_graft($clusterData,$curcluster,$distmap)};
				for my $cur_map (@tmpDistMap) {	#####					
					my ($tmpclusterData,$tmpMaxTreeHeight);
					($tmpclusterData,$cur_map,$tmpMaxTreeHeight) = @{distmatrix2cluster($cur_map)};
					my $newprob = get_tree_probability2($strainpair2Prob_hash,$cur_map);
					if ($newprob>$treeProb) {
						push @changeProb_IncreaseProb,exp(($newprob-$treeProb)/$MCtemperature);
						push @changeData_IncreaseProb,[$tmpclusterData,$cur_map,$tmpMaxTreeHeight,$newprob];
						$loopCompleted = 1;
					} else {
						push @changeProb_DecreaseProb,exp(($newprob-$treeProb)/$MCtemperature);
						push @changeData_DecreaseProb,[$tmpclusterData,$cur_map,$tmpMaxTreeHeight,$newprob];
						$loopCompleted = 1;
					}
					last if $loopCompleted;
				}
				last if $loopCompleted;
			}
			$f_update_tree_topology->();
		}
		
		# the follow step decide whether or not to move an internal upwards for downwards		
		do {
			# pick on internal node to move or not, the chosen internal node has index $qToMv
			my @clusterToMove = shuffle grep {scalar(keys %{$clusterData->{$_}{tipnode}})>1} keys %$clusterData;
			my @clusterWeight = (1) x scalar(@clusterToMove);			
			my $qToMv = choose_index_from_prob_array(@clusterWeight);
			
			# pick one move on one of the chosen internal node
			for my $cluster ($clusterToMove[$qToMv]) {	#####
				$curStep++;					
				my $tmpdist = $clusterData->{$cluster}{height};
				
				# save all possible moves in the @moves array
				my @moves;
				
				# if a tree is not given when algorithm starts, then we can mutate the tree
				if ($mutateTreeTopology) {
					# consider moving a node downward by 1 step, if the height of the node is greater >0
					if (($tmpdist-1>=0) && (scalar(@{$clusterData->{$cluster}{clusterTipPair}})==1)) {
						my $f = sub {
							my $tmpDistMataa = f_new_dist_map($distmap,$tmpdist-1,$clusterData->{$cluster}{clusterTipPair}[0]);
							my ($tmpclusterData,$tmpMaxTreeHeight);
							($tmpclusterData,$tmpDistMataa,$tmpMaxTreeHeight) = @{distmatrix2cluster($tmpDistMataa)};
							my $newprob = get_tree_probability2($strainpair2Prob_hash,$tmpDistMataa);
							if ($newprob>$treeProb) {
								push @changeProb_IncreaseProb,exp(($newprob-$treeProb)/$MCtemperature);
								push @changeData_IncreaseProb,[$tmpclusterData,$tmpDistMataa,$tmpMaxTreeHeight,$newprob];
							} else {
								push @changeProb_DecreaseProb,exp(($newprob-$treeProb)/$MCtemperature);
								push @changeData_DecreaseProb,[$tmpclusterData,$tmpDistMataa,$tmpMaxTreeHeight,$newprob];
							}
						};
						push @moves,$f;
					}
						
					# consider moving a node upward by 1 step, if the height of the node is greater >0
					if (($tmpdist+1<=$#$arrTime) && (scalar(@{$clusterData->{$cluster}{clusterTipPair}})==1)) {
						my $f1 = sub {
							my $tmpDistMataa = f_new_dist_map($distmap,$tmpdist+1,$clusterData->{$cluster}{clusterTipPair}[0]);
							my ($tmpclusterData,$tmpMaxTreeHeight);
							($tmpclusterData,$tmpDistMataa,$tmpMaxTreeHeight) = @{distmatrix2cluster($tmpDistMataa)};
							my $newprob = get_tree_probability2($strainpair2Prob_hash,$tmpDistMataa);
							if ($newprob>$treeProb) {
								push @changeProb_IncreaseProb,exp(($newprob-$treeProb)/$MCtemperature);
								push @changeData_IncreaseProb,[$tmpclusterData,$tmpDistMataa,$tmpMaxTreeHeight,$newprob];
							} else {
								push @changeProb_DecreaseProb,exp(($newprob-$treeProb)/$MCtemperature);
								push @changeData_DecreaseProb,[$tmpclusterData,$tmpDistMataa,$tmpMaxTreeHeight,$newprob];
							}
						};
						push @moves,$f1;
					}
				
					# consider the case of breaking down an internal node with >2 child-nodes into two, and move one downwards, if the current node height is >0
					if (($tmpdist-1>=0) && (scalar(@{$clusterData->{$cluster}{clusterTipPair}})>1)) {
						my $f = sub {
							my @tmp = shuffle 0..$#{$clusterData->{$cluster}{clusterTipPair}};		
							for my $q ($tmp[0]) {	#####			
								my @arrTipPair = @{$clusterData->{$cluster}{clusterTipPair}[$q]};
								my $tmpDistMat = f_new_dist_map($distmap,$tmpdist-1,\@arrTipPair);
								my ($tmpclusterData,$tmpMaxTreeHeight);
								($tmpclusterData,$tmpDistMat,$tmpMaxTreeHeight) = @{distmatrix2cluster($tmpDistMat)};
								my $newprob = get_tree_probability2($strainpair2Prob_hash,$tmpDistMat);
								if ($newprob>$treeProb) {
									push @changeProb_IncreaseProb,exp(($newprob-$treeProb)/$MCtemperature);
									push @changeData_IncreaseProb,[$tmpclusterData,$tmpDistMat,$tmpMaxTreeHeight,$newprob];
								} else {
									push @changeProb_DecreaseProb,exp(($newprob-$treeProb)/$MCtemperature);
									push @changeData_DecreaseProb,[$tmpclusterData,$tmpDistMat,$tmpMaxTreeHeight,$newprob];
								}
							}
						};
						push @moves,$f;
					}
					
					# consider the case of breaking down an internal node with >2 child-nodes into two, and move one upwards
					if (($tmpdist+1<=$#$arrTime) && (scalar(keys %{$clusterData->{$cluster}{downstreamCluster}})>2)) {
						my $f2 = sub {
							my @dsnodes = shuffle keys %{$clusterData->{$cluster}{downstreamCluster}};
							for my $q (0) {
								my $branchToMv = $dsnodes[$q];
								my @branchesNotMv = @dsnodes;
								splice @branchesNotMv,$q,1;
								my @arrTipPair;
								for my $tipA (keys %{$clusterData->{$branchToMv}{tipnode}}) {
									for my $branchNotMv (@branchesNotMv) {
										for my $tipB (keys %{$clusterData->{$branchNotMv}{tipnode}}) {
											push @arrTipPair,[$tipA,$tipB];
										}
									}							
								}
								my $tmpDistMat = f_new_dist_map($distmap,$tmpdist+1,\@arrTipPair);
								my ($tmpclusterData,$tmpMaxTreeHeight);
								($tmpclusterData,$tmpDistMat,$tmpMaxTreeHeight) = @{distmatrix2cluster($tmpDistMat)};		
								my $newprob = get_tree_probability2($strainpair2Prob_hash,$tmpDistMat);
								if ($newprob>$treeProb) {
									push @changeProb_IncreaseProb,exp(($newprob-$treeProb)/$MCtemperature);
									push @changeData_IncreaseProb,[$tmpclusterData,$tmpDistMat,$tmpMaxTreeHeight,$newprob];
								} else {
									push @changeProb_DecreaseProb,exp(($newprob-$treeProb)/$MCtemperature);
									push @changeData_DecreaseProb,[$tmpclusterData,$tmpDistMat,$tmpMaxTreeHeight,$newprob];
								}
							}
						};
						push @moves,$f2;
					}
				} 
				
				do { 
					my @arrAllTipPairs;
					for my $q (0..$#{$clusterData->{$cluster}{clusterTipPair}}) {
						my @arrTipPair = @{$clusterData->{$cluster}{clusterTipPair}[$q]};
						push @arrAllTipPairs,@arrTipPair if ($tmpdist+1 <= $#$arrTime);
					}
					
					# consider moving a node downward by 1 step, if the height of the node is greater >0
					if ($tmpdist-1>=0) {
						my $f = sub {
							my $canDo = 1;
							my @down_nodes = keys %{$clusterData->{$cluster}{downstreamCluster}};
							for my $dn (@down_nodes) {
								if ($tmpdist-1 <= $clusterData->{$dn}{height}) {
									$canDo = 0;
									last;
								}
							}									
							if ($canDo) {
								my $tmpDistMataa = f_new_dist_map($distmap,$tmpdist-1,\@arrAllTipPairs);
								my ($tmpclusterData,$tmpMaxTreeHeight);
								($tmpclusterData,$tmpDistMataa,$tmpMaxTreeHeight) = @{distmatrix2cluster($tmpDistMataa)};
								my $newprob = get_tree_probability2($strainpair2Prob_hash,$tmpDistMataa);
								if ($newprob>$treeProb) {
									push @changeProb_IncreaseProb,exp(($newprob-$treeProb)/$MCtemperature);
									push @changeData_IncreaseProb,[$tmpclusterData,$tmpDistMataa,$tmpMaxTreeHeight,$newprob];
								} else {
									push @changeProb_DecreaseProb,exp(($newprob-$treeProb)/$MCtemperature);
									push @changeData_DecreaseProb,[$tmpclusterData,$tmpDistMataa,$tmpMaxTreeHeight,$newprob];
								}
							}
						};
						push @moves,$f;
					}
						
					# consider moving a node upward by 1 step, if the height of the node is greater >0
					if ($tmpdist+1<=$#$arrTime) {
						my $f1 = sub {
							my $canDo = 1;
							if (exists $clusterData->{$cluster}{upstreamCluster}) {
								my $up_node = $clusterData->{$cluster}{upstreamCluster};
								$canDo = 0 if $tmpdist+1 >= $clusterData->{$up_node}{height};
							}
							if ($canDo) {
								my $tmpDistMataa = f_new_dist_map($distmap,$tmpdist+1,\@arrAllTipPairs);
								my ($tmpclusterData,$tmpMaxTreeHeight);
								($tmpclusterData,$tmpDistMataa,$tmpMaxTreeHeight) = @{distmatrix2cluster($tmpDistMataa)};
								my $newprob = get_tree_probability2($strainpair2Prob_hash,$tmpDistMataa);
								if ($newprob>$treeProb) {
									push @changeProb_IncreaseProb,exp(($newprob-$treeProb)/$MCtemperature);
									push @changeData_IncreaseProb,[$tmpclusterData,$tmpDistMataa,$tmpMaxTreeHeight,$newprob];
								} else {
									push @changeProb_DecreaseProb,exp(($newprob-$treeProb)/$MCtemperature);
									push @changeData_DecreaseProb,[$tmpclusterData,$tmpDistMataa,$tmpMaxTreeHeight,$newprob];
								}
							}
						};
						push @moves,$f1;
					}
				};
				
				@moves = shuffle @moves;
				$moves[0]->();				
				$f_update_tree_topology->();				
			}
		};
		# update the termination point
		if ($treeProb>$sav_score+1) {
			($sav_score,$savTree,@savPara) = ($treeProb,gen_newick($clusterData),($mu,$rho,$theta,$deltaTE));
			$step_to_stop = ($curStep+$N_step_relaxation>$minMCMCstep) ? $curStep+$N_step_relaxation : $minMCMCstep;
		}
		if ($treeProb>$maxx_score) {
			($maxx_score,$maxxTree,@maxxPara) = ($treeProb,gen_newick($clusterData),($mu,$rho,$theta,$deltaTE));			
		}
		
		if ($curStep%50000==0) {
			$MCtemperature = $MCtemperature/2;
		}
	}
	
	my $printstep = $step0 + $curStep;
print "step=$printstep\tlikelihood=$maxx_score\tmu=$maxxPara[0]\trho=$maxxPara[1]\ttheta=$maxxPara[2]\tdeltaTE=$maxxPara[3]\tdt=$dt\n";
print "$maxxTree;\n";
	print $fout_score "$printstep\t$maxx_score\t$maxxPara[0]\t$maxxPara[1]\t$maxxPara[2]\t$maxxPara[3]\t$dt\n";
	print $fout_tree "$maxxTree;\n";
	$fout_score->autoflush;
	$fout_tree->autoflush;
}



# calculate the log posterior probability
# inputs:
#		negative-cross-entropy for each strain-pair at different coalescent time-step: 
#			$strainpair2Prob_hash->{$strain1}{$strain2}[$timeStep] = $negativeCrossEntropy
#		distance matrix of strain pairs:
#			$distmap->{$strain1}{$strain2} = $timeStep
sub get_tree_probability2 {
	my ($strainpair2Prob_hash,$distmap) = @_;	
	my $totHeightProb = 0;	
	my $totTopologyProb = 0;
	my $totTopologyProb_NinternalNode = 0;
	
	my $strain = [keys %$distmap];
	my $totprob = 0;
	for my $q (0..$#$strain-1) {
		my $sn1 = $strain->[$q];
		for my $r ($q+1..$#$strain) {
			my $sn2 = $strain->[$r];
			my $tmpdist = $distmap->{$sn1}{$sn2};
			if ($tmpdist>$#{$strainpair2Prob_hash->{$sn1}{$sn2}}) {
				print STDERR "ererere $tmpdist $#{$strainpair2Prob_hash->{$sn1}{$sn2}} $sn1 $sn2\n";
			} else {
				$totprob += $strainpair2Prob_hash->{$sn1}{$sn2}[$tmpdist];
			}
		}
	}
	return $totprob;
}



# weighted random integer generator
sub choose_index_from_prob_array {
	my @prob = @_;
	my @cumsum = @prob;
	for (1..$#cumsum) {
		$cumsum[$_] = $cumsum[$_] + $cumsum[$_-1];
	}
	my $F = rand()*sum(@prob);
	my @ii = grep {$F<=$cumsum[$_]} 0..$#cumsum;
	return $ii[0];
}



# convert a distance matrix has into a ultrametric tree cluster object, the input is a distance matrix, the conversion is done by single linkage clustering
# format of input
#	$input_dist_matrix_hash->{$strain1}{$strain2} = $distance
# format of output (ultrametric tree cluster object)
# 	$inputClusterData->{$currentNode}{upstreamCluster} = $parentNode
#	$inputClusterData->{$currentNode}{downstreamCluster}{$cildNode1} = 1
#	$inputClusterData->{$currentNode}{downstreamCluster}{$cildNode2} = 1
#	$inputClusterData->{$currentNode}{tipnode}{$leafNode1} = 1
#	$inputClusterData->{$currentNode}{tipnode}{$leafNode2} = 1
#	$inputClusterData->{$currentNode}{height} = $heightOfNode
#	$inputClusterData->{$currentNode}{clusterTipPair}[$indexA][$indexB] = [$leafNode1,$leafNode2]
#		currentNode is the MRCA of leafNode1 and leafNode2
#		$indexA is the index for the childNode pairs, if currentNode has c childNodes, then $indexA ranges from 0 to c(c-1)/2
#		$indexB ranges from 0 to d-1, where d is the number of leaf-node-pairs that pass through child-node-pair $indexA
sub distmatrix2cluster {
	my ($input_dist_matrix_hash) = @_;
	#print "aa$input_dist_matrix_hash bb\n";
	my $distmat = dclone $input_dist_matrix_hash;
	my @tipnode = keys %$distmat;
	my $clusterData = {};
	for my $strain (@tipnode) {
		$clusterData->{$strain}{tipnode}{$strain} = 1;
		$clusterData->{$strain}{height} = 0;
	}
#print join("\t",'tipnode',@tipnode)."\n";
	
	my $get_min_element = sub {
		my ($dist__mat) = @_;
		my @cluster = keys %$dist__mat;
		my $minVal = [{},1e50];
		for my $q (0..$#cluster-1) {
			for my $r ($q+1..$#cluster) {
				if ($dist__mat->{$cluster[$q]}{$cluster[$r]} < $minVal->[1]) {
					$minVal->[1] = $dist__mat->{$cluster[$q]}{$cluster[$r]};
					$minVal->[0] = {};
					$minVal->[0]{$cluster[$q]} = 1;
					$minVal->[0]{$cluster[$r]} = 1;					
				}
			}
		}
		while (1) {
			my $hasAdded = 0;
			for my $node(@cluster) {
				next if exists $minVal->[0]{$node};
				for my $node2 (keys %{$minVal->[0]}) {
					if ($dist__mat->{$node}{$node2} == $minVal->[1]) {
						$minVal->[0]{$node} = 1;
						$hasAdded = 1;
						last;	# loop again if $hasAdded
					}
				}
				last if $hasAdded;	# loop again if $hasAdded
			}
			
			last unless $hasAdded;	# break the loop is not $hasAdded
		}
		return $minVal;
	};
	
	my $clusterID = 0;
	while (scalar(keys %$distmat)>1) {
		my $clusterHeight = $get_min_element->($distmat);
		my $curClusterName = 'c'.$clusterID++;
		foreach my $dsnode (keys %{$clusterHeight->[0]}) {
			$clusterData->{$dsnode}{upstreamCluster} = $curClusterName;
			$clusterData->{$curClusterName}{downstreamCluster}{$dsnode} = 1;
			foreach my $tipstrain (keys %{$clusterData->{$dsnode}{tipnode}}) {
				$clusterData->{$curClusterName}{tipnode}{$tipstrain} = 1;
			}
		}

		$clusterData->{$curClusterName}{height} = $clusterHeight->[1];
		
		foreach my $node2modify (grep {not exists $clusterHeight->[0]{$_}} keys %$distmat) {
			my $minValue = 1e50;
			foreach my $node2delete (keys %{$clusterHeight->[0]}) {
				if ($minValue > $distmat->{$node2modify}{$node2delete}) {
					$minValue = $distmat->{$node2modify}{$node2delete};
				}
				delete $distmat->{$node2modify}{$node2delete};
			}
			$distmat->{$node2modify}{$curClusterName} = $minValue;
			$distmat->{$curClusterName}{$node2modify} = $minValue;
		}
		foreach my $node2delete (keys %{$clusterHeight->[0]}) {
			delete $distmat->{$node2delete};
		}
	}
	
	my $outdistmat;
	my $maxHeight = 0;
	for my $curcluster (keys %$clusterData) {
		if (exists $clusterData->{$curcluster}{downstreamCluster}) {
			my @downstreamCluster = keys %{$clusterData->{$curcluster}{downstreamCluster}};
			for my $q (0..$#downstreamCluster-1) {
				for my $r ($q+1..$#downstreamCluster) {
					my @clusterTipPair;
					for my $sq (keys %{$clusterData->{$downstreamCluster[$q]}{tipnode}}) {
						for my $sr (keys %{$clusterData->{$downstreamCluster[$r]}{tipnode}}) {
							my $curHeight = $clusterData->{$curcluster}{height};
							$outdistmat->{$sq}{$sr} = $curHeight;
							$outdistmat->{$sr}{$sq} = $curHeight;
							push @clusterTipPair,[$sq,$sr];
							if ($curHeight>$maxHeight) {
								$maxHeight = $curHeight;
							}
						}
					}
					push @{$clusterData->{$curcluster}{clusterTipPair}},\@clusterTipPair;			
				}
			}
		}
	}
	
	return [$clusterData,$outdistmat,$maxHeight];
}





# convert a ultrametric tree cluster object into newick format
# format of ultrametric tree cluster object
# 	$inputClusterData->{$currentNode}{upstreamCluster} = $parentNode
#	$inputClusterData->{$currentNode}{downstreamCluster}{$cildNode1} = 1
#	$inputClusterData->{$currentNode}{downstreamCluster}{$cildNode2} = 1
#	$inputClusterData->{$currentNode}{tipnode}{$leafNode1} = 1
#	$inputClusterData->{$currentNode}{tipnode}{$leafNode2} = 1
#	$inputClusterData->{$currentNode}{height} = $heightOfNode
#	$inputClusterData->{$currentNode}{clusterTipPair}[$indexA][$indexB] = [$leafNode1,$leafNode2]
#		currentNode is the MRCA of leafNode1 and leafNode2
#		$indexA is the index for the childNode pairs, if currentNode has c childNodes, then $indexA ranges from 0 to c(c-1)/2
#		$indexB ranges from 0 to d-1, where d is the number of leaf-node-pairs that pass through child-node-pair $indexA
sub gen_newick {
	my ($inputClusterData) = @_;
	my $clusterData = dclone $inputClusterData;
	
	# fix a bug, as every internal node within $clusterData has its height shifted by 1
	my @internalNodes = sort {$a cmp $b} grep {exists $clusterData->{$_}{downstreamCluster}} keys %$clusterData;
	foreach (@internalNodes) {
		$clusterData->{$_}{height}++;
	}
	
	while (1) {
		my @lowestCluster = sort {$a cmp $b} grep {not exists $clusterData->{$_}{downstreamCluster}} keys %$clusterData;
		foreach my $curCluster (@lowestCluster) {
			if (not exists $clusterData->{$curCluster}{upstreamCluster}) {
				return '('.join(',',@{$clusterData->{$curCluster}{curNewick}}).'):0.0';
			}
			my $higherCluster = $clusterData->{$curCluster}{upstreamCluster};
			my $curBranchLength = $clusterData->{$higherCluster}{height} - $clusterData->{$curCluster}{height};
			my $tmptext = (exists $clusterData->{$curCluster}{curNewick}) ? '('.join(',',@{$clusterData->{$curCluster}{curNewick}}).')' : $curCluster;
			$tmptext .= ':'.$curBranchLength;
			push @{$clusterData->{$higherCluster}{curNewick}},$tmptext;
			delete $clusterData->{$curCluster};
			delete $clusterData->{$higherCluster}{downstreamCluster}{$curCluster};
			delete $clusterData->{$higherCluster}{downstreamCluster} if (scalar(keys %{$clusterData->{$higherCluster}{downstreamCluster}})==0);

		}
	}
}




