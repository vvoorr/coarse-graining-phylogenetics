package CGP_lib_SSP_theory;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);
use threads;
use strict;
use warnings;
use PDL;
$PDL::BIGPDL = 1;
use PDL::Ufunc;

use Exporter qw(import);
our @EXPORT_OK = qw/get_parameter_probability2 gen_evolution_matrix f_new_dist_map/;




my $dNstep = 2;








sub get_parameter_probability2 {
	my ($nsnpMatrix,$mu,$rho,$dN,$dTE,$endAvgSnp,$dt,$matrixSize,$use_rho_as_rhobare,$max_time_step) = @_;
	$matrixSize = 100 unless $matrixSize;
	my $pdlnsnp = pdl($nsnpMatrix)->transpose;
	
	$max_time_step+=2 if $max_time_step;
	my ($log_M_NSNP_timestep_pdf,$arrTime,$M_NSNP_timestep_cdf) = @{gen_evolution_matrix($mu,$rho,$dN,$dTE,$endAvgSnp,$dt,$matrixSize,$use_rho_as_rhobare,$max_time_step)};

	my $prob = $log_M_NSNP_timestep_pdf x $pdlnsnp;	
	$prob = $prob->transpose; # dimension of $prob is     N_distribution_Vectors x N_time_step
	my @arrayProb;	
	for my $q (0..$#$nsnpMatrix) {
		push @arrayProb,[list $prob->range([0,$q],[scalar(@$arrTime),1])];
	}
	return [\@arrayProb,$arrTime,$log_M_NSNP_timestep_pdf,$M_NSNP_timestep_cdf];
}



sub gen_evolution_matrix {
	my ($mu,$rho,$dN,$dTE,$endAvgSnp,$dt,$matrixSize,$use_rho_as_rhobare,$ending_time_step) = @_;
	$matrixSize = 100 unless $matrixSize;
	
	my $tmp_resolution = 5;	
	my $tmp_dt = ($dt>1/$mu/$tmp_resolution) ? 1/$mu/$tmp_resolution : $dt;
	if ($rho>0) {
		$tmp_dt = ($dt>1/$rho/$tmp_resolution) ? 1/$rho/$tmp_resolution : $tmp_dt;
	}
	if ($tmp_dt<$dt) {
		print STDERR "werwerwer\tdt too big, should be $tmp_dt, now is $dt\n";
	}
	
	# prepare the recombination matrix
	my $mm = xvals $matrixSize+1,$matrixSize+1;
	my $nn = yvals $matrixSize+1,$matrixSize+1;
	my $A = zeros($matrixSize+1,$matrixSize+1);
	my $f0 = exp(-1*yvals(1,$matrixSize+1)/$dN);
	$f0 = $f0 / sum($f0);
	my $F0 = cumusumover($f0->xchg(0,1))->transpose;
	
	$F0 = $F0/$F0->range([0,$matrixSize]);
	for my $q (0..$matrixSize-1) {
		$A->range([$q+1,$q],[$matrixSize-$q,1]) .= exp(-1*(xvals($matrixSize-$q,1)+$q+1)/$dTE)*(1-$F0->range([0,$q])->transpose)*$f0->range([0,$q])->transpose;
		if (not $use_rho_as_rhobare) {
			$A->range([$q,$q]) .= sum((1-$F0->range([0,0],[1,$q+1])) * $f0->range([0,0],[1,$q+1]) * exp(-1*yvals(1,$q+1)/$dTE));
		}
		$A->range([$q,$q+1],[1,$matrixSize-$q]) .= exp(-1*(yvals(1,$matrixSize-$q)+$q+1)/$dTE)*(1-$F0->range([0,$q]))*$f0->range([0,$q+1],[1,$matrixSize-$q]);
	}
	if (not $use_rho_as_rhobare) {
		do {
			my $q = $matrixSize;
			$A->range([$q,$q]) .= sum((1-$F0->range([0,0],[1,$q+1])) * $f0->range([0,0],[1,$q+1]) * exp(-1*yvals(1,$q+1)/$dTE));
		};
	}	
	if (not $use_rho_as_rhobare) {
		for my $q (0..$matrixSize) {
			for my $r (0..$matrixSize) {
				if ($A->range([$q,$r])<0) {
					$A->range([$q,$r]) .= 0;
				}			
			}
		}
		for my $q (0..$matrixSize) {
			$A->range([$q,0],[1,$matrixSize+1]) .= $A->range([$q,0],[1,$matrixSize+1]) / sum($A->range([$q,0],[1,$matrixSize+1]));
		}
	} else {
		for my $q (0..$matrixSize) {
			$A->range([$q,$q]) .= 1 - sum($A->range([$q,0],[1,$matrixSize+1]));
		}
		for my $q (0..$matrixSize) {
			for my $r (0..$matrixSize) {
				if ($A->range([$q,$r])<0) {
					$A->range([$q,$r]) .= 0;
				}			
			}
		}
	}	
	
	# prepare the mutation matrix
	my $B = zeros $matrixSize+1,$matrixSize+1;
	for my $q (0..$matrixSize-1) {
		$B->range([$q,$q+1]) .= 1;
	}
	$B->range([$matrixSize,$matrixSize]) .= 1;
	
	# identity matrix
	my $IDmatrix = identity $matrixSize+1;
	
	# evolutionary matrix
	my $RRho = $dt*(1-(1-$rho)**2);
	my $MMu = $dt*(1-(1-$mu)**2);
	die "mu:$mu and rho:$rho too big!!\n" if (1-$MMu-$RRho) < 0.5;
	my $evoM = (1-$MMu-$RRho) * $IDmatrix + $RRho * $A + $MMu * $B;
		
	my @time;
	my @record_vector;
	my $t = $dt;
	my $curM = $evoM;
	# print log $curM;
	
	
	my $v0 = zeros 1,1+$matrixSize;
	$v0->range([0,0]) .= 1;
	my $xval = xvals 1+$matrixSize,1;
	while (1) {
		push @time,$t;
		my $tmpM = $curM->range([0,0],[1,1+$matrixSize]);
		push @record_vector,$tmpM;
		my $avgDelta = $xval x $curM x $v0;
		
		#if (defined($ending_time_step) && ($t>$ending_time_step*$dt) && ($avgDelta>$endAvgSnp)) {
		if (defined($ending_time_step) && ($t>$ending_time_step*$dt)) {
			last;
		} elsif (not(defined($ending_time_step)) && ($avgDelta>$endAvgSnp)) {
			last;
		} else {
			$curM = $curM x $evoM;
			$t += $dt;
		}
		last if $t > 10000;
	}
	my $outputM = zeros scalar(@record_vector),1+$matrixSize;
	for my $q (0..$#record_vector) {
		$outputM->range([$q,0],[1,1+$matrixSize]) .= $record_vector[$q]->range([0,0],[1,1+$matrixSize]);
	}	
	$outputM = $outputM->transpose();

	my $output_M_cdf = dcumusumover $outputM;	
	
	return [log($outputM),\@time,$output_M_cdf];
}



sub f_new_dist_map {
	my ($old_dist_map,$new_dist,$arr_tip_pair) = @_;
	my $tmpDistMap = dclone $old_dist_map;
	for my $tmp (@$arr_tip_pair) {
		my ($s1,$s2) = @$tmp;
		$tmpDistMap->{$s1}{$s2} = $new_dist;
		$tmpDistMap->{$s2}{$s1} = $new_dist;
	}
	return $tmpDistMap;
}











1;
