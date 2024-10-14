#!/usr/bin/perl
## here does the P(miss) - Fisher's test on NA counts (no calls or fails DPGQ)
## check last 8 columns in sampleGenotype sorted file are as follows:
## Control_NA      Control_0       Control_1       Control_2       Case_NA Case_0  Case_1  Case_2
##
## NOTE: check the columns
##
## - added in modified Haldane-Anscombe correction to get the odds ratio which involves adding 0.5 to each cell count
## whenever there is a 0-value cell in 2x2 table. Ensures no division by 0.
#######################################################################################################################################

use Text::NSP::Measures::2D::Fisher::twotailed; 
use strict;
use warnings;

open (INPUT1, "$ARGV[0]") or die "Cannot open $ARGV[0]\n";
open (OUTPUT, ">$ARGV[1]") or die "Cannot open $ARGV[1]\n";

## reading in the input, extracting info of LAST 8 columns to store in hash in following manner
## Control_NA	Control_0	Control_1	Control_2	Case_NA	Case_0	Case_1	Case_2

## filling up the hash like
#     		cases   |  controls
# NA count 	  n11		n12   | n1p 
# no NA count	  n21		n22   | n2p
# 		---------------------------
#		  np1		np2     npp 

my %fisher_hash_variant;

my $header=<INPUT1>;
chomp($header);
print OUTPUT "$header\tP_miss\tP_miss_OddRatio\tP_miss_Log_OR\tP_miss_abs_logOR\ttotal\n";


while (my $line1=<INPUT1>) {
     chomp($line1);
     my @v=split(/\t/,$line1);

     my $control_NA = $v[-8];        ##genotype code: "NA"
     my $control_hommajor = $v[-7];  ##genotype code: "0"
     my $control_het = $v[-6];	     ##genotype code: "1"
     my $control_homminor = $v[-5];  ##genotype code: "2"
     my $case_NA = $v[-4];
     my $case_hommajor = $v[-3];
     my $case_het = $v[-2];
     my $case_homminor = $v[-1];

     my $total_case=$case_NA + $case_hommajor + $case_het + $case_homminor;
     my $total_control=$control_NA + $control_hommajor + $control_het + $control_homminor;

     my $case_NA_count = $case_NA;
     my $case_nonNA_count = $total_case - $case_NA_count;
     my $control_NA_count = $control_NA;
     my $control_nonNA_count = $total_control - $control_NA_count;

     my $total_NA_count = $case_NA_count + $control_NA_count;
     my $total_nonNA_count = $case_nonNA_count + $control_nonNA_count;
     my $total_n= $total_NA_count + $total_nonNA_count;


     ##calculating the fisher's test p-value, at allele level for non-zero 2x2 cell values
     my $p_value="NA";
     
     if($case_NA_count ==0 and $control_NA_count ==0 and $case_nonNA_count ==0 and $control_nonNA_count ==0) {
        $p_value=1;
     }         	
     else{
	$fisher_hash_variant{"n11"} = $case_NA_count;
     	$fisher_hash_variant{"n12"} = $control_NA_count;
     	$fisher_hash_variant{"n1p"} = $total_NA_count;

     	$fisher_hash_variant{"n21"} = $case_nonNA_count;
     	$fisher_hash_variant{"n22"} = $control_nonNA_count;
     	$fisher_hash_variant{"n2p"} = $total_nonNA_count;

     	$fisher_hash_variant{"np1"} = $case_NA_count + $case_nonNA_count;
     	$fisher_hash_variant{"np2"} = $control_NA_count + $control_nonNA_count;
     	$fisher_hash_variant{"npp"} = $total_n;

     	$p_value = calculateStatistic(%fisher_hash_variant);
     }

     ##uses a modified Haldane-Anscombe correction to get the odds ratio
     ##if there is any 0 value in 2x2 table of counts, add 0.5 to each cell count
     ##else odds ratio calculation remains the same
     my $oddratio=0;
     my $logodd=0;
     my $abs_log=0;

     if( $case_NA_count==0 || $control_nonNA_count==0 || $case_nonNA_count==0 || $control_NA_count==0 ) {
	my $case_NA_count_corr=$case_NA_count+0.5;
	my $control_nonNA_count_corr=$control_nonNA_count+0.5;
	my $case_nonNA_count_corr=$case_nonNA_count+0.5;
	my $control_NA_count_corr=$control_NA_count+0.5;

	$oddratio=($case_NA_count_corr * $control_nonNA_count_corr) / ($case_nonNA_count_corr * $control_NA_count_corr);
	$logodd=log($oddratio);
	$abs_log=abs($logodd);
     }
     else{
	$oddratio=($case_NA_count * $control_nonNA_count) / ($case_nonNA_count * $control_NA_count);
	$logodd=log($oddratio);
	$abs_log=abs($logodd);
     }

     print OUTPUT "$line1\t$p_value\t$oddratio\t$logodd\t$abs_log\t$total_n\n";
}
close(INPUT1);
close(OUTPUT);

