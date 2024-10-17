#!/usr/bin/perl
## script here appends hwe pval (controls) to samplegenotype_sorted_pmiss
## hwe output i.e. all-variants-fisher-hwe.hwe
#######################################################################################

use strict;
use warnings;

open (INPUT1, "$ARGV[0]") or die "Cannot open $ARGV[0]\n";  ##plink hwe results
open (INPUT2, "$ARGV[1]") or die "Cannot open $ARGV[1]\n";  ##samplegenotype_sorted_pmiss 
open (OUTPUT1,">$ARGV[2]") or die "Cannot open >$ARGV[2]\n";   ##combined variant metric file ie samplegenotype_sorted_pmiss_hwe

my %hash_hwe;


## reads in plink hwe results and only store hwe control pval in hash
## format: 	CHR	SNP	TEST	A1	A2	GENO	O(HET)	E(HET)	P
##	 1	indelIDb1	ALL	AC	A	0/1/3474	0.0002878	0.0002877	1
## 	1	indelIDb1	AFF	AC	A	0/0/0	nan	nan	NA
## 	1	indelIDb1	UNAFF	AC	A	0/1/3474	0.0002878	0.0002877	1    -> take this pval (for controls)
## there is a hidden additional space in front of each line

while(my $line1=<INPUT1>) {
	chomp($line1);

	##skips header line
        if($line1=~/^CHR/) {
                next;
        }
	else{
                my @t=split(/\s{1,}/,$line1);
                my $snpid=$t[2];  ##there is a hidden additional space in front of each line
                my $test=$t[3];
		my $hwe_p=$t[9];

		#print $test,"\n";
		#print $hwe_p,"\n";

		##only stores hwe p-val of controls
		if($test=~/^UNAFF/){
                	if(!exists $hash_hwe{$snpid}) {
                        	$hash_hwe{$snpid}=$hwe_p;
                	}
			else{;}
		}
                else{;}
        }
}
close(INPUT1);


## now reads in samplegenotype_sorted_pmiss file and does matches by snpid (format chr:pos:ref:alt)
## retrieves corresponding hwe pval control and append col behind
my $header=<INPUT2>;
chomp($header);
print OUTPUT1 "$header\tHWE_ctrl_pval\n";

my $nosnpid_ct=0;
my $nohwepval_ct=0;

while(my $line2=<INPUT2>) {
	chomp($line2);
	my @v=split(/\t/,$line2);
	my $chr=$v[0];
	$chr=~s/^chr//;
	my $pos=$v[1];
	my $ref=$v[2];
	my $alt=$v[3];

	my $key="$chr:$pos:$ref:$alt";

		
	if(exists $hash_hwe{$key}) {
		print OUTPUT1 "$line2\t$hash_hwe{$key}\n";
	}
	else{
		$nohwepval_ct++;
		print OUTPUT1 "$line2\tNA\n";
	}
}
close(INPUT2);
close(OUTPUT1);


##checks
print "No. of variants without hwe_pval:$nohwepval_ct\n";


