#!/usr/bin/perl -w
## calculates %het_snps, singleton counts for each sample - sampleQC stage
##
##########################################################################################

use strict;
use warnings;
no warnings 'uninitialized';
no warnings 'numeric';


use Pod::Usage;
use Getopt::Long;

my $tsvFile = $ARGV[0]; #variantannot_samplegeno_final_all_results_filtered.tsv
my $gvFile = $ARGV[1]; #samplegeno_filtered.tsv
my $location =  $ARGV[2]; #output location
my $outfile = $location."all_samples_summary";

our ($help);
GetOptions ( 'help|h'=>\$help) or pod2usage ();
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);

die "Usage: $0 PATTERN [variantannot_samplegeno_final_all_results_filtered.tsv samplegeno_filtered.tsv outputFolder/]\n" if @ARGV < 3;
die "$tsvFile file not exists\n" if (!(-f $tsvFile));
die "$gvFile file not exists\n" if (!(-f $gvFile));

unless(-d $location){
    mkdir $location or die;
}

my %hashTableVariantsLIST = ();
my %hashTableGNLIST=();
my $totalSample=0;
my @sampleId = ();
my $afCol;
my %hashTableSingletons_no_dbsnpid=();
my $tsv_header;
my %hashTableTSV=();

open (OUTSAMPLE_SUMMARY, ">$outfile") || die "can't open $outfile\n";

open GVFILE, $gvFile or die $!;
	while (<GVFILE>) 
	{
		chomp($_);

		if ($.==1)
		{
				my @t = split (/\t/, $_);
				my $size = @t;
				$afCol = &getNumberOfAFColumns(\@t);
				#$totalSample=$size-(int($afCol)+4); #from $totalSample=$size-11
				$totalSample=$size-(int($afCol)+5);  ##added transcriptid in GV cols; 02/03/2022

				for(my $i = 5; $i < $size-$afCol; $i++) ##changed starting idx from 4 to 5; added transcriptid in GV cols; 02/03/2022
				{
					print OUTSAMPLE_SUMMARY "\t".$t[$i];
					push (@sampleId, $t[$i]);
				}
				print OUTSAMPLE_SUMMARY "\n";
				next;
		}
		&getHashTableGVFile($_);
	}

close(GVFILE);	

sub getNumberOfAFColumns
{
	my ($aref) = @_;
	my ($total) = 0;
	my $keyword ="Allele_Freq_Control|Allele_Freq_Case|Allele_Freq_All_Samples|Missing_Call_Rate_Control|Missing_Call_Rate_Case|Missing_Call_Rate_All|Filter|Control_NA|Control_0|Control_1|Control_2|Case_NA|Case_0|Case_1|Case_2|P_miss|P_miss_OddRatio|P_miss_Log_OR|P_miss_abs_logOR|total|HWE_ctrl_pval";
	foreach (@$aref) { 
		chomp($_);
		if($_ =~ m/($keyword)/i){
			$total=$total+1;
		}
	}
	return $total;
}

sub getHashTableGVFile()
{
	chomp($_);
	my @t = split (/\t/, $_);
	my $chr = $t[0];
	$chr =~ s/chr//;
	my $pos= $t[1];
	my $ref= $t[2];
	my $alt= $t[3];

	my $arraySize = @t;
	my $str="";
	my $count=0;
	my $sample=0;
	#for(my $i = 4; $i < $arraySize-$afCol; $i++) #arraySize not counting 0
	for(my $i = 5; $i < $arraySize-$afCol; $i++)   ##changed starting indx from 4 to 5; added transcriptid in GV cols; 02/03/2022
	{
		$str .=$t[$i]."@";
		#added for singleton
		if($t[$i]==1 or $t[$i]==2)
		{
			$count++;
			#$sample=$i-4;  
			$sample=$i-5;   ##added transcriptid in GV cols; 02/03/2022
		}
		########
	}
	#add for singleton	
	if($count==1)
	{
		$str .=$sample."@";
	}
	else
	{
		$str .="-1@";
	}
########
	push @{$hashTableGNLIST{$chr."@".$pos."@".$ref."@".$alt}},  $str;
	
}

open TSVFILE, $tsvFile or die $!;

	while (<TSVFILE>) 
	{
		chomp($_);
		if ($.==1) {
		    $tsv_header=$_;
		    next;
		}
		&getTSVVarFile($_);
	}

close(TSVFILE);	

&computeSummarize;


sub computeSummarize 
{
	my @arr_syno_1 = ();
	my @arr_syno_2 = ();
	my @arr_non_syno_1 = ();
	my @arr_non_syno_2 = ();

	my @arr_singleton= ();
	my @arr_singleton_no_dbsnpid=();

 	foreach my $mykey (keys %hashTableGNLIST)
	{
		chomp($mykey);
			
			my @hpuVal = @{$hashTableGNLIST{$mykey}};
			if ( exists $hashTableVariantsLIST{$mykey} )
			{
				my @t = split (/@/, $hpuVal[0]);
				my $index=0;
				#print "hpuval0:$hpuVal[0]\n";
				my $intSize=@t;

				foreach(@t)
				{
					my @val =  @{$hashTableVariantsLIST{$mykey}};
					my @criteria = split (/@/, $val[0]);
					
					my $prev_syno_1=$arr_syno_1[$index];
					if (!defined($prev_syno_1)){$prev_syno_1=0;$arr_syno_1[$index]=0;}
					
					my $prev_syno_2=$arr_syno_2[$index];
					if (!defined($prev_syno_2)){$prev_syno_2=0;$arr_syno_2[$index]=0;}
						
					my $prev_non_syno_1=$arr_non_syno_1[$index];
					if (!defined($prev_non_syno_1)){$prev_non_syno_1=0;$arr_non_syno_1[$index]=0;}
										
					my $prev_non_syno_2=$arr_non_syno_2[$index];
					if (!defined($prev_non_syno_2)){$prev_non_syno_2=0;	$arr_non_syno_2[$index]=0;}
				

					my $prev_singleton=$arr_singleton[$index];
					if (!defined($prev_singleton)){$prev_singleton=0;$arr_singleton[$index]=0;}

					my $prev_singleton_no_dbsnpid=$arr_singleton_no_dbsnpid[$index];
					if (!defined($prev_singleton_no_dbsnpid)){$prev_singleton_no_dbsnpid=0;$arr_singleton_no_dbsnpid[$index]=0;}		

					##consider all consequences for synonymous, missense, start_loss, stop_gain, stop_loss
					if ( ( ($criteria[0] =~ m/^synonymous/) or ($criteria[0] =~ m/^missense/) or ($criteria[0] =~ m/^start_loss/) or ($criteria[0] =~ m/^stop_gained/) or ($criteria[0] =~ m/^stop_loss/) ) and ($intSize != $index)) #not to include last column(singleton info)
					{

						if ($criteria[0] =~ m/^synonymous/) 
						{
							
							if($t[$index] eq '1') {$arr_syno_1[$index]=int($prev_syno_1)+1;}
							if($t[$index] eq '2') {$arr_syno_2[$index]=int($prev_syno_2)+1;}

							#print "after assignment index:$index,arr_syno_1:$arr_syno_1[$index]\n";

						}

						if ( ($criteria[0]  =~ m/^missense/) or ($criteria[0] =~ m/^start_loss/) or ($criteria[0] =~ m/^stop_gained/) or ($criteria[0] =~ m/^stop_loss/) )
						{
							if($t[$index] eq '1') {$arr_non_syno_1[$index]=int($prev_non_syno_1)+1;}
							if($t[$index] eq '2') {$arr_non_syno_2[$index]=int($prev_non_syno_2)+1;}
						}
						
					}



					##also counts singletons that are potentially novel (without snpid assigned);
					##singletons count only for synonymous, missense, start_lost, stop_gain, stop_loss
					if ( ($criteria[0] =~ m/^synonymous/) or ($criteria[0] =~ m/^missense/) or ($criteria[0] =~ m/^start_loss/) or ($criteria[0] =~ m/^stop_gained/) or ($criteria[0] =~ m/^stop_loss/) ) {

						if($t[$#t] == $index)
						{	
							$arr_singleton[$index]=int($prev_singleton)+1;

							if($criteria[1] =~ m/^NA/) {
						     		$arr_singleton_no_dbsnpid[$index]=int($prev_singleton_no_dbsnpid)+1;
						     		push @{$hashTableSingletons_no_dbsnpid{$index}},$mykey;
 							}
						}
					}
					$index++;
				}
			}
	}


	##PRINT OUT THE RESULTS
	print OUTSAMPLE_SUMMARY "NonSynonymous_1";
	for(my $f=0; $f<$totalSample; $f++) { print OUTSAMPLE_SUMMARY "\t".$arr_non_syno_1[$f];}
	print OUTSAMPLE_SUMMARY "\n";

	print OUTSAMPLE_SUMMARY "NonSynonymous_2";
	for(my $f=0; $f<$totalSample; $f++) { print OUTSAMPLE_SUMMARY "\t".$arr_non_syno_2[$f];}
	print OUTSAMPLE_SUMMARY "\n";

	print OUTSAMPLE_SUMMARY "Synonymous_1";
	for(my $f=0; $f<$totalSample; $f++) { print OUTSAMPLE_SUMMARY "\t".$arr_syno_1[$f];}
	print OUTSAMPLE_SUMMARY "\n";

	print OUTSAMPLE_SUMMARY "Synonymous_2";
	for(my $f=0; $f<$totalSample; $f++) { print OUTSAMPLE_SUMMARY "\t".$arr_syno_2[$f];}
	print OUTSAMPLE_SUMMARY "\n";

	
	##for %het snps calculation
	print OUTSAMPLE_SUMMARY "%Het_snps";	
	for(my $f=0; $f<$totalSample; $f++) {
		my $numerator=int($arr_non_syno_1[$f])+int($arr_syno_1[$f]);
		my $denominator=int($arr_non_syno_1[$f])+int($arr_non_syno_2[$f])+int($arr_syno_1[$f])+int($arr_syno_2[$f]);
		if ($numerator != 0 and $denominator != 0){
			print OUTSAMPLE_SUMMARY "\t".($numerator)/($denominator);
		}
		else {
			print OUTSAMPLE_SUMMARY "\t0";
		}

	}
	print OUTSAMPLE_SUMMARY "\n";


	print OUTSAMPLE_SUMMARY "SINGLETON";
	for(my $f=0; $f<$totalSample; $f++) { print OUTSAMPLE_SUMMARY "\t".$arr_singleton[$f];}
	print OUTSAMPLE_SUMMARY "\n";

	print OUTSAMPLE_SUMMARY "SINGLETON(without_snpid)";
        for(my $f=0; $f<$totalSample; $f++) { print OUTSAMPLE_SUMMARY "\t".$arr_singleton_no_dbsnpid[$f];}
        print OUTSAMPLE_SUMMARY "\n";


}

sub getTSVVarFile()
{
	chomp($_);
	my @t = split (/\t/, $_);
	my $chr = $t[0];
	$chr =~ s/chr//;
	my $pos= $t[1];
	my $ref= $t[2];
	my $alt= $t[3];	

	##consequence: col 26, dbsnp: col 34
	my $snpType = $t[25];
	my $dbsnp= $t[33];

	if($dbsnp =~/^\-/) {
		$dbsnp="NA";
	}
	else{;}
	

	my $res ="";
	$res = $res.$snpType."@".$dbsnp;
	

	push @{$hashTableVariantsLIST{$chr."@".$pos."@".$ref."@".$alt}},  $res;
	push @{$hashTableTSV{$chr."@".$pos."@".$ref."@".$alt}}, $_;
}



exit (0);

=head1 SYNOPSIS

 perl sample_summarize.pl [arguments] <variantannot_samplegeno_final_all_results_filtered.tsv> <samplegeno_filtered.tsv> <OutputFolder/>

 Optional arguments:
 	-h, --help			print help message
 	
 Function: Calculates %het_snps, singleton counts per sample in preQC dataset. Part of sampleQC step.
 
 Example: perl sample_summarize.pl <preQC_CDS dir/variantannot_samplegeno_final_all_results_filtered.tsv> <preQC_CDS dir/samplegeno_filtered.tsv> <OutputFolder/>
 
 

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit
