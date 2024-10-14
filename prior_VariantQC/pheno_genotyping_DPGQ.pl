#!/usr/bin/perl -w
## script generates samplegenotype file
########################################################################

no warnings;
use Pod::Usage;
use Getopt::Long;

my $snpFile = $ARGV[0]; #variant.vcf file
my $sampleFile = $ARGV[1]; #sample file
my $location =  $ARGV[2]; #output location
my $outfilePhenotype = $location."phenotyping_variants";
my $outfileGenotype = $location."samplegenotype_unsorted";

our ($help);
GetOptions ( 'help|h'=>\$help) or pod2usage ();
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);

my %hashTableSAMPLELIST = ();
my %hashTableGenotype = ();
my @phenotype = ();
open (OUTPHENOTYPE, ">$outfilePhenotype") || die "can't open $outfilePhenotype";
open (OUTGENOTYPE, ">$outfileGenotype") || die "can't open $outfileGenotype";

open SAMPLEFILE, $sampleFile or die $!;

	while (<SAMPLEFILE>) 
	{
			&getHashTableSampleFile($_);
	}

close(SAMPLEFILE);	

print OUTGENOTYPE "#Chr\tPosition\tRef\tAlt";

open SNPFILE, $snpFile or die $!;

	while (<SNPFILE>) 
	{
		if ($.==1 and /^# autoFile/ )
		{
			next;
		}
		if (/^##/){
			next;
		}
		if (/^\"/){
			next;
		}
		
		if(/^#/){
			my @lines = split (/\t/, $_);
			my @samples = @lines[9..$#lines];
			my $index=0;
			foreach (@samples) 
				{ 
					chomp($_);
					my $sname=$_;
					
					#my $gSname = substr $sname, 0, 6;
					my $gSname =$sname;
					#print "gSname:$gSname\n";
					foreach my $key (keys %hashTableSAMPLELIST)
					{
					
						if ($gSname =~ $key) {
					#			print "key=$key,gsname=$gSname\n";
							my @hpuVal = @{$hashTableSAMPLELIST{$key}};
							print OUTPHENOTYPE $key,"\t".$hpuVal[0]."\n";
							print OUTGENOTYPE "\t$key:$hpuVal[0]";
							$phenotype[$index]=$hpuVal[0];	
							$index++;
							last;
						}
					}

				}
				print OUTGENOTYPE "\tAllele_Freq_Control\tAllele_Freq_Case\tAllele_Freq_All_Samples\tMissing_Call_Rate_Control\tMissing_Call_Rate_Case\tMissing_Call_Rate_All\tFilter\tControl_NA\tControl_0\tControl_1\tControl_2\tCase_NA\tCase_0\tCase_1\tCase_2\n";
				next;
		}
		
		&getGenotype($_);

	}

close(SNPFILE);	


while ((my $k, my $v) = each(%hashTableGenotype)){     
	my @keyVal = split(/l/, $k);
	print OUTGENOTYPE "$keyVal[0]\t$keyVal[1]\t$keyVal[2]\t$keyVal[3]";

	foreach (@{$v}) 
	{ 
		#Using an array reference as an array see perlreftut        
		print OUTGENOTYPE "\t$_";     	
	}
	print OUTGENOTYPE "\n";
					
}


sub getHashTableSampleFile()
{
	my @t = split (/[ \t]+/, $_);
	my $sampleId = $t[0];
	my $phenotype = $t[1];
	chomp($phenotype);
	
	if( $phenotype eq "control")
	{
		push @{$hashTableSAMPLELIST{$sampleId}},  0;
	}
	else
	{
		push @{$hashTableSAMPLELIST{$sampleId}},  1;
	}
}


sub getGenotype()
{
	
	my @cols = split (/\t/, $_);
	my $chr = $cols[0];
	my $pos = $cols[1];
	my $ref = $cols[3];
	my $alt = $cols[4];
	my $filter= $cols[6]; #VQSR Filter status
	
	my $format = $cols[8];
	my @samples = @cols[9..$#cols];
		
	my $totalSample=0;
	my $alleleControl=0;
	my $alleleCase=0;
	my $totalCase=0;
	my $totalControl=0;
	my $alleleFreqCountCase=0;
	my $alleleFreqCountControl=0;
	my $alleleFreqForAllSample=0;

	my $missingSampleControl=0;
	my $missingSampleCase=0;
	my $phenoControl=0;
	my $phenoCase=0;
	

	my $percentageOfMissingControl=0;
	my $percentageOfMissingCase=0;
	my $percentageOfMissingAll=0;
	my $alleleFreqControl=0;
	my $alleleFreqCase=0;


	my $countNAControl=0;   ##count # of "NA" nocalls/failDPGQ genotype in controls
	my $count0Control=0;    ##count # of "0" hom for ref allele in controls
	my $count1Control=0;    ##count # of "1" het genotype in controls
	my $count2Control=0;    ##count # of "2" hom for alt allele in controls
	my $countNACase=0;    ##count # of "NA" nocalls/failDPGQ genotype in cases
	my $count0Case=0;    ##count # of "0" hom for ref allele in cases
	my $count1Case=0;   ##count # of "1" het genotype in cases
	my $count2Case=0;   ##count # of "2" hom for alt allele in cases
	

	my $colDepth=2;
	my $colQS=3;
	chomp($format);
	my @formatlist = split (/:/, $format);
	my $index=0;
	my $i=0;

	foreach(@formatlist)
	{
		if($_ eq 'DP')
		{
			$colDepth=$index;
		}
		if($_ eq 'GQ')
		{
			$colQS=$index;
		}
		$index++;
	}

	foreach (@samples) 
	{
		chomp;
		my $sample =$_;
		$totalSample++;
		my $pheno = $phenotype[$i];

			my @individualdata = split(/:/, $sample);
			my $homozygousity="NA";
			
			if ($pheno ==0) 
			{
				$phenoControl++;
			}
			if ($pheno ==1) 
			{
				$phenoCase++;
			}


			##also does DP, GQ filter
			if (int($individualdata[$colDepth]) >=8 and int($individualdata[$colQS]) >=20)
			{
				
				if ((substr $individualdata[0],0,3)  eq "0/0" or (substr $individualdata[0],0,3)  eq "0|0") 
				{
					$homozygousity= "0"; #match
					#next;
					if($pheno == 1)
					{
						$totalCase++;
					}
					else
					{
						$totalControl++;
					}
				} 
				
				if ((substr $individualdata[0],0,3)  eq "1/1" or (substr $individualdata[0],0,3) eq "1|1")
				{
					$homozygousity= "2"; #hom
					
					if($pheno == 1)
					{
						$alleleFreqCountCase =$alleleFreqCountCase+2;
						$totalCase++;
					}
					else
					{
						$alleleFreqCountControl =$alleleFreqCountControl+2;
						$totalControl++;
					}
				}
			
				if ((substr $individualdata[0],0,3) eq "0/1" or (substr $individualdata[0],0,3) eq "0|1")
				{
					$homozygousity = "1"; #het		
					
					if($pheno == 1)
					{
						$alleleFreqCountCase =$alleleFreqCountCase+1;
						$totalCase++;
					}
					else
					{
						$alleleFreqCountControl =$alleleFreqCountControl+1;
						$totalControl++;
					}
				}
				
			}

			##counting no of samples with "NA" genotype
			if ($homozygousity eq "NA" && $pheno ==0) 
			{
				$missingSampleControl++;
				$countNAControl++;
			}
			if ($homozygousity eq "NA" && $pheno ==1) 
			{
				$missingSampleCase++;
				$countNACase++;
			}


			##counting no of samples with "0" hom ref allele genotype
			if ($homozygousity eq "0" && $pheno ==0)
			{
				$count0Control++;
			}
			if ($homozygousity eq "0" && $pheno ==1)
			{	
				$count0Case++;
			}


			##counting no of samples with "1" het genotype
			if ($homozygousity eq "1" && $pheno ==0)
                	{
                        	$count1Control++;
                	}
                	if ($homozygousity eq "1" && $pheno ==1)
                	{
                        	$count1Case++;
                	}


			##counting no of samples with "2" hom alt allele genotype
			if ($homozygousity eq "2" && $pheno ==0)
			{
				$count2Control++;
			}
			if ($homozygousity eq "2" && $pheno ==1)
                	{
                        	$count2Case++;
                	}

			
			push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $homozygousity;
			$i++;
	}
	if($missingSampleControl != 0 and $phenoControl !=0)
	{
		$percentageOfMissingControl = $missingSampleControl/$phenoControl;
	}
	if($missingSampleCase != 0 and $phenoCase !=0)
	{
		$percentageOfMissingCase = $missingSampleCase/$phenoCase;
	}	
	if($missingSampleCase != 0 or $missingSampleControl !=0)
	{
		$percentageOfMissingAll = ($missingSampleCase + $missingSampleControl) / ($phenoCase + $phenoControl );
	}
	if($alleleFreqCountCase != 0 and $totalCase !=0)
	{
	
		$alleleFreqCase=$alleleFreqCountCase/($totalCase * 2);
	}
	if($alleleFreqCountControl != 0 and $totalControl !=0)
	{
	
		$alleleFreqControl=$alleleFreqCountControl/($totalControl * 2);
	}

	if (($alleleFreqCountControl + $alleleFreqCountCase) != 0 and ($totalControl+$totalCase) !=0)
	{
		$alleleFreqForAllSample = ($alleleFreqCountControl + $alleleFreqCountCase) /  (($totalControl+$totalCase)*2);
	}


	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $alleleFreqControl;
	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $alleleFreqCase;	
	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $alleleFreqForAllSample;	
	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $percentageOfMissingControl;
	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $percentageOfMissingCase;
	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $percentageOfMissingAll;
	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $filter;

	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $countNAControl;
	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $count0Control;	
	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $count1Control;	
	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $count2Control;
	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $countNACase;
	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $count0Case;
	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $count1Case;
	push @{$hashTableGenotype{$chr."l".$pos."l".$ref."l".$alt}},  $count2Case;

}
exit (0);

=head1 SYNOPSIS

 pheno_genotyping_DPGQ.pl [arguments] <Variants.vcf> <sample_file.txt> <OutputFolder/>

 Optional arguments:
 	-h, --help			print help message
 	
 Function: To produce Sample Genotype file (each variant and its genotype per-sample site)
 
 Example: perl pheno_genotyping_DPGQ.pl Variants.vcf sample.txt OutputFolder/
 Argument2 -> Sample file (tab delimited) Format: SampleName tab case or control"
 

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit
