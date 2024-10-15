#!/usr/bin/perl
## internal script that does actual filtering of variant quality & variant annotation
#############################################################################################################################################


my $annotationResults = $ARGV[0]; ##ie variantannot_samplegeno.tsv
my $seqVariants = $ARGV[1]; ##samplegeno.tsv
my $parameter	=  $ARGV[2]; ##parameters file (refer to README criteria format)
my $outdir = $ARGV[3];
my $outAnnotations = $outdir . "variantannot_samplegeno_final_all_results_filtered.tsv";
my $outGenotyping = $outdir . "samplegeno_filtered.tsv";

open (OUTANNOTATION, ">$outAnnotations") || die "can't open $outAnnotations";
open (OUTGENOTYPING, ">$outGenotyping") || die "can't open $outGenotyping";

my %hashTableVariant = ();
my %hashTableFiltering = ();
my %hashTableAnnoHeader = ();
my %hashTableAnnoFilter=();
my @phenotype=();

my $afCol;


push @{$hashTableAnnoHeader{"GNOMAD_EXOMES_AF"}},7;
push @{$hashTableAnnoHeader{"GNOMAD_EXOMES_AFR"}},8;
push @{$hashTableAnnoHeader{"GNOMAD_EXOMES_AMR"}},9;
push @{$hashTableAnnoHeader{"GNOMAD_EXOMES_EAS"}},10;
push @{$hashTableAnnoHeader{"GNOMAD_EXOMES_SAS"}},11;
push @{$hashTableAnnoHeader{"GNOMAD_EXOMES_NFE"}},12;
push @{$hashTableAnnoHeader{"GNOMAD_EXOMES_FIN"}},13;

push @{$hashTableAnnoHeader{"GNOMAD_GENOMES_AF"}},14;
push @{$hashTableAnnoHeader{"GNOMAD_GENOMES_AFR"}},15;
push @{$hashTableAnnoHeader{"GNOMAD_GENOMES_AMR"}},16;
push @{$hashTableAnnoHeader{"GNOMAD_GENOMES_EAS"}},17;
push @{$hashTableAnnoHeader{"GNOMAD_GENOMES_SAS"}},18;
push @{$hashTableAnnoHeader{"GNOMAD_GENOMES_NFE"}},19;
push @{$hashTableAnnoHeader{"GNOMAD_GENOMES_FIN"}},20;

push @{$hashTableAnnoHeader{"CONSEQUENCE"}},25;
push @{$hashTableAnnoHeader{"INTRON_VARIANT"}},26;

push @{$hashTableAnnoHeader{"POLYPHEN2_HDIV_pred"}},40;
push @{$hashTableAnnoHeader{"POLYPHEN2_HVAR_pred"}},41;
push @{$hashTableAnnoHeader{"LOF"}},42;

push @{$hashTableAnnoHeader{"CANONICAL"}},45;

push @{$hashTableAnnoHeader{"AF_CTRL"}},46;
push @{$hashTableAnnoHeader{"AF_CASE"}},47;
push @{$hashTableAnnoHeader{"AF_ALL"}},48;
push @{$hashTableAnnoHeader{"MCR_CTRL"}},49;
push @{$hashTableAnnoHeader{"MCR_CASE"}},50;
push @{$hashTableAnnoHeader{"MCR_ALL"}},51;
push @{$hashTableAnnoHeader{"FILTER"}},52;

push @{$hashTableAnnoHeader{"P_MISS"}},61;
push @{$hashTableAnnoHeader{"HWE_CTRL_PVAL"}},66;


open FILTERING, $parameter or die $!;

	while (<FILTERING>) 
	{	chomp($_);
			&getFiltering($_);
	}

close(FILTERING);	

sub getFiltering()
{
	my @t = split(/\s+/, $_);
	$t[0] =~ s/(^\s+|\s+$)//g;
	
	if ( exists $hashTableAnnoHeader{$t[0]}) 
	{
		my $criteria="";
		my $index=0;
		foreach(@t)
		{
			++$index;
			$_ =~ s/(^\s+|\s+$)//g;
			if($_ ne 'or' and $index > 1)
			{
				$criteria = $criteria.$_."%";
			}
		}
		#print "insert key:$t[0]\n";
		#print "criteria:$criteria\n";
		chomp($t[0]);
		chomp($criteria);
		if ((length($t[0]) >0) and (length($criteria) >0))
		{
			push @{$hashTableFiltering{$t[0]}},  $criteria;	
		}
	}
}

open FINALANNOTATION, $annotationResults or die $!;

	while (<FINALANNOTATION>) 
	{
		chomp($_);
		if ($.==1)
		{
			#not to print out heaader
			#print OUTANNOTATION $_;
			next;
		}
		&getFilteredAnnotationFile($_);
	}

close(FINALANNOTATION);	

open GENOTYPINGANNO, $seqVariants or die $!;

	while (<GENOTYPINGANNO>) 
	{
			if ($.==1)
		{
			my @lines = split (/\t/, $_);
			my $size=$#lines;
			$afCol = &getNumberOfAFColumns(\@lines);
			#print "afCol:$afCol\n";
			my $last=$size-$afCol;
			my @samples = @lines[5..$last];
			my $index=0;
			foreach (@samples) 
			{ 
					chomp($_);
					my @sampleid= split (/:/, $_);
					$phenotype[$index]=$sampleid[1];	
					$index++;
			}
		
			
			next;
		}
		&getGenotypingFile($_);
	}

close(GENOTYPINGANNO);	

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


sub getGenotypingFile()
{
	chomp($_);
	my @t = split (/\t/, $_);
	my $chr  = $t[0];
	$chr =~ s/chr//;
	my $pos = $t[1];
	my $ref=$t[2];
	my $alt=$t[3];
	my $transcriptid=$t[4];

	if ( exists $hashTableVariant{$chr."l".$pos."l".$ref."l".$alt."l".$transcriptid} )
	{
		my $size=$#t;
		my $last=$size-$afCol;
		my @samples = @t[5..$last];
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
		my $col_filter=$size-14;
		my $filter=$t[$col_filter];

		my $countNAControl=0;
		my $count0Control=0;   
		my $count1Control=0;
		my $count2Control=0;
		my $countNACase=0;
		my $count0Case=0;
		my $count1Case=0;
		my $count2Case=0;


		my $col_pmiss=$size-5;
		my $col_pmissOR=$size-4;
		my $col_pmisslogOR=$size-3;
		my $col_pmissabslogOR=$size-2;
		my $col_pmisstotal=$size-1;
		my $col_hwectrlP=$size;

		my $pmiss=$t[$col_pmiss];
		my $pmissOR=$t[$col_pmissOR];
		my $pmisslogOR=$t[$col_pmisslogOR];
		my $pmissabslogOR=$t[$col_pmissabslogOR];
		my $pmisstotal=$t[$col_pmisstotal];
		my $hwectrlP=$t[$col_hwectrlP];


		my @hpuVal = @{$hashTableVariant{$chr."l".$pos."l".$ref."l".$alt."l".$transcriptid}};
		my @tsv_annotation_val = split (/\t/, $hpuVal[0]);
		$size=@tsv_annotation_val;

		my $i=0;
		foreach (@samples) 
		{
			chomp($_);
			$totalSample++;
			my $pheno = $phenotype[$i];

			my $homozygousity=$_;

		#	print "homoz:$homozygousity,pheno:$pheno\n";
			if ($pheno ==0) 
			{
				$phenoControl++;
			}
			if ($pheno ==1) 
			{
				$phenoCase++;
			}

			if ($homozygousity eq "0" or $homozygousity eq "1" or $homozygousity eq "2") 
			{
				
				if($pheno == 1)
				{

					if($homozygousity eq "0")
					{
						$count0Case++;
					}
					if($homozygousity eq "1")
					{
						$alleleFreqCountCase =$alleleFreqCountCase+1;
						$count1Case++;
					}
					if($homozygousity eq "2")
					{
						$alleleFreqCountCase =$alleleFreqCountCase+2;
						$count2Case++;
					}
					$totalCase++;
				}
				else
				{
					if($homozygousity eq "0")
					{
						$count0Control++;
					}
					if($homozygousity eq "1")
					{
						$alleleFreqCountControl =$alleleFreqCountControl+1;
						$count1Control++;
					}
					if($homozygousity eq "2")
					{
						$alleleFreqCountControl =$alleleFreqCountControl+2;
						$count2Control++;
					}
					$totalControl++;
				}
			}

			if ($homozygousity eq "NA" and $pheno ==0) 
			{
				$missingSampleControl++;
				$countNAControl++;
			}
			if ($homozygousity eq "NA" and $pheno ==1) 
			{
				$missingSampleCase++;
				$countNACase++;
			}
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

		for(my $f=0; $f<=$last; $f++)
		{
			print OUTGENOTYPING "$t[$f]\t";
		}
		
		my $line="$alleleFreqControl\t$alleleFreqCase\t$alleleFreqForAllSample\t$percentageOfMissingControl\t$percentageOfMissingCase\t$percentageOfMissingAll\t$filter\t$countNAControl\t$count0Control\t$count1Control\t$count2Control\t$countNACase\t$count0Case\t$count1Case\t$count2Case\t$pmiss\t$pmissOR\t$pmisslogOR\t$pmissabslogOR\t$pmisstotal\t$hwectrlP";
		print OUTGENOTYPING "$line\n";

		
		for($i = 0; $i < $size-$afCol; $i++) 
		{
			print OUTANNOTATION "$tsv_annotation_val[$i]\t";
		}
		print OUTANNOTATION "$line\n";
	}
}



sub getFilteredAnnotationFile()
{
	my @t = split (/\t/, $_);
	$t[0] =~ s/chr//;
	my $pos=$t[1];
	my $ref=$t[2];
	my $alt=$t[3];
	my $transcriptid=$t[4];

	
	my $overall_flag=0;
	foreach my $key (keys %hashTableFiltering)
		{
			my $hpukey = $key;
			my @hpuVal = @{$hashTableFiltering{$hpukey}};
			my @parameters =  split (/%/, $hpuVal[0]);
			my @headCol = @{$hashTableAnnoHeader{$hpukey}};
			my $col=$headCol[0];
			my $status =0;
			my $sizeParam = @parameters;
			
			foreach(@parameters)
			{
				chomp($_);
				$_ =~ s/^\s+|\s+$//g;

				my $lt = length($_);
				my $compare =substr($_,0,1);
				my $compare2 =substr($_,0,2);
				#print "here: $_,compare:$compare,col:$col,tcol:$t[$col],key:$hpukey\n";
				if ($compare eq '<')
				{
					
					#print "small t col:$t[$col], criteria:$_\n";
					if ($compare2 eq '<=')
					{
						$_=substr($_,2,$lt);
						chomp($_);
						#print "small t col:$t[$col], criteria:$_\n";
						if($t[$col] <= $_ and $t[$col] ne 'NA')
						{
							$status=1;						
						}
						if($t[$col] > $_ and $t[$col] ne 'NA')
						{
							$status=0;
							last;
						}
					}
					else
					{
						$_=substr($_,1,$lt);
						chomp($_);
						if($t[$col] < $_ and $t[$col] ne 'NA')
						{
							$status=1;						
						}
						if($t[$col] >= $_ and $t[$col] ne 'NA')
						{
							$status=0;
							last;
						}
					}
				}
				if ( $compare eq '>') 
				{
				
					#print "big t col:$t[$col], criteria:$_\n";
					if ($compare2 eq '>=')
					{
						$_=substr($_,2,$lt);
						chomp($_);
						if($t[$col] >= $_ and $t[$col] ne 'NA')
						{
							$status=1;

						}
						if($t[$col] < $_ and $t[$col] ne 'NA')
						{
							$status=0;
							last;
						}
					}
					else
					{
						$_=substr($_,1,$lt);
						chomp($_);
						if($t[$col] > $_ and $t[$col] ne 'NA')
						{
							$status=1;

						}
						if($t[$col] <= $_ and $t[$col] ne 'NA')
						{
							$status=0;
							last;
						}
					}
				}
				
				
				if ($_ =~ "-") 
				{
				
					my @freq = split (/-/, $_);
					
					#print "big t col:$t[$col], criteria:$_\n";
					
					if($t[$col] >= $freq[0] and $t[$col] <= $freq[1]  and $t[$col] ne 'NA')
					{
						$status=1;

					}
					else
					{
						$status=0;
						last;
					}				
				}

				if ($compare ne '<' and $compare ne '>' and $compare2 ne '<=' and $compare2 ne '>=' and ($_ !~ "-") ) 			
				{
					#$var= substr($t[$col], 0, $lt);
					my $var=$t[$col];
					#print "var:$var,here dollar:$_\n";
					chomp($var);
				        
					#print "filtering_parameter:$_\n";
					#print "SNPTYPE read:$var\n";
	
					#if ((length($var) >0) && (length($_) >0) && $var eq $_)
					#if($t[$col] =~ m/$_/)
					if($_ eq "LowQual" or $_ eq "PASS"  )
					{
						if ((length($var) >0) and ($var eq $_))
						{
							$status=1;
						}
					}
					else{
						if ((length($var) >0) and ($var =~ m/^$_/))
						{
							$status=1;
							#print "SNPTYPE status=1 is $var\n";
							#print "parameter compared is $_\n";
						}
						if ((length($var) >0) and ($var !~ m/^$_/) && $sizeParam < 2)
						{
							$status=0;
							last;
						}
					}
				}

				

			}
			
			#print "overall status:$status\n";
				
			if ($status ==1) {
				$overall_flag=1;
			}
			else
			{
				$overall_flag=0;
				last;
			}
		}

	if ($overall_flag == 1) 
	{
	#	print OUTANNOTATION $_;
		push @{$hashTableVariant{$t[0]."l".$t[1]."l".$t[2]."l".$t[3]."l".$t[4]}},  $_;
	}
}

#################################################


exit (0);
