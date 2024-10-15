#!/usr/bin/perl
## script here extracts specific (field) information from VEP annot table output
##################################################################################################################

use strict;
use warnings;

##USAGE: perl VEPannot_extract.pl [input] <*DP8GQ20_filtered.uniqID.vep.table.out> [output] <*DP8GQ20_filtered.uniqID.vep.table.filtered.out>
##example: perl VEPannot_extract.pl variants.DP8GQ20_filtered.uniqID.vep.table.out variants.DP8GQ20_filtered.uniqID.vep.table.filtered.out

open (INPUT1, "$ARGV[0]") or die "Cannot open INPUT1\n"; ##VEP input (tab-delimited)
open (OUTPUT,">$ARGV[1]") or die "Cannot open OUTPUT1\n";  ##extracted output


##format 
##Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation
#EXON
#INTRON
#HGVSc
#HGVSp
#IMPACT
#STRAND
#BIOTYPE
#VARIANT_CLASS
#Ensembl_transcriptid
#Uniprot_acc_Polyphen2
#Polyphen2_HDIV_pred
#Polyphen2_HVAR_pred
#LoF
#LoF_filter
#LoF_info
#CANONICAL


print OUTPUT "#Chr\tPosition\tRef\tAlt\tTranscriptID\tGeneID\tGeneName\t#Uploaded_variation\tLocation\tAllele\tTranscript\tConsequence\tIntron_variant\tExon_number\tIntron_number\tProtein_pos\tAmino_acids\tHGVSc\tHGVSp\tSnpID\tImpact\tStrand\tVariant_class\tBiotype\tPolyphen2_Transcriptid\tUniprot_acc_Polyphen2\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_pred\tLoF\tLoF_filter\tLoF_info\tCANONICAL\n";

while (my $line1=<INPUT1>) {
	chomp($line1);

	my $variation="NA";
	my $location="NA";
	my $allele="NA";
	my $ensg="NA";
	my $enst="NA";
	my $consequence="NA";
	my $intron_variant="NA";
	my $exon_number="NA";
	my $intron_number="NA";

	my $aa_pos="NA";
	my $aa_change="NA";
	my $genename="NA";

	my $HGVSc="NA";
	my $HGVSp="NA";
	my $snpid="NA";
	my $impact="NA";
	my $strand="NA";
	my $variant_class="NA";
	my $biotype="NA";
	
	my $polyphen2_transcriptid="NA";
	my $uniprot_acc_polyphen2="NA";
	my $polyphen2_hdiv="NA";
	my $polyphen2_hvar="NA";
	
	my $LoF="NA";
	my $LoF_filter="NA";
	my $LoF_info="NA";
	
	my $canonical="NA";

	my @t=split(/\t/,$line1);

	if ($line1=~/^#/) {
		next;
	}


	##extracting variant info from #Uploaded_variation col
	my @var_array=split(/\:/,$t[0]);
	my $chr=$var_array[0];
	my $pos=$var_array[1];
	my $ref=$var_array[2];
	my $alt=$var_array[3];
	#print "$chr\t$pos\t$ref\t$alt\t$t[22]\n";

	$variation=$t[0];
	$location=$t[1];
	$allele=$t[2];
	$enst=$t[4];
	$ensg=$t[3];
	$genename=$t[18];

	$consequence=$t[6];
	$exon_number=$t[23];
	$intron_number=$t[24];
	$aa_pos=$t[9];
	$aa_change=$t[10];

	$HGVSc=$t[25];
	$HGVSp=$t[26];
	$snpid=$t[12];
	$impact=$t[13];
	$strand=$t[15];
	$variant_class=$t[17];
	$biotype=$t[21];

	$polyphen2_transcriptid=$t[31];
	$uniprot_acc_polyphen2=$t[34];
	$polyphen2_hdiv=$t[32];
	$polyphen2_hvar=$t[33];

	$LoF=$t[35];
	$LoF_filter=$t[36];
	$LoF_info=$t[38];

	$canonical=$t[22];

	my @extract;


	##edits exon number annotation
	if($exon_number=~/^-$/) {
		$exon_number="NA";
	}
	else {
		$exon_number=~s/\//\/\//;    ##replace 3/3 to 3//3, else dated when list exported to excel
	}


	##edits intron number annotation
	if($intron_number=~/^-$/) {
		$intron_number="NA";
	}
	else{
		$intron_number=~s/\//\/\//;    ##replace 3/3 to 3//3, else dated when list exported to excel
	}


	##identifies if variant is in intron
        if( ($intron_number!~/^NA$/) && ($exon_number=~/^NA$/) ){
		$intron_variant="YES";
        }
        else {
                $intron_variant="NO";
        }	


	##re-categorise non-canonical transcript annotation
	if($canonical=~/^-$/) {
		$canonical="NO";
	}
	else{;}  


	##edit synonymous HGVSp annotations; remove "%3D" to "="
        ##not sure if is issue with v104 cache database
        if( ($consequence=~/synonymous/) && ($HGVSp=~/\%3D$/) ) {
                $HGVSp=~ s/\%3D/=/;
        }


	##edits annotations with '-' to NA
	if($LoF=~/^-$/) {
		$LoF="NA";
	}
	else{;}

	if($polyphen2_transcriptid=~/^-$/) {
                $polyphen2_transcriptid="NA";
        }
        else{;}


	if($uniprot_acc_polyphen2=~/^-$/) {
                $uniprot_acc_polyphen2="NA";
        }
        else{;}


	if($polyphen2_hdiv=~/^-$/) {
                $polyphen2_hdiv="NA";
        }
        else{;}


	if($polyphen2_hvar=~/^-$/) {
                $polyphen2_hvar="NA";
        }
        else{;}


	##push variables in specific order, exclude up/downstream consequences annotations
	if( ($consequence=~/upstream_gene_variant/) || ($consequence=~/downstream_gene_variant/) ) {
		next;
	}

	else{
		push(@extract,$chr);
		push(@extract,$pos);
		push(@extract,$ref);
		push(@extract,$alt);

		push(@extract,$enst);
		push(@extract,$ensg);
		push(@extract,$genename);

		push(@extract,$variation);		
		push(@extract,$location);
		push(@extract,$allele);
		push(@extract,$enst);
		push(@extract,$consequence);
		push(@extract,$intron_variant);
		push(@extract,$exon_number);
		push(@extract,$intron_number);
		push(@extract,$aa_pos);
		push(@extract,$aa_change);
		push(@extract,$HGVSc);
		push(@extract,$HGVSp);
		push(@extract,$snpid);
		push(@extract,$impact);
		push(@extract,$strand);
		push(@extract,$variant_class);
		push(@extract,$biotype);
		push(@extract,$polyphen2_transcriptid);
		push(@extract,$uniprot_acc_polyphen2);
		push(@extract,$polyphen2_hdiv);
		push(@extract,$polyphen2_hvar);
		push(@extract,$LoF);
		push(@extract,$LoF_filter);
		push(@extract,$LoF_info);
		push(@extract,$canonical);

		print OUTPUT join("\t",@extract),"\n";
	}
}
close(INPUT1);
close(OUTPUT);
