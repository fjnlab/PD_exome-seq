#!/usr/bin/perl
## script here merges variantannot.tsv with samplegenotype_sorted_pmiss_hwe
############################################################################

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;


#open (INPUT1,"$ARGV[0]") or die "Cannot open INPUT1\n";   ##ie variantannot.tsv
#open (INPUT2,"$ARGV[1]") or die "Cannot open INPUT2\n";   ##ie samplegenotype_sorted_pmiss_hwe
#open (OUTPUT1,">$ARGV[2]") or die "Cannot open OUTPUT1\n"; ##ie variantannot_samplegeno.tsv  
#open (OUTPUT2,">$ARGV[3]") or die "Cannot open OUTPUT2\n"; ##ie samplegeno.tsv, with variant tid col added. Required for variant filtering later


my $variantannot=$ARGV[0];      ##ie variantannot.tsv
my $samplegeno=$ARGV[1];        ##ie samplegenotype_sorted_pmiss_hwe
my $outfile1 = "variantannot_samplegeno.tsv";
my $outfile2 = "samplegeno.tsv";       ##ie sample genotype with variant tid col added. Required for variant filtering later

our ($help);
GetOptions ( 'help|h'=>\$help) or pod2usage ();
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);

die "Usage: merge_annot_geno-pmiss-hwe.pl variantannot.tsv samplegenotype_sorted_pmiss_hwe\n" if @ARGV < 2;
die "$variantannot file not exists\n" if (!(-f $variantannot));
die "$samplegeno file not exists\n" if (!(-f $samplegeno));



open (INPUT1, $variantannot) || die "cannot open $variantannot\n";
open (INPUT2, $samplegeno) || die "cannot open $samplegeno\n";
open (OUTPUT1, ">$outfile1") || die "can't open $outfile1\n";
open (OUTPUT2, ">$outfile2") || die "can't open $outfile2\n";


my %hash_geno;
my %hash_genometric;
my @header_geno;
my @header_genometric;


##reads in samplegenotype_sorted, store in hash
while(my $line1=<INPUT2>) {
	chomp($line1);

	if ($.==1) {
		@header_geno = split(/\t/, $line1);
		@header_genometric = @header_geno[-21 .. -1];
		#print join("\t",@header_genometric),"\n";
		next;
	}

	my @t=split(/\t/, $line1);
	my $chr=$t[0];
	my $pos=$t[1];
	my $ref=$t[2];
	my $alt=$t[3];
	my $key="$chr:$pos:$ref:$alt";


	my @samplegeno = @t[4 .. $#t];
	my $samplegeno_line=join("\t",@samplegeno);
	#print "$samplegeno_line","\n";


	my @genometric = @t[-21 .. -1];
	my $genometric_line=join("\t",@genometric);

	if(!exists $hash_geno{$key}) {
		$hash_geno{$key}=$samplegeno_line;
	}
	else{;}


	if(!exists $hash_genometric{$key}) {
		$hash_genometric{$key}=$genometric_line;
	}
	else{;}

}
close(INPUT2);


##header re-arrangements
my $header_variantannot=<INPUT1>;
chomp($header_variantannot);
my $h_genometric=join("\t",@header_genometric);
my @newheader_geno2=@header_geno[4 .. $#header_geno];
my $h_newgeno2=join("\t",@newheader_geno2);
print OUTPUT1 "$header_variantannot\t$h_genometric\tChr\tPos\tRef\tAlt\tTranscriptID\t$h_newgeno2\n";
print OUTPUT2 "#Chr\tPos\tRef\tAlt\tTranscriptID\t$h_newgeno2\n";


##reads in variantannot.tsv, does the match & merge corresponding samplegeno
while(my $line2=<INPUT1>) {
	chomp($line2);

	my @k=split(/\t/, $line2);
	my $chr=$k[0];
        my $pos=$k[1];
        my $ref=$k[2];
        my $alt=$k[3];
        my $key="$chr:$pos:$ref:$alt";

	my $transcriptID=$k[4];


	##if matches $hash_geno{$key}, will also exists match for $hash_genometric{$key}
	if(exists $hash_geno{$key}) {
		#print $hash_geno{$key},"\n";
		print OUTPUT1 "$line2\t$hash_genometric{$key}\t$chr\t$pos\t$ref\t$alt\t$transcriptID\t$hash_geno{$key}\n";
		print OUTPUT2 "$chr\t$pos\t$ref\t$alt\t$transcriptID\t$hash_geno{$key}\n";
	}
	else{;}

}
close(INPUT1);
close(OUTPUT1);
close(OUTPUT2);



=head1 SYNOPSIS

 perl merge_annot_geno-pmiss-hwe.pl [arguments] <variantannot.tsv> <samplegenotype_sorted_pmiss_hwe> [output files] <variantannot_samplegeno.tsv> <samplegeno.tsv>

 Optional arguments:
 	-h, --help			print help message

 Function: To merge variantannot.tsv with samplegenotype_sorted_pmiss_hwe

 Example: perl merge_annot_geno-pmiss-hwe.pl variantannot.tsv samplegenotype_sorted_pmiss_hwe variantannot_samplegeno.tsv samplegeno.tsv



=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

