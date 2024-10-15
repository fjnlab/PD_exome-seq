#!/usr/perl
## Script here extracts the benign filtered dataset
## Rare pathogenic set/ rare deleterious set
############################################################

use strict;
no warnings;
use Pod::Usage;
use Getopt::Long;

my $rare_variantannot=$ARGV[0];      ##ie rare pathogenic set; variantannot_samplegeno_final_all_results_filtered.tsv
my $del_variantannot=$ARGV[1];        ##ie rare deleterious set; variantannot_samplegeno_final_all_results_filtered.tsv
my $outdir = $ARGV[2];    ##benign dataset output dir
my $outAnnotations = $outdir . "variantannot_samplegeno_final_all_results_filtered.tsv";   ##benign dataset FULL annot+samplegeno,
my $outGenotyping = $outdir . "samplegeno_filtered.tsv";	##benign dataset samplegeno


our ($help);
GetOptions ( 'help|h'=>\$help) or pod2usage ();
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);

die "Usage: perl extract_benign.pl ./rare_qced/variantannot_samplegeno_final_all_results_filtered.tsv ./rare_qced_del/variantannot_samplegeno_final_all_results_filtered.tsv ./rare_qced_benign/ \n" if @ARGV < 3;
die "$rare_variantannot file not exists\n" if (!(-f $rare_variantannot));
die "$del_variantannot file not exists\n" if (!(-f $del_variantannot));
die "$outdir directory not exists\n" if (!(-d $outdir));



open(INPUT1, $rare_variantannot) || die "can't open $rare_variantannot";
open(INPUT2, $del_variantannot) ||  die "can't open $del_variantannot";
open (OUTANNOTATION, ">$outAnnotations") || die "can't open $outAnnotations";
open (OUTGENOTYPING, ">$outGenotyping") || die "can't open $outGenotyping";


my %hash_del;
my $count_del=0;
my $count_rare=0;
my $count_unmatch=0;
my $count_benign=0;


##reads in rare deleterious set, store variant in hash
##key: chr:pos:ref:alt:enstid
while(my $line2=<INPUT2>) {
	chomp($line2);
	my @t=split(/\t/,$line2);

	if($.==1) {
		#print $line2,"\n";
		next;
	}	

	$count_del++;
	my $chr=$t[0];
	$chr=~s/^chr//;
	my $pos=$t[1];
	my $ref=$t[2];
	my $alt=$t[3];
	my $tid=$t[4];
	my $key="$chr:$pos:$ref:$alt:$tid";

	if(!exists $hash_del{$key}) {
		$hash_del{$key}=1;
	}

}
close(INPUT2);


##reads in rare pathogenic set, prints variant not in rare_del hash
##key: chr:pos:ref:alt:enstid
while(my $line1=<INPUT1>) {
	chomp($line1);
	my @s=split(/\t/,$line1);

	my $last=$#s;
	#print "$last\n";

	if($.==1) {
		print OUTANNOTATION $line1,"\n";
		my @header_geno=@s[67 .. $last];
		print OUTGENOTYPING "#",join("\t",@header_geno),"\n";
		next;
	}

	$count_rare++;
	my $chr=$s[0];
        $chr=~s/^chr//;
        my $pos=$s[1];
        my $ref=$s[2];
        my $alt=$s[3];
        my $tid=$s[4];
        my $key="$chr:$pos:$ref:$alt:$tid";

	if(!exists $hash_del{$key}) {
		$count_unmatch++;
		print OUTANNOTATION $line1,"\n";
		my @variant_geno=@s[67 .. $last];
                print OUTGENOTYPING "#",join("\t",@variant_geno),"\n";
	}
	else{;}

}
close(INPUT1);
close(OUTANNOTATION);
close(OUTGENOTYPING);

$count_benign=$count_rare - $count_del;


## checks
print "No. of variants in rare pathogenic set:$count_rare\n";
print "No. of variants in rare deleterious set:$count_del\n";
print "No. of variants to retrieve in rare benign set:$count_benign\n";
print "No. of rare benign variants retrieved in list:$count_unmatch\n";


=head1 SYNOPSIS

 perl extract_benign.pl [arguments] <rare pathogenic variantannot_samplegeno_final_all_results_filtered.tsv> <rare deleterious variantannot_samplegeno_final_all_results_filtered.tsv> [output directory] <out dir/>

 Optional arguments:
 	-h, --help			print help message

 Function: To extract rare benign filtered dataset (rare pathogenic - rare deleterious)

 Example: perl extract_benign.pl ./rare_qced/variantannot_samplegeno_final_all_results_filtered.tsv ./rare_qced_del/variantannot_samplegeno_final_all_results_filtered.tsv ./rare_qced_benign/ 



=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

