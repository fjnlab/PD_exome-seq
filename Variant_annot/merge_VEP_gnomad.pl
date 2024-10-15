#!/usr/bin/perl
## Script here merges VEP annot + ANNOVAR gnomad *_dropped => variant annotation list
########################################################################################

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;


#open (INPUT1, "$ARGV[0]") or die "Cannot open INPUT1\n";    ##ie VEP filtered output
#open (INPUT2, "$ARGV[1]") or die "Cannot open INPUT2\n";    ##ie ANNOVAR gnomad exome *_dropped output
#open (INPUT3, "$ARGV[2]") or die "Cannot open INPUT3\n";    ##ie ANNOVAR gnomad genome *_dropped output
#open (OUTPUT, ">$ARGV[3]") or die "Cannot open OUTPUT\n";


my $vep = $ARGV[0];     ##ie VEP filtered output
my $gnomad_exome_dropped = $ARGV[1];    ##ie ANNOVAR gnomad exome *_dropped output
my $gnomad_genome_dropped = $ARGV[2];   ##ie ANNOVAR gnomad genome *_dropped output
my $outfile = "variantannot.tsv";

our ($help);
GetOptions ( 'help|h'=>\$help) or pod2usage ();
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);

die "Usage: merge_VEP_gnomad.pl VEP.out.filtered variants.avlist.hg19_gnomad211_exome_dropped variants.avlist.hg19_gnomad211_genome_dropped\n" if @ARGV < 3;
die "$vep file not exists\n" if (!(-f $vep));
die "$gnomad_exome_dropped file not exists\n" if (!(-f $gnomad_exome_dropped));
die "$gnomad_genome_dropped file not exists\n" if (!(-f $gnomad_genome_dropped));



open (INPUT1, $vep) || die "cannot open $vep\n";
open (INPUT2, $gnomad_exome_dropped) || die "cannot open $gnomad_exome_dropped\n";
open (INPUT3, $gnomad_genome_dropped) || die "cannot open $gnomad_genome_dropped\n";
open (OUTPUT, ">$outfile") || die "can't open $outfile\n";


my %hashTableDBGNOMAD_EXOME;
my %hashTableDBGNOMAD_GENOME;


## Reads in ANNOVAR gnomad exome file, store in hash
while(my $line1=<INPUT2>){
    chomp($line1);
    my @t = split(/\t/, $line1);

    my $key_gnomadexome=$t[10];
    $key_gnomadexome =~ s/chr//;

    ##not using all the gnomad_exomes MAF of annovar
    my @k = split(",", $t[1]);
    my $gnomad_exomes_AF=$k[0];
    my $gnomad_exomes_AFR_AF=$k[5];
    my $gnomad_exomes_AMR_AF=$k[7];
    my $gnomad_exomes_EAS_AF=$k[8];
    my $gnomad_exomes_SAS_AF=$k[6];
    my $gnomad_exomes_NFE_AF=$k[9];
    my $gnomad_exomes_FIN_AF=$k[10];

    push @{$hashTableDBGNOMAD_EXOME{$key_gnomadexome}},  $gnomad_exomes_AF;
    push @{$hashTableDBGNOMAD_EXOME{$key_gnomadexome}},  $gnomad_exomes_AFR_AF;
    push @{$hashTableDBGNOMAD_EXOME{$key_gnomadexome}},  $gnomad_exomes_AMR_AF;
    push @{$hashTableDBGNOMAD_EXOME{$key_gnomadexome}},  $gnomad_exomes_EAS_AF;
    push @{$hashTableDBGNOMAD_EXOME{$key_gnomadexome}},  $gnomad_exomes_SAS_AF;
    push @{$hashTableDBGNOMAD_EXOME{$key_gnomadexome}},  $gnomad_exomes_NFE_AF;
    push @{$hashTableDBGNOMAD_EXOME{$key_gnomadexome}},  $gnomad_exomes_FIN_AF;

}
close(INPUT2);



## Reads in ANNOVAR gnomad genome file, store in hash
while (my $line2=<INPUT3>){
    chomp($line2);
    my @t = split(/\t/, $line2);

    my $key_gnomadgenome=$t[10];
    $key_gnomadgenome =~ s/chr//;

    ##not using all the gnomad_exomes MAF of annovar
    my @k = split(",", $t[1]);
    my $gnomad_genomes_AF=$k[0];
    my $gnomad_genomes_AFR_AF=$k[5];
    my $gnomad_genomes_AMR_AF=$k[7];
    my $gnomad_genomes_EAS_AF=$k[8];
    my $gnomad_genomes_SAS_AF=$k[6];
    my $gnomad_genomes_NFE_AF=$k[9];
    my $gnomad_genomes_FIN_AF=$k[10];

    push @{$hashTableDBGNOMAD_GENOME{$key_gnomadgenome}},  $gnomad_genomes_AF;
    push @{$hashTableDBGNOMAD_GENOME{$key_gnomadgenome}},  $gnomad_genomes_AFR_AF;
    push @{$hashTableDBGNOMAD_GENOME{$key_gnomadgenome}},  $gnomad_genomes_AMR_AF;
    push @{$hashTableDBGNOMAD_GENOME{$key_gnomadgenome}},  $gnomad_genomes_EAS_AF;
    push @{$hashTableDBGNOMAD_GENOME{$key_gnomadgenome}},  $gnomad_genomes_SAS_AF;
    push @{$hashTableDBGNOMAD_GENOME{$key_gnomadgenome}},  $gnomad_genomes_NFE_AF;
    push @{$hashTableDBGNOMAD_GENOME{$key_gnomadgenome}},  $gnomad_genomes_FIN_AF;
}
close(INPUT3);



print OUTPUT "#Chr\tPosition\tRef\tAlt\tTranscriptID\tGeneID\tGeneName\tgnomAD2.1.1_exomes_AF\tgnomAD2.1.1_exomes_AFR\tgnomAD2.1.1_exomes_AMR\tgnomAD2.1.1_exomes_EAS\tgnomAD2.1.1_exomes_SAS\tgnomAD2.1.1_exomes_NFE\tgnomAD2.1.1_exomes_FIN\tgnomAD2.1.1_genomes_AF\tgnomAD2.1.1_genomes_AFR\tgnomAD2.1.1_genomes_AMR\tgnomAD2.1.1_genomes_EAS\tgnomAD2.1.1_genomes_SAS\tgnomAD2.1.1_genomes_NFE\tgnomAD2.1.1_genomes_FIN\t#Uploaded_variation\tLocation\tAllele\tTranscript\tConsequence\tIntron_variant\tExon_number\tIntron_number\tProtein_pos\tAmino_acids\tHGVSc\tHGVSp\tSnpID\tImpact\tStrand\tVariant_class\tBiotype\tPolyphen2_Transcriptid\tUniprot_acc_Polyphen2\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_pred\tLoF\tLoF_filter\tLoF_info\tCANONICAL\n";


## Reads in VEP annot list, merge with corresponding gnomad exome, genome pop AF
while(my $line3=<INPUT1>) {
	chomp($line3);

	if ($.==1) {
		next;
	}

	my @t = split(/\t/, $line3);

	my $chr = $t[0];
   	my $pos = $t[1];
    	my $ref = $t[2];
    	my $alt = $t[3];

	my $key=$chr.":".$pos.":".$ref.":".$alt;


	my @part1=@t[0 .. 6];
	my @part2=@t[7 .. 31];

	
	##gnomad exomes and genomes MAF; r2.1.1
    	my $gnomadexomes_AF="NA";
    	my $gnomadexomes_AFR_MAF="NA";
    	my $gnomadexomes_AMR_MAF="NA";
    	my $gnomadexomes_EAS_MAF="NA";
    	my $gnomadexomes_SAS_MAF="NA";
    	my $gnomadexomes_NFE_MAF="NA";
    	my $gnomadexomes_FIN_MAF="NA";

    	my $gnomadgenomes_AF="NA";
    	my $gnomadgenomes_AFR_MAF="NA";
    	my $gnomadgenomes_AMR_MAF="NA";
    	my $gnomadgenomes_EAS_MAF="NA";
    	my $gnomadgenomes_SAS_MAF="NA";
    	my $gnomadgenomes_NFE_MAF="NA";
    	my $gnomadgenomes_FIN_MAF="NA";


	if ( exists $hashTableDBGNOMAD_EXOME{$key} ) {
	       my @hpuGNOMADEXOMEVal = @{$hashTableDBGNOMAD_EXOME{$key}};
	       $gnomadexomes_AF=$hpuGNOMADEXOMEVal[0];
	       $gnomadexomes_AFR_MAF=$hpuGNOMADEXOMEVal[1];
	       $gnomadexomes_AMR_MAF=$hpuGNOMADEXOMEVal[2];
	       $gnomadexomes_EAS_MAF=$hpuGNOMADEXOMEVal[3];
	       $gnomadexomes_SAS_MAF=$hpuGNOMADEXOMEVal[4];
	       $gnomadexomes_NFE_MAF=$hpuGNOMADEXOMEVal[5];
	       $gnomadexomes_FIN_MAF=$hpuGNOMADEXOMEVal[6];
	}

	if ( exists $hashTableDBGNOMAD_GENOME{$key} ) {
               my @hpuGNOMADGENOMEVal = @{$hashTableDBGNOMAD_GENOME{$key}};
               $gnomadgenomes_AF=$hpuGNOMADGENOMEVal[0];
               $gnomadgenomes_AFR_MAF=$hpuGNOMADGENOMEVal[1];
               $gnomadgenomes_AMR_MAF=$hpuGNOMADGENOMEVal[2];
               $gnomadgenomes_EAS_MAF=$hpuGNOMADGENOMEVal[3];
	       $gnomadgenomes_SAS_MAF=$hpuGNOMADGENOMEVal[4];
	       $gnomadgenomes_NFE_MAF=$hpuGNOMADGENOMEVal[5];
	       $gnomadgenomes_FIN_MAF=$hpuGNOMADGENOMEVal[6];
        }


	if ($gnomadexomes_AF =~/^\.$/) {$gnomadexomes_AF="NA";}
	if ($gnomadexomes_AFR_MAF =~/^\.$/) {$gnomadexomes_AFR_MAF="NA";}
	if ($gnomadexomes_AMR_MAF =~/^\.$/) {$gnomadexomes_AMR_MAF="NA";}
	if ($gnomadexomes_EAS_MAF =~/^\.$/) {$gnomadexomes_EAS_MAF="NA";}
	if ($gnomadexomes_SAS_MAF =~/^\.$/) {$gnomadexomes_SAS_MAF="NA";}
	if ($gnomadexomes_NFE_MAF =~/^\.$/) {$gnomadexomes_NFE_MAF="NA";}
	if ($gnomadexomes_FIN_MAF =~/^\.$/) {$gnomadexomes_FIN_MAF="NA";}
	if ($gnomadgenomes_AF =~/^\.$/) {$gnomadgenomes_AF="NA";}
	if ($gnomadgenomes_AFR_MAF =~/^\.$/) {$gnomadgenomes_AFR_MAF="NA";}
	if ($gnomadgenomes_AMR_MAF =~/^\.$/) {$gnomadgenomes_AMR_MAF="NA";}
	if ($gnomadgenomes_EAS_MAF =~/^\.$/) {$gnomadgenomes_EAS_MAF="NA";}
	if ($gnomadgenomes_SAS_MAF =~/^\.$/) {$gnomadgenomes_SAS_MAF="NA";}
	if ($gnomadgenomes_NFE_MAF =~/^\.$/) {$gnomadgenomes_NFE_MAF="NA";}
	if ($gnomadgenomes_FIN_MAF =~/^\.$/) {$gnomadgenomes_FIN_MAF="NA";}



	my $part1_line=join("\t",@part1);
	my $part2_line=join("\t",@part2);


	print OUTPUT "$part1_line\t$gnomadexomes_AF\t$gnomadexomes_AFR_MAF\t$gnomadexomes_AMR_MAF\t$gnomadexomes_EAS_MAF\t$gnomadexomes_SAS_MAF\t$gnomadexomes_NFE_MAF\t$gnomadexomes_FIN_MAF\t$gnomadgenomes_AF\t$gnomadgenomes_AFR_MAF\t$gnomadgenomes_AMR_MAF\t$gnomadgenomes_EAS_MAF\t$gnomadgenomes_SAS_MAF\t$gnomadgenomes_NFE_MAF\t$gnomadgenomes_FIN_MAF\t$part2_line\n";

}
close(INPUT1);
close(OUTPUT);



=head1 SYNOPSIS

 perl merge_VEP_gnomad.pl [arguments] <VEP.out.filtered> <ANNOVAR *gnomad211_exome_dropped> <ANNOVAR *gnomad211_exome_dropped> [output file] <variantannot.tsv>

 Optional arguments:
 	-h, --help			print help message

 Function: To merge VEP annot with ANNOVAR gnomad *_dropped => variant annotation list

 Example: perl merge_VEP_gnomad.pl VEP.out.filtered variants.avlist.hg19_gnomad211_exome_dropped variants.avlist.hg19_gnomad211_genome_dropped



=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

