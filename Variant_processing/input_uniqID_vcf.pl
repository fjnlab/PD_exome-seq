#!/usr/bin/perl
## here inputs unique id for variants with unknown snpid (marked as "." in 3rd col)
## 
## formated id to chr:pos:ref:alt
############################################################################################

use strict;
use warnings;

open(VCF_IN,"$ARGV[0]") or die ("Cannot open VCF_IN file\n");  ##read in VCF
open(VCF_OUT,">$ARGV[1]") or die ("Cannot open OUT file\n");  ##edited VCF

my $count=0;

while (my $line1=<VCF_IN>) {
     chomp($line1);

     ##print out VCF headers
     if ($line1=~/^#/){
         print VCF_OUT "$line1\n";
     }
     else {
         my @line1_edit=split(/\t/,$line1);
	 my $chr=$line1_edit[0];
	 my $pos=$line1_edit[1];
	 my $ref=$line1_edit[3];
	 my $alt=$line1_edit[4];

	 my $dbsnpid="$chr:$pos:$ref:$alt";

         ##for variants with unknown snpid (marked as "." in 3rd col) - create unique ID
         if($line1_edit[2] =~/^\./) {
            $count++;
            my $tmp=$dbsnpid;
            $line1_edit[2] =~ s/^\./$tmp/;

            print VCF_OUT join("\t",@line1_edit),"\n";            
         }
         ##else do nothing, just print the variant information
         else {
            print VCF_OUT $line1,"\n";
         }
     }
}
close(VCF_IN);
close(VCF_OUT);

