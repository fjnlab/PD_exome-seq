#!/bin/bash
source ~/.bash_profile

BFILE_PREFIX='exome_biallelic_snps'

##Remove snps with high missing genotype rate ie >5%, prior to IBD step
plink --noweb --bfile ${bfile_prefix} --mind 0.05 --make-bed --out ${bfile_prefix}-mind0.05

