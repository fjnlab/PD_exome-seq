#!/bin/bash
## using vcftools to filter for DP8GQ20 - to comply with variant annotation list
## individual sites GT converted to ./. if fails DP8GQ20 criteria
## reference: https://www.biostars.org/p/194918
##
## purpose: genotype level QC
######################################################################################
source ~/.bash_profile

VCF_IN='variants.vcf.gz'
VCF_OUT_PREFIX='variants.filtered_DPGQ'

## vcftools cannot generate out vcf.gz file
## output file <file_prefix>.recode.vcf
vcftools --gzvcf ${VCF_IN} --minGQ 20 --minDP 8 --recode --recode-INFO-all --out ${VCF_OUT_PREFIX}
mv ${VCF_OUT_PREFIX}.recode.vcf ${VCF_OUT_PREFIX}.vcf

bgzip --threads 10 -c ${VCF_OUT_PREFIX}.vcf > ${VCF_OUT_PREFIX}.vcf.gz
tabix -p vcf ${VCF_OUT_PREFIX}.vcf.gz

