#!/bin/bash

###################################################
## for plink gender imputation only
## here Ycount used
###################################################

## reference: https://www.cog-genomics.org/plink/1.9/basic_stats
## reference: https://groups.google.com/forum/#!topic/plink2-users/28LESfNj64A


source ~/.bash_profile

## 1. Convert VCF to PLINK format AND for sex chromosomes only
VCF_IN='exome_biallelic_snps.vcf.gz'
BFILE_PREFIX='exome_biallelic_snps'
BINARY_OUT=${BFILE_PREFIX}.chrXY-only

plink --vcf ${VCF_IN} --make-bed --double-id --out ${BINARY_OUT} --chr 23-24



## 2. split-x; X chromosome's pseudo-autosomal region numeric coded as chr 25
plink --bfile ${BFILE_PREFIX}.chrXY-only --split-x b37 --make-bed --out ${BFILE_PREFIX}.chrXY-only.split-x --allow-no-sex


## 3. LD prune data 
plink --bfile ${BFILE_PREFIX}.chrXY-only.split-x --indep-pairphase 20000 2000 0.5 --chr 23-24 --allow-no-sex
plink --bfile ${BFILE_PREFIX}.chrXY-only.split-x --extract plink.prune.in --make-bed --out ${BED_PREFIX}.chrXY-only.split-x.ld_pruned_xy --allow-no-sex


## 4. gender checks; using Ycount
plink --bfile ${BFILE_PREFIX}.chrXY-only.split-x.ld_pruned_xy --check-sex ycount --out ${BFILE_PREFIX}.chrXY-only.split-x.ld_pruned_xy.ycount


## 5. use R to construct histogram; find gap between 2 distribution first BEFORE running the next PLINK command


## 6. Use --impute-sex with the gap position, ie 0.6
## then look at fam file for imputed gender, Sex code ('1' = male, '2' = female, '0' = unknown)
plink --bfile ${BFILE_PREFIX}.chrXY-only.split-x.ld_pruned_xy --impute-sex y-only 0.6 0.6 --make-bed --out ${BFILE_PREFIX}.chrXY-only.split-x.ld_pruned_xy.ycount.imputed_sex


