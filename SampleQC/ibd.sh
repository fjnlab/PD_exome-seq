#!/bin/bash
source ~/.bash_profile

BFILE_PREFIX='exome_biallelic_snps'
OUT_PREFIX='exome_samples_removed'


##remove samples failing upstream sampleQC
plink --noweb --bfile ${BFILE_PREFIX} --remove samplesfail_QC1to4.list --make-bed --out ${OUT_PREFIX}


## Prune first ##
## REMOVE SNPs IN LD FOR PCA and use only common SNPs; ie only include SNPs with MAF >= 0.01
plink --noweb --bfile ${OUT_PREFIX} --maf 0.01 --indep-pairwise 500 50 0.1 --out ${OUT_PREFIX}.common-prune
plink --noweb --bfile ${OUT_PREFIX} --maf 0.01 --extract ${OUT_PREFIX}.common-prune.prune.in --make-bed --out ${OUT_PREFIX}.common-pruned



##Sieve out CONSERVED REGION SNPS (CRS)
bfile_prefix='exome_samples_removed.common-pruned'

plink --noweb --bfile ${bfile_prefix} --chr 6 --from-bp 24892021 --to-bp 36892022 --write-snplist --out temp1
plink --noweb --bfile ${bfile_prefix} --chr 11 --from-bp 45043424 --to-bp 57243424 --write-snplist --out temp2
plink --noweb --bfile ${bfile_prefix} --chr 5 --from-bp 43964243 --to-bp 51464243 --write-snplist --out temp3
plink --noweb --bfile ${bfile_prefix} --chr 8 --from-bp 10000 --to-bp 12655629 --write-snplist --out temp4
plink --noweb --bfile ${bfile_prefix} --chr 11 --from-bp 84322352 --to-bp 86322352 --write-snplist --out temp5

cat temp*snplist > snps-conserved-correlated.txt

plink --noweb --bfile ${bfile_prefix} --exclude snps-conserved-correlated.txt --make-bed --out ${bfile_prefix}-recode



##remove non-autosomal snps
##also remove snps with high missing genotype rate ie >5%
cat ${bfile_prefix}-recode.bim | awk '{if($1~/^23$/)print $2}' > chrX_snps.list
cat ${bfile_prefix}-recode.bim | awk '{if($1~/^24$/)print $2}' > chrY_snps.list
cat chrX_snps.list chrY_snps.list > non-autosomal_snps.txt

plink --noweb --bfile ${bfile_prefix}-recode --geno 0.05 --exclude non-autosomal_snps.txt --make-bed --out ${bfile_prefix}-recode-qced



##run PLINK Identity-By-Descent
plink --noweb --file ${bfile_prefix}-recode-qced --genome gz full --read-freq ${bfile_prefix}-recode-qced-hwe.hwe --out ${bfile_prefix}-recode-qced.genome
unpigz -p 10 ${bfile_prefix}-recode-qced.genome.gz
plinkgenome2sar ${bfile_prefix}-recode-qced.genome

