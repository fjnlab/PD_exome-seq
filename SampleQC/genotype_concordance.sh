#!/bin/bash
source ~/.bash_profile

FILE1_PREFIX='gwas'
FILE2_PREFIX='exome'
MERGE_OUT_PREFIX='exome_gwas_merge'


##merge 2 PLINK formatted dataset
plink --file ${FILE1_PREFIX} --bmerge ${FILE2_PREFIX} --make-bed --double-id --out ${MERGE_OUT_PREFIX} --biallelic-only strict --allow-no-sex



## Prune first ##
## REMOVE SNPs IN LD FOR PCA and use only common SNPs; ie only include SNPs with MAF >= 0.01
plink --noweb --bfile ${MERGE_OUT_PREFIX} --maf 0.01 --indep-pairwise 500 50 0.1 --out ${MERGE_OUT_PREFIX}.common-prune --allow-no-sex
plink --noweb --bfile ${MERGE_OUT_PREFIX} --maf 0.01 --extract ${MERGE_OUT_PREFIX}.common-prune.prune.in --allow-no-sex --make-bed --out ${MERGE_OUT_PREFIX}.common-pruned



##Sieve out CONSERVED REGION SNPS (CRS)
bfile_prefix='exome_gwas_merge.common-pruned'

plink --noweb --bfile ${bfile_prefix} --chr 6 --from-bp 24892021 --to-bp 36892022 --write-snplist --out temp1 --allow-no-sex
plink --noweb --bfile ${bfile_prefix} --chr 11 --from-bp 45043424 --to-bp 57243424 --write-snplist --out temp2 --allow-no-sex
plink --noweb --bfile ${bfile_prefix} --chr 5 --from-bp 43964243 --to-bp 51464243 --write-snplist --out temp3 --allow-no-sex
plink --noweb --bfile ${bfile_prefix} --chr 8 --from-bp 10000 --to-bp 12655629 --write-snplist --out temp4 --allow-no-sex
plink --noweb --bfile ${bfile_prefix} --chr 11 --from-bp 84322352 --to-bp 86322352 --write-snplist --out temp5 --allow-no-sex

cat temp*snplist > snps-conserved-correlated.txt

plink --noweb --bfile ${bfile_prefix} --allow-no-sex --exclude snps-conserved-correlated.txt --make-bed --out ${bfile_prefix}-recode




##convert to ped
plink --bfile ${bfile_prefix}-recode --recode --out ${bfile_prefix}-recode



## obtain the gt file using plink2gt (fratools package)
plink2gt ${bfile_prefix}-recode



## run fcmp, using merged exome-gwas file 
fcmp -c comparelist1.txt ${bfile_prefix}-recode.gt
