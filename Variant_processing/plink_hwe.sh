#!/bin/bash
source ~/.bash_profile

FILE_PREFIX='variants.filtered_DPGQ.uniqID'

source ~/.bash_profile

## convert VCF to PLINK format
plink --vcf ${FILE_PREFIX}.vcf.gz --make-bed --double-id --out ${FILE_PREFIX} --allow-no-sex


## edit phenotype & gender for each sample first, BEFORE running the next PLINK command
## vim ${FILE_PREFIX}.fam

## run HWE
plink --bfile ${FILE_PREFIX} --hardy --out ${FILE_PREFIX}-hwe  ## then extract p-val of UNAFF
