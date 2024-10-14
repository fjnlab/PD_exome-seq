#!/bin/bash
## run pca using fpca from fratools
########################################

source ~/.bash_profile

## prepare ibd-pruned PLINK file for fpca run
PREFIX='exome_samples_removed.common-pruned-recode-qced'

##convert bed file to map & ped file
plink --noweb --bfile ${PREFIX} --recode --out ${PREFIX}

##Convert ped file to GT and transpose
plink2gt ${PREFIX}
ftranspose ${PREFIX}.gt

##Remove samples failing QC step 1-5 before running fpca
fsieve -s samplesfail_QC1to5.txt -c -m ${PREFIX}.tg
fpca -i sieved-${PREFIX}.tg


