#!/bin/bash
################################################################
## create unique variant identifier of format: chr:pos:ref:alt
################################################################

source ~/.bash_profile

## manually assign unique id
input_vcf='variants.filtered_DPGQ.vcf'
output='variants.filtered_DPGQ.uniqID.vcf'

unpigz -p 10 ${input_vcf}.gz
perl input_uniqID_vcf.pl ${input_vcf} ${output}

bgzip --threads 10 -c ${input_vcf} > ${input_vcf}.gz
tabix -p vcf ${input_vcf}.gz

bgzip --threads 10 -c ${output} > ${output}.gz
tabix -p vcf ${output}.gz
