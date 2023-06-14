# PD exome-seq Discovery dataset analyses

This repository contains codes and analyses conducted for the Discovery dataset of paper: 


Codes & details for following steps can be found in the respective repository folders: 
- Sample quality checks, `SampleQC`
- Steps prior to variant quality filtering, `prior_VariantQC`
- Variant annotation,`Variant_annot`
- Filtering: Variant Quality and Variant Pathogenecity prediction, `Filtering`
- Gene-based test: Stratified Cochran-Mantel-Haenszel (CMH), `Gene_based_test`


## Summary of pipeline
Steps for coverage, variant calling, joint calling, VQSR, split to different allele types follows GATK recommended pipeline documentation (v3.7) (https://github.com/broadinstitute/gatk-docs/tree/master/gatk3-tooldocs/3.7-0)


![image](https://github.com/fjnlab/PD_exome-seq/assets/58157134/fcfc0263-a844-46c7-a1bf-0ea6fd4fdfa2)


## Contact
For questions, you may contact:

Foo Jia Nee, jianee.foo@ntu.edu.sg

