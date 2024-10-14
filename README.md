# PD exome-seq Discovery dataset analyses

This repository contains codes and analyses conducted for the Discovery dataset of paper: 


Codes & details for following steps can be found in the respective repository folders: 
- Sample quality checks, `SampleQC`
- Steps prior to variant quality filtering, `prior_VariantQC`
- Variant annotation,`Variant_annot`
- Filtering: Variant Quality and Variant Pathogenicity prediction, `Filtering`
- Gene-based test: Stratified Cochran-Mantel-Haenszel (CMH), `Gene_based_test`


## Summary of pipeline
Steps for coverage, variant calling, joint calling, VQSR, split to different allele types follows GATK recommended pipeline documentation (v3.7) (https://github.com/broadinstitute/gatk-docs/tree/master/gatk3-tooldocs/3.7-0).

Our pipeline is based on GRCh37 (hg19) genomic coordinates.


![image](https://github.com/user-attachments/assets/3c0da211-004d-498b-a94d-75f715a56c1f)




## Contact
For questions, you may contact:

Foo Jia Nee, jianee.foo@ntu.edu.sg

