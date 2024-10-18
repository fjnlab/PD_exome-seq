# PD exome-seq Discovery dataset analyses

This repository contains codes and analyses conducted for the Discovery dataset.

Codes & details for following steps can be found in the respective repository folders: 
1. Sample quality checks, `SampleQC`
2. Steps prior to variant quality filtering, `Variant_processing`
3. Variant annotation,`Variant_annot`
4. Filtering: Variant Quality and Variant Pathogenicity prediction, `Filtering`
5. Gene-based test: Stratified Cochran-Mantel-Haenszel (CMH), `Gene_based_test`


## Summary of pipeline
Steps for coverage, variant calling, joint calling, VQSR, split to different allele types follows GATK recommended pipeline documentation (v3.7) (https://github.com/broadinstitute/gatk-docs/tree/master/gatk3-tooldocs/3.7-0).

Our pipeline is based on GRCh37 (hg19) genomic coordinates.


![image](https://github.com/user-attachments/assets/3c0da211-004d-498b-a94d-75f715a56c1f)




## Contact
For questions, you may contact:

Jia Nee Foo, jianee.foo@ntu.edu.sg
Elaine Chew, elaine.chew@ntu.edu.sg
Michelle Lian, michellelian@ntu.edu.sg

