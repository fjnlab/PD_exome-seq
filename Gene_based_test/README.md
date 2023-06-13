# Stratified Cochran-Mantel-Haenszel gene-based test
For the discovery cohort, a stratified Cochran-Mantel-Haenszel (CMH) test was used to evaluate gene-based burden across the exomes of the study participants from the 5 countries studied (Singapore, Malaysia, Hong Kong, South Korea, Taiwan). Singapore and Malaysia samples were considered as one stratum due to the similarity in genetic background. Exome-wide significance was preset at P < 2.5 × 10−6 (two-tailed), taking into account multiple hypothesis testing correction for an estimated 20,000 protein-coding genes in the human genome. 

## Requirements
- Conda or Miniconda3 installed (https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh)

## Installation
Import conda environment and load required packages

```bash
conda env create -n essentials_env --file essentials_env.yml
```

## Input files required
### 1)	annotated variant list `rare_del_variantlist.tsv`
- based on filtered annotation output of `filter_annotation.pl` script
- contains matrix of all rare deleterious variants with variant annotations (VEP annotations, pathogenicity predictions, etc..) and per sample genotype calls.

### 2) sample information `sample_country-condition.txt`
- a tab-delimited text file of Sample.ID, Library.ID, strata, case or control status
  
  
## Usage
First, activate environment `essentials_env`

```bash
conda activate essentials_env
```

To run R script for stratified CMH test
```bash
Rscript --verbose --no-save CMH_rare-del_gene-based.R >& stderr-out_rare-del_gene-based.txt
```

Output files: 
- rare_del_gene-based_Fishers-CMH.txt: Gene-based test results (see example output)
- stderr-out_rare-del_gene-based.txt: Script running errors, warning and message (see example output)


Column description of gene-based test results output file `rare_del_gene-based_Fishers-CMH.txt`:
```
1. ENSGid_gene: ENSG_ID and Gene Symbol
2. SGMAL.case.carrier: Number of Singapore and Malaysia cases that carry rare deleterious variants in gene
3. SGMAL.case.noncarrier: Number of Singapore and Malaysia cases that do not carry rare deleterious variants in gene
4. SGMAL.control.carrier: Number of Singapore and Malaysia controls that carry rare deleterious variants in gene
5. SGMAL.control.noncarrier: Number of Singapore and Malaysia controls that do not carry rare deleterious variants in gene
6. HK.case.carrier: Number of Hong Kong cases that carry rare deleterious variants in gene
7. HK.case.noncarrier: Number of Hong Kong cases that do not carry rare deleterious variants in gene
8. HK.control.carrier: Number of Hong Kong controls that carry rare deleterious variants in gene
9. HK.control.noncarrier: Number of Hong Kong controls that do not carry rare deleterious variants in gene
10. KR.case.carrier: Number of South Korea cases that carry rare deleterious variants in gene
11. KR.case.noncarrier: Number of South Korea cases that do not carry rare deleterious variants in gene
12. KR.control.carrier: Number of South Korea controls that carry rare deleterious variants in gene
13. KR.control.noncarrier: Number of South Korea controls that do not carry rare deleterious variants in gene
14. TW.case.carrier: Number of Taiwan cases that carry rare deleterious variants in gene
15. TW.case.noncarrier: Number of Taiwan cases that do not carry rare deleterious variants in gene
16. TW.control.carrier: Number of Taiwan controls that carry rare deleterious variants in gene
17. TW.control.noncarrier: Number of Taiwan controls that do not carry rare deleterious variants in gene
18. ALL.case.carrier: Total number of Discovery cases that carry rare deleterious variants in gene
19. ALL.case.noncarrier: Total number of Discovery cases that do not carry rare deleterious variants in gene
20. ALL.control.carrier: Total number of Discovery controls that carry rare deleterious variants in gene
21. ALL.control.noncarrier: Total number of Discovery cases that do not carry rare deleterious variants in gene
22. Fisher_SGMAL_pval: Fisher’s test P for Singapore and Malaysia
23. Fisher_SGMAL_estimate: Fisher’s test odds ratio for Singapore and Malaysia 
24. Fisher_SGMAL_ci1: Fisher’s test lower 95% confidence interval for Singapore and Malaysia
25. Fisher_SGMAL_ci2: Fisher’s test upper 95% confidence interval for Singapore and Malaysia
26. Fisher_HK_pval: Fisher’s test P for Hong Kong
27. Fisher_HK_estimate: Fisher’s test odds ratio for Hong Kong
28. Fisher_HK_ci1: Fisher’s test lower 95% confidence interval for Hong Kong
29. Fisher_HK_ci2: Fisher’s test upper 95% confidence interval for Hong Kong
30. Fisher_KR_pval: Fisher’s test P for South Korea
31. Fisher_KR_estimate: Fisher’s test odds ratio for South Korea
32. Fisher_KR_ci1: Fisher’s test lower 95% confidence interval for South Korea
33. Fisher_KR_ci2: Fisher’s test upper 95% confidence interval for South Korea 
34. Fisher_TW_pval: Fisher’s test P for Taiwan 
35. Fisher_TW_estimate: Fisher’s test odds ratio for Taiwan 
36. Fisher_TW_ci1: Fisher’s test lower 95% confidence interval for Taiwan 
37. Fisher_TW_ci2: Fisher’s test upper 95% confidence interval for Taiwan 
38. CMH_SGMAL_HK_KR_TW_pvalue: Stratified Cochran-Mantel-Haenszel gene-based test P for Discovery cohort
39. CMH_SGMAL_HK_KR_TW_OR: Stratified Cochran-Mantel-Haenszel gene-based test odds ratio for Discovery cohort
40. CMH_SGMAL_HK_KR_TW_confint: Stratified Cochran-Mantel-Haenszel gene-based test lower 95% confidence interval for Discovery cohort
41. CMH_SGMAL_HK_KR_TW_confint2: Stratified Cochran-Mantel-Haenszel gene-based test upper 95% confidence interval for Discovery cohort
42. CMH_SGMAL_HK_KR_TW_statistic: Stratified Cochran-Mantel-Haenszel gene-based test statistics for Discovery cohort
```
