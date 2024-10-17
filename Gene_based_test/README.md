# Stratified Cochran-Mantel-Haenszel gene-based test
For the discovery cohort, a stratified Cochran-Mantel-Haenszel (CMH) test was used to evaluate gene-based burden across the exomes of the study participants from the 5 countries studied (Singapore, Malaysia, Hong Kong, South Korea, Taiwan). Singapore and Malaysia samples were considered as one stratum due to the similarity in genetic background. Exome-wide significance was preset at P < 2.5 × 10−6 (two-tailed), taking into account multiple hypothesis testing correction for an estimated 20,000 protein-coding genes in the human genome. 

# Requirements
- Conda or Miniconda3 installed (https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh)

# Installation
Import conda environment and load required packages

```bash
conda env create -n essentials_env --file essentials_env.yml
```

# Input files required
## 1)	filtered annotated variant list `variantannot_samplegeno_final_all_results_filtered.tsv`
- based on filtered annotation output of `filter_annotation.pl` script
- contains matrix of rare deleterious variants with variant annotations (VEP annotations, pathogenicity predictions, etc..) and per sample genotype calls.


## 2) sample information [`sample_country-condition.txt`](/Gene_based_test/sample_country-condition.txt)
- a tab-delimited text file of Sample.ID, LibraryID, strata, status ("case" or "control)
  
  
# Usage
First, activate environment `essentials_env`

```bash
conda activate essentials_env
```

To run R script for stratified CMH test
```bash
Rscript --verbose --no-save CMH_rare-del_gene-based.R >& stderr-out_rare-del_gene-based.txt
```

Output files: 
- rare_del_gene-based_Fishers-CMH.txt: Gene-based test results (see [`rare_del_gene-based_Fishers-CMH.txt.descriptor`](/Gene_based_test/rare_del_gene-based_Fishers-CMH.txt.descriptor))
- stderr-out_rare-del_gene-based.txt: Script running errors, warning and message (see example output)

