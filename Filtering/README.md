# Filtering: Variant Quality and Variant Pathogenicity predictions

***Variant quality filtering***

Variant quality score recalibration (VQSR, GATK v3.7) was performed to exclude low quality SNP calls. Low quality genotype calls (DP <8 and GQ <20) were excluded. Variants with low call rates in cases, controls and all samples respectively (<95%), failed tests for Hardy-Weinberg equilibrium (P<1×10-6) and high differential missingness between cases and controls (P<1×10-4) were excluded from further analyses.

***Variant pathogenicity filtering***

Functional predictions for missense variants were derived from PolyPhen-2 (v2.2.3r406) while the VEP plugin LOFTEE (v1.0.3) was used for functional prediction of stop gain or stop loss, splice site disrupting and frameshift variants. Each variant was classified as damaging based on PolyPhen2 (probably or possibly damaging) and LOFTEE (high confidence) predictions.

Genome Aggregation Database (gnomAD) population-level frequencies (v2.1.1) were added using ANNOVAR and used for classifying rare variants (MAF ≤1%).


# Requirements
- PERL

# Usage
``` bash
## filtering_criteria_parameters.txt: a tab-delimited text file, see details for parameters of filtering criteria
## output directory with backslash “/” included, ie “output/” without quotes
perl filter_annotation.pl variantannot_samplegeno.tsv samplegeno.tsv <filtering_criteria_parameters.txt> <output directory/>
```

Input files: variantannot_samplegeno.tsv, samplegeno.tsv, filtering_criteria_parameters.txt

Filtered lists in output directory: variantannot_samplegeno_final_all_results_filtered.tsv, samplegeno_filtered.tsv

## (1) Filter for rare variants

``` bash
mkdir rare_qced
perl filter_annotation.pl variantannot_samplegeno.tsv samplegeno.tsv rare_qced_params.txt rare_qced/
```

Input files: variantannot_samplegeno.tsv, samplegeno.tsv, rare_qced_params.txt

Output files: rare_qced/variantannot_samplegeno_final_all_results_filtered.tsv, rare_qced/samplegeno_filtered.tsv


## (2) Filter for rare pathogenic variants

``` bash
mkdir rare_qced_del
perl filter_annotation.pl variantannot_samplegeno.tsv samplegeno.tsv rare_qced_del_params.txt rare_qced_del/
```

Input files: variantannot_samplegeno.tsv, samplegeno.tsv, rare_qced_del_params.txt

Output files: rare_qced_del/variantannot_samplegeno_final_all_results_filtered.tsv, rare_qced_del/samplegeno_filtered.tsv


## (3) Filter for rare benign variants

``` bash
mkdir rare_qced_benign
perl extract_benign.pl rare_qced/variantannot_samplegeno_final_all_results_filtered.tsv rare_qced_del/variantannot_samplegeno_final_all_results_filtered.tsv rare_qced_benign/

```

Input files: rare_qced/variantannot_samplegeno_final_all_results_filtered.tsv, rare_qced_del/variantannot_samplegeno_final_all_results_filtered.tsv

Output files: rare_qced_benign/variantannot_samplegeno_final_all_results_filtered.tsv, rare_qced_benign/samplegeno_filtered.tsv


