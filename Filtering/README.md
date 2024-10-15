# Filtering: Variant Quality and Variant Pathogenicity predictions

***Variant quality filtering***

Variant quality score recalibration (GATK v3.7 VQSR) was performed to exclude low quality SNP calls. Low quality genotype calls (DP <8 and GQ <20) were excluded. Variants with low call rates in cases, controls and all samples respectively (<95%), failed tests for Hardy-Weinberg equilibrium (P<1×10-6) and high differential missingness between cases and controls (P<1×10-4) were excluded from further analyses.

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

### (a) Rare pathogenic variants based on VEP canonical transcripts were filtered with the following criteria (criteria list is tab-delimited, see example file `rare_qced_params.txt`)
```
CONSEQUENCE	missense_variant	or	stop_lost	or	start_lost	or	coding_sequence_variant	or	protein_altering_variant	or	stop_gained	or	frameshift_variant	or	inframe_insertion	or	inframe_deletion	or	splice_region_variant	or	splice_donor_5th_base_variant	or	splice_donor_region_variant	or	splice_polypyrimidine_tract_variant
INTRON_VARIANT	NO
CANONICAL	YES
GNOMAD_EXOMES_EAS	NA	or	<=0.01
GNOMAD_GENOMES_EAS	NA	or	<=0.01
FILTER	PASS
MCR_CASE	<0.05
MCR_CTRL	<0.05
MCR_ALL	<=0.05
P_MISS	>=0.001
HWE_CTRL_PVAL	>=0.000001
```

### (b) Filter for deleterious variants (see example file `rare_qced_del_params.txt`)
Stating criteria in following manner (separated by "#####") will enable script to filter each set of criteria independently, then join (& take the unique of) each filtered set into 1 combined filtered dataset

- Filter for damaging or probably damaging based on PolyPhen-2 -> set 1
- Filter for highly confidence based on LOFTEE -> set 2
- Join both set 1 & 2, and obtain the unique variant list to form rare pathogenic & deleterious filtered set
```
POLYPHEN2_HDIV_pred	D	or	P
#####
LOF	HC
```

### (c) Filter for benign variants (Rare pathogenic set ***minus*** rare deleterious set)
``` bash
## Rare pathogenic variants – rare deleterious
## output directory with backslash “/” included, ie “output/” without quotes
mkdir <rare_qced_benign>
perl extract_benign.pl <rare_qced dir>/variantannot_samplegeno_final_all_results_filtered.tsv <rare_qced_del dir>l/variantannot_samplegeno_final_all_results_filtered.tsv <rare_qced_benign/>

```

Input files: variantannot_samplegeno_final_all_results_filtered.tsv in rare_qced, rare_qced_del directories respectively

Filtered lists in output directory: variantannot_samplegeno_final_all_results_filtered.tsv, samplegeno_filtered.tsv



### (d) Filter for CDS variants (to obtain preQC biallelic snps to conduct sampleQC, see example file `preQC_CDS_params.txt`)
```
INTRON_VARIANT	NO
CANONICAL	YES
FILTER	PASS
```

