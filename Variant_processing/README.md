# Variant processing
Variant quality control procedures were performed using VCFtools (version 0.1.16), BCFtools (version 1.9) and PLINK (version 1.9). Low quality genotype calls (DP <8 and GQ <20) were excluded. At the filtering step (conducted later), variants with low call rates in cases, controls and all samples respectively (<95%), failed tests for Hardy-Weinberg equilibrium (P<1×10-6) and high differential missingness between cases and controls (P<1×10-4) were excluded from further analyses.

Code here were conducted before variant quality filtering and based on GRCh37 (hg19) genomic coordinates. 
At the filtering stage, variant QC was conducted together with variant pathogenicity prediction filtering. 


# Requirements
- bgzip (as part of [BCFtools](http://www.htslib.org/download/), HTSlib utils)
- tabix (as part of [BCFtools](http://www.htslib.org/download/), HTSlib utils)
- [VCFtools](https://vcftools.github.io/index.html)
- [PLINK v1.90](https://www.cog-genomics.org/plink/)
- PERL
- PERL module [Text::NSP::Measures::2D::Fisher::twotailed](https://metacpan.org/pod/Text::NSP::Measures::2D::Fisher::twotailed)


# Usage
## 1. Generate VCF filtered file based on read depth per sample site (DP) and genotype quality (GQ)
``` bash
./filter_DP8GQ20_vcf.sh
```

Input file: variants.vcf.gz (see [GATK VCF format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format)) \
Output file: variants.filtered_DPGQ.vcf.gz


## 2. Generate unique identifier (also for VEP annotation later) from `variants.filtered_DPGQ.vcf`
``` bash
./run_uniqueID_vcf.sh
```

Input file: variants.filtered_DPGQ.vcf.gz \
Output file: variants.filtered_DPGQ.uniqID.vcf.gz


## 3. Calculate Hardy-Weinberg Equilibrium (HWE) using PLINK (v1.90)
See also https://zzz.bwh.harvard.edu/plink/summary.shtml#hardy

``` bash
plink --vcf ${FILTERED_VCF_PREFIX}.vcf.gz --make-bed --double-id --out ${FILTERED_VCF_PREFIX} --allow-no-sex

##edit fam file for each sample’s phenotype and gender first, BEFORE running the next PLINK command
vim ${FILTERED_VCF_PREFIX}.fam  

plink --bfile ${FILTERED_VCF_PREFIX} --hardy --out ${FILTERED_VCF_PREFIX}-hwe  ## then extract p-val of UNAFF 

```
Input files: variants.filtered_DPGQ PLINK formatted files (*.bed, *.bim, *.fam) \
Output file: <$FILTERED_VCF_PREFIX>-hwe.hwe


## 4. Calculate quality control metrics (missing call rates) and genotype count summary
Genotype codes:
- “NA” – no calls or fail DP, GQ filter
- “0” – homozygous for reference allele
- “1” – heterozygous
- “2” – homozygous for alternate allele

``` bash
perl pheno_genotyping_DPGQ.pl <VCF> <samplelist.txt> <output dir/>
cat samplegenotype_unsorted | awk 'NR==1; NR > 1 {print $0 | "sort -V -k1,1 -k2,2n"}' > samplegenotype_sorted   ##sort by genomic position

<samplelist.txt> - A tab-delimited sample list (without header) with format: SampleName Phenotype (as "case" or "control")
```

Input files: variants.filtered_DPGQ.vcf, samplelist.txt \
Output file: samplegenotype_sorted


## 5. Calculate P(miss) for each variant
- Calculated columns (“Pmiss”, “Pmiss_OddRatio”, Pmiss_Log_OR”, Pmiss_abs_logOR”, “total”) added to SampleGenotype_sorted
- NOTE: ensure last 8 columns of input file `samplegenotype_sorted` are the genotype count summary (“Control_NA”, “Control_0”, “Control_1”, “Control_2”, “Case_NA”, “Case_0”, “Case_1”, “Case_2”)

``` bash
perl calculate_pmiss.pl samplegenotype_sorted samplegenotype_sorted_pmiss
``` 

Input file: samplegenotype_sorted \
Output file: samplegenotype_sorted_pmiss


## 6. Combine all variant metrics: Merge `samplegenotype_sorted_pmiss` and `*.hwe`
``` bash
perl merge_geno_pmiss_hwe.pl ${FILTERED_VCF_PREFIX}-hwe.hwe samplegenotype_sorted_pmiss samplegenotype_sorted_pmiss_hwe
```

Input: ${FILTERED_VCF_PREFIX}-hwe.hwe, samplegenotype_sorted_pmiss \
Output file: samplegenotype_sorted_pmiss_hwe (see example files `samplegenotype_sorted_pmiss_hwe`, `samplegenotype_sorted_pmiss_hwe.descriptor`)


```
1. Chr : chromosomal position of variant
2. Position : variant genomic coordinates based on GRCh37 (hg19)
3. Ref : reference allele
4. Alt : alternate allele
5 - 9814. Per-sample genotype (x N sample columns, i.e. 9810 total samples in Discovery cohort): for genotype codes see previous step #4 
9815. Allele_Freq_Control : minor allele frequencies in control samples
9816. Allele_Freq_Case : minor allele frequencies in case samples
9817. Allele_Freq_All_Samples : minor allele frequencies in all samples
9818. Missing_Call_Rate_Control : calculated from “NA” genotype calls in control samples; for variant quality QC
9819. Missing_Call_Rate_Case : calculated from “NA” genotype calls in case samples; for variant quality QC
9820. Missing_Call_Rate_All : calculated from “NA” genotype calls in all samples; for variant quality QC
9821. Filter : VQSR filter status; for variant quality QC
9822. Control_NA : no. of “NA” (no calls or fail variant quality QC) genotype calls in control samples
9823. Control_0 : no. of “0” (homozygous for reference allele) genotype calls in control samples
9824. Control_1 : no. of “1” (heterozygous) genotype calls in control samples
9825: Control_2 : no. of “2” (homozygous for alternate allele) genotype calls in control samples
9826: Case_NA : no. of “NA” (no calls or fail variant quality QC) genotype calls in case samples
9827. Case_0 : no. of “0” (homozygous for reference allele) genotype calls in case samples
9828: Case_1 : no. of “1” (heterozygous) genotype calls in case samples
9829: Case_2 : no. of “2” (homozygous for alternate allele) genotype calls in case samples
9830: Pmiss : p-value for differential missingness between cases and controls; for variant quality QC
9831: Pmiss_OddRatio : odds ratio for differential missingness between cases and controls; for variant quality QC
9832: Pmiss_Log_OR : (log base 10) of odds ratio for differential missingness between cases and controls; for variant quality QC
9833: Pmiss_abs_logOR : absolute value of (log base 10) odds ratio for differential missingness between cases and controls; for variant quality QC
9834: total : total as calculated in 2 x 2 Chi-squared for pmiss calculation; for variant quality QC
9835: HWE_ctrl_pval : p-value of Hardy-Weinberg equilibrium test for controls; for variant quality QC
```





