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

<br/>

## 2. Generate unique identifier (also for VEP annotation later) from `variants.filtered_DPGQ.vcf`
``` bash
./run_uniqueID_vcf.sh
```

Input file: variants.filtered_DPGQ.vcf.gz \
Output file: variants.filtered_DPGQ.uniqID.vcf.gz

<br/>

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

<br/>

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

<br/>

## 5. Calculate P(miss) for each variant
- Calculated columns (“Pmiss”, “Pmiss_OddRatio”, Pmiss_Log_OR”, Pmiss_abs_logOR”, “total”) added to SampleGenotype_sorted
- NOTE: ensure last 8 columns of input file `samplegenotype_sorted` are the genotype count summary (“Control_NA”, “Control_0”, “Control_1”, “Control_2”, “Case_NA”, “Case_0”, “Case_1”, “Case_2”)

``` bash
perl calculate_pmiss.pl samplegenotype_sorted samplegenotype_sorted_pmiss
``` 

Input file: samplegenotype_sorted \
Output file: samplegenotype_sorted_pmiss

<br/>

## 6. Combine all variant metrics: Merge `samplegenotype_sorted_pmiss` and `*.hwe`
``` bash
perl merge_geno_pmiss_hwe.pl ${FILTERED_VCF_PREFIX}-hwe.hwe samplegenotype_sorted_pmiss samplegenotype_sorted_pmiss_hwe
```

Input: ${FILTERED_VCF_PREFIX}-hwe.hwe, samplegenotype_sorted_pmiss \
Output file: samplegenotype_sorted_pmiss_hwe (see `samplegenotype_sorted_pmiss_hwe.descriptor`)






