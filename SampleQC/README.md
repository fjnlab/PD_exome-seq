# Sample Quality Checks
We removed patients with sex discordance between genetically inferred sex obtained from exome sequencing with patient records. We also removed samples showing poor genotyping concordance (<95% concordance) between exome sequence data and previous fingerprinting with genome-wide association arrays. Samples with excessive heterozygosity (defined as >3.5 standard deviations from the mean), unusually high exome-wide singleton counts, and low genotyping completion rate (<95%) were also excluded. We used identity-by-descent analysis to identify related sample pairs up to the third degree (IBD >0.125). For these pairs, the sample with the lower call rate was removed. Lastly, samples with outlying genetic ancestry in ancestry principal component analysis were also removed. We took care to treat variant identification in all samples equally regardless of affection status. 

Code here was conducted at the preQC stage and based on GRCh37 (hg19) genomic coordinates.

# Requirements
- [PLINK v1.90](https://www.cog-genomics.org/plink/)
- PERL
- PERL scripts, [fratools](https://github.com/atks/fratools) repository
- R or Rstudio, and RColorBrewer library

# Usage
Input: Pre-QC biallelic snps dataset filtered for DP, GQ, PASS VQSR, CDS region, canonical

| Step | Scripts |
| -----| ------- |
| ***Sample quality control step 1: Gender discrepancies*** </br> Gender on patient records discordant with sex genetically inferred from exome sequencing data | plink_gendercheck_chrXY_Xcount.sh or plink_gendercheck_chrXY_Ycount.sh[^1] |
| ***Sample quality control step 2: Genotype discordance*** </br> Sample exome genotype discordant (<95%) with genome-wide association assay (GWAS) data | genotype_concordance.sh |
| ***Sample quality control step 3*** </br> Removed samples with: </br> (a) high singleton counts (>500) and/or </br> (b) high % het snps (>3.5 standard deviations from mean) | sample_summarize.pl <br/><br/> Usage: perl sample_summarize.pl <preQC_CDS dir/variantannot_samplegeno_final_all_results_filtered.tsv> <preQC_CDS dir/samplegeno_filtered.tsv> <output_dir/> <br/><br/> - include "/" backslash in output dir name| 
| ***Sample quality control step 4*** </br> Samples with high genotype missing calls (>95%) | prior_ibd.sh |
| ***Sample quality control step 5*** </br> Removed samples that: </br> (a) high genotype missing calls in steps prior to [Identity-By-Descent (IBD)](https://zzz.bwh.harvard.edu/plink/ibdibs.shtml) test (MAF >1%, excluded SNPs in linkage disequilibrium, genotype missing calls >95%) </br> (b) fail IBD test (IBD >0.125) | ibd.sh |
| ***Sample quality control step 6*** </br> Removed samples that are PC1-PC4 outliers in ancestry principal component analysis | fpca.sh, fpca-plots.r |


[^1]: PLINK users Google-group topic ["On gender calling"](https://groups.google.com/forum/#!topic/plink2-users/28LESfNj64A)





