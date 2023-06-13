# Variant Annotation
Code here are based on GRCh37 (hg19) genomic coordinates

## Requirements
- PERL
- [Ensembl VEP](https://asia.ensembl.org/info/docs/tools/vep/script/vep_download.html) (v104) and [cache](https://asia.ensembl.org/info/docs/tools/vep/script/vep_cache.html) database for GRCh37
- [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP) for VEP plugin
- [LOFTEE](https://github.com/konradjk/loftee) for VEP plugin
- [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) to extract Genome Aggregation Database (gnomAD) population-level frequences
- [PolyPhen-2](http://genetics.bwh.harvard.edu/pph2/) tool


## Usage
### 1. Generate unique identifier for VEP submission from `variants.filtered_DPGQ.vcf.gz`
```bash
./run_uniqueID_vcf.sh
```
Input: variants.filtered_DPGQ.vcf.gz

Output: variants.filtered_DPGQ.uniqID.vcf.gz



### 2. Variants annotation was carried out using VEP
```bash
${VEP} -i ${VCF_IN} --cache --dir_cache ${CACHE_DB} --offline --cache_version 104 --use_given_ref \
 --assembly GRCh37 \
 --canonical \
 --numbers \
 --variant_class --hgvs --check_existing \
 --fasta ${REFERENCE} --symbol \
 --plugin dbNSFP,${DBNSFP},${REPLACE_LOGIC},Ensembl_transcriptid,Uniprot_acc_Polyphen2,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
 --dir_plugins ${VEP_plugin} --plugin LoF,loftee_path:${LOFTEE},human_ancestor_fa:${LOFTEE_DB}/human_ancestor.fa.gz,conservation_file:${LOFTEE_DB}/phylocsf_gerp.sql,gerp_bases:${LOFTEE_DB}/GERP_scores.final.sorted.txt.gz,gerp_exons:${LOFTEE_DB}/GERP_scores.exons.txt.gz \
 --output_file ${OUT} –tab –verbose
```

Input: variants.filtered_DPGQ.uniqID.vcf.gz

Output files: variants.DP8GQ20_filtered.uniqID.vep.table.out*



### 3. Filter VEP output & extract relevant columns 
```bash
## to extract required VEP annotation columns
perl VEPannot_extract.pl variants.DP8GQ20_filtered.uniqID.vep.table.out variants.DP8GQ20_filtered.uniqID.vep.table.filtered.out
```

Input: variants.DP8GQ20_filtered.uniqID.vep.table.out

Output: variants.DP8GQ20_filtered.uniqID.vep.table.filtered.out


``` bash
1. Chr : chromosomal position of variant
2. Position : variant genomic coordinates based on GRCh37 (hg19)
3. Ref : reference allele
4. Alt : alternate allele
5. TranscriptID : ENSTID format (Ensembl transcript ID)
6. GeneID : ENSGID format (Ensembl gene ID)
7. GeneName : Gene symbol
8. Uploaded_variation : format takes chr:position:ref:alt; from VEP output
9. Location : variant genomic coordinates based on GRCh37 (hg19); from VEP output
10. Allele : alternate allele; from VEP output
11. Transcript : ENSTID; from VEP output
12. Consequence : predicted effects of allele on transcript, see [Ensembl Calculated consequences](https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html) for the exhaustive list on consequence types
13. Intron_variant : "YES" or "NO" if variant is in intronic region
14. Exon_number : Affected exon numbering, format takes exon number/total; VEP output
15. Intron_number : Affected intron numbering, format takes exon number/total; VEP output
16. Protein_pos : Relative position of amino acid in protein; VEP output
17. Amino_acids : Reference and variant amino acids; VEP output
18. HGVSc : HGVS coding sequence name; VEP output
19. HGVSp : HGVS protein sequence name; VEP output
20. SnpID : Identifier(s) of co-located known variants; VEP output
21. Impact : Subjective impact classification of consequence type; VEP output
22. Strand : 1/-1 for strand; VEP output
23. Variant_class : SO variant class; VEP output
24. Biotype : Biotype of transcript or regulatory feature; VEP output
25. Polyphen2_Transcriptid : Ensembl transcript ids (Multiple entries separated by ";"); dbNSFP VEP plugin output
26. Uniprot_acc_Polyphen2 : Uniprot accession number provided by PolyPhen2 (Multiple entries separated by ";"); dbNSFP VEP plugin output
27. Polyphen2_HDIV_pred : Polyphen2 prediction based on HumDiv Model ("D" probably damaging, "P" possibly damaging, "B" benign); dbNSFP VEP plugin output
28. Polyphen2_HVAR_pred : Polyphen2 prediction based on HumVar Model ("D" probably damaging, "P" possibly damaging, "B" benign); dbNSFP VEP plugin output
29. LoF : Loss-of-function annotation (HC = High Confidence; LC = Low Confidence); LOFTEE VEP plugin output
30. LoF_filter : Reason for LoF not being HC; LOFTEE VEP plugin output
31. LoF_info :  Info used for LoF annotation; LOFTEE VEP plugin output
32. CANONICAL : if annotation is based on Ensembl canonical transcript
```



### 4. Genome Aggregation Database (gnomAD) population-level frequencies (v2.1.1) were added using ANNOVAR
Format input file `variants.filtered_DPGQ.uniqID.vcf` for ANNOVAR according to https://annovar.openbioinformatics.org/en/latest/user-guide/input/ with modifications

```bash
## we have modified output of convert2annovar.pl, to facilitate merging of annotation output later
## convert2annovar.pl is provided as part of ANNOVAR package
perl ${$ANNOVAR_DIR}/convert2annovar.pl -format vcf4 -allsample -withfreq -includeinfo variants.DP8GQ20_filtered.uniqID.vcf > variants.avlist.tmp
cut -f1-8,11 variants.avlist.tmp > variants.avlist

## using ANNOVAR to extract gnomAD exome, gnomAD genome population AF
## to use *_dropped for combine step later
perl ${ANNOVAR_DIR}/annotate_variation.pl -build hg19 -filter -dbtype gnomad211_exome $outdir/variants.avlist -otherinfo $ANNOVAR_DIR/humandb/
perl ${ANNOVAR_DIR}/annotate_variation.pl -build hg19 -filter -dbtype gnomad211_genome $outdir/variants.avlist -otherinfo $ANNOVAR_DIR/humandb/
```



### 5. PolyPhen-2 predictions
PolyPhen-2 HDIV predictions for missense variants were [batch queried](http://genetics.bwh.harvard.edu/pph2/bgi.shtml) with the options of: 
- Classified model: HumDiv
- Genome assembly: GRCh37/hg19
- Transcripts: All
- Annotations: All

Input format example chr1:1267483 G/A

Predictions corresponding to the canonical transcripts (ENST_ID) were extracted based on UniProt_ID (UniProtKB/Swiss-Prot ID and UniProtKB/TrEMBL ID), UCSC_ID (UCSC Stable ID) or RefSeq_ID (RefSeq peptide ID, accession prefix “NP_”) obtained from [Ensembl BioMart for GRCh37](https://grch37.ensembl.org/biomart/martview).

We have also provided the matched IDs of our annotated variants based on canonical transcripts, see `resources/allcanonicalvariants_ID_PolyPhen2.list`



### 6. Combine variant annotated list with variant genotypes
Combine VEP annotation, population level gnomAD frequencies, PolyPhen-2 prediction together with genotypes of pass QC variants to give the matrix (see example files `variantannot_samplegeno.tsv`, `samplegeno.tsv`) with the following headers:

``` bash
## Combine variant annotation outputs
## Output: variantannot.tsv
perl merge_VEP_gnomad.pl variants.DP8GQ20_filtered.uniqID.vep.table.filtered.out variants.avlist.hg19_gnomad211_exome_dropped variants.avlist.hg19_gnomad211_genome_dropped

## Combine annotation.tsv + samplegenotype_sorted_pmiss_hwe
## Output files: variantannot_samplegeno.tsv, samplegeno.tsv
perl merge_annot_geno-pmiss-hwe.pl variantannot.tsv samplegenotype_sorted_pmiss_hwe
```

- Chr
- Position
- Ref
- Alt
- TranscriptID
- GeneID
- GeneName
- gnomAD2.1.1_exomes_AF
- gnomAD2.1.1_exomes_AFR
- gnomAD2.1.1_exomes_AMR
- gnomAD2.1.1_exomes_EAS
- gnomAD2.1.1_exomes_SAS
- gnomAD2.1.1_exomes_NFE
- gnomAD2.1.1_exomes_FIN
- gnomAD2.1.1_genomes_AF
- gnomAD2.1.1_genomes_AFR
- gnomAD2.1.1_genomes_AMR
- gnomAD2.1.1_genomes_EAS
- gnomAD2.1.1_genomes_SAS
- gnomAD2.1.1_genomes_NFE
- gnomAD2.1.1_genomes_FIN
- #Uploaded_variation
- Location
- Allele
- Transcript
- Consequence
- Intron_variant
- Exon_number
- Intron_number
- Protein_pos
- Amino_acids
- HGVSc
- HGVSp
- SnpID
- Impact
- Strand
- Variant_class
- Biotype
- Polyphen2_Transcriptid
- Uniprot_acc_Polyphen2
- Polyphen2_HDIV_pred
- Polyphen2_HVAR_pred
- LoF
- LoF_filter
- LoF_info
- CANONICAL
- Allele_Freq_Control
- Allele_Freq_Case
- Allele_Freq_All_Samples
- Missing_Call_Rate_Control
- Missing_Call_Rate_Case
- Missing_Call_Rate_All
- Filter
- Control_NA
- Control_0
- Control_1
- Control_2
- Case_NA
- Case_0
- Case_1
- Case_2
- Pmiss
- Pmiss_OddRatio
- Pmiss_Log_OR
- Pmiss_abs_logOR
- total
- HWE_ctrl_pval
- Chr
- Position
- Ref
- Alt
- TranscriptID
- Per-sample genotype (x N columns)
- Allele_Freq_Control
- Allele_Freq_Case
- Allele_Freq_All_Samples
- Missing_Call_Rate_Control
- Missing_Call_Rate_Case
- Missing_Call_Rate_All
- Filter
- Control_NA
- Control_0
- Control_1
- Control_2
- Case_NA
- Case_0
- Case_1
- Case_2
- Pmiss
- Pmiss_OddRatio
- Pmiss_Log_OR
- Pmiss_abs_logOR
- total
- HWE_ctrl_pval
