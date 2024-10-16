#' ---
#' title: "CMH_rare_del"
#' output: html_document
#' ---
#' 

# load libraries
library(purrr)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(psych)
library(data.table)
library(broom)
library(dplyr)
library(tidyr)
library(tidyverse)
library(janitor)

writeLines("\n===== Libraries loaded =====\n")

setname <- "rare_del"    

 
##### 1) Load variant matrix ######
# ---- 1.1) Load annotated variant matrix with per sample per variant genotype ------------------------------------------------------------------------------------------------------------------
writeLines("\n# 1.1) Load annotated variant matrix \n")

tsv <- read_tsv(file="variantannot_samplegeno.tsv", col_names = TRUE)

# prepare keys of ENSG_id to gene symbols and chr no
key_ENSGid_gene <- tsv %>% 
  select(c("GeneID","GeneName")) %>% unique() %>%
  mutate(ENSGid_gene = paste(GeneID,GeneName, sep = '_')) %>%
  relocate(ENSGid_gene)


 
# # ----1.2) Summarise sample-level carrier counts for each gene------------------------------------------------------------------------------------------------------------------
writeLines("\n# 1.2) Summarise sample-level carrier counts for each gene \n")
genocount_tmp <- tsv %>% 
  mutate(ENSGid_gene = paste(GeneID,GeneName, sep = '_')) %>%      # create unique geneid
  select(c("ENSGid_gene","CHH31909_SM:0":"SGH_P_067_SM:1")) %>%    # select unique geneid and all sample genotypes
  mutate(across(c("SAMPLE1:0":"SAMPLE9810:1"), ~ case_when(. == 1 ~ 1, . == 2 ~ 1, TRUE ~ 0))) %>%   # convert genotypes to variant carrier count by recoding heterozygous and homozygous alternative calls, based on genotypes from the first to last sample
  group_by(ENSGid_gene) %>% summarise(across(everything(), sum, na.rm = TRUE), .groups = 'drop') # sum up variant carrier counts per sample per gene

genocount <-  genocount_tmp %>%
   mutate(across(c("SAMPLE1:0":"SAMPLE9810:1"), ~ case_when(. > 0 ~ 1, TRUE ~ 0)))   %>%  # convert variant-level carrier sum to gene-level carrier sum
   t() %>% as.data.frame() %>% row_to_names(row_number = 1) %>% 
   rownames_to_column(var = "LibraryID")
  


# # ----1.3) Annotate gene-level carrier matrix with country of origin and case/control status------------------------------------------------------------------------------------------------------------------
writeLines("\n# 1.3) Annotate gene-level carrier matrix with strata info and case/control status \n")

sample_info <- read.table("sample_country-condition.txt", header = TRUE, sep = "\t")  ## Edit file name of sample demographic file if needed. 

genocount_info <- as.data.table(genocount)[sample_info, on = "LibraryID"] %>%
  select(c("Sample.ID", "Country.condition", contains("ENSG"))) %>%
  mutate(across(starts_with("ENSG"), as.integer))


# splitting into smaller dataframes for processing

genocount_subset_size <- 1000
df_list <- list()

for (i in 1:ceiling(nrow(genocount_info)/genocount_subset_size)) {
  start <- (i-1)*genocount_subset_size + 1
  end <- min(i*genocount_subset_size, nrow(genocount_info))
  df <- genocount_info[start:end, ]
  
  # Print the class of the data frame
  print(class(df))
  print((paste0("df", i," generated")))
  
  # Gather columns and filter by n > 0
  df_tmp <- gather(df, "ENSGid_gene", "n", 3:ncol(df)) %>% 
    filter(n > 0)
  
  # Append the temporary data frame to the list
  df_list[[i]] <- df_tmp
}

genocount_all <- bind_rows(df_list)




#' 
#' # 2) Summarise counts by population
## ----Population and condition summary--------------------------------------------------------------------------------------------------------------
writeLines("\n# 2) Summarise by population \n")

# no of case/control 
HKcase_total <-	70
HKcontrol_total <- 586
KRcase_total <-	1417
KRcontrol_total <- 1040
MALcase_total <- 656
MALcontrol_total <- 114
SGcase_total <- 1955
SGcontrol_total <- 3630
TWcase_total <- 200
TWcontrol_total <- 142

SGMALcase_total <- SGcase_total + MALcase_total
SGMALcontrol_total <- SGcontrol_total + MALcontrol_total

ALLcase_total <- HKcase_total + KRcase_total + MALcase_total + SGcase_total + TWcase_total
ALLcontrol_total <- HKcontrol_total + KRcontrol_total + MALcontrol_total + SGcontrol_total + TWcontrol_total

conditions <- as.matrix(unique(unlist(genocount_info[,2])))
conditions


# # split list to each pop and disease status #
	# 0 -> non-carrier, 1 -> carrier
lvls <- c("0","1")

# "SGMAL case"
gene_SGMALcase <- genocount_all %>% filter(Country.condition %in% c("SG case","MAL case"))
gene_SGMALcase_count_final <- gene_SGMALcase %>% group_by(ENSGid_gene) %>% summarise(SGMAL.case.carrier = sum(n))

# "SGMAL control"
gene_SGMALcont <- genocount_all %>% filter(Country.condition %in% c("SG control","MAL control"))
gene_SGMALcont_count_final <- gene_SGMALcont %>% group_by(ENSGid_gene) %>% summarise(SGMAL.control.carrier = sum(n)) 

# "HK case" 
gene_HKcase <- genocount_all %>% filter(Country.condition ==  "HK case")
gene_HKcase_count_final <- gene_HKcase %>% group_by(ENSGid_gene) %>% summarise(HK.case.carrier = sum(n)) 

# "HK control"
gene_HKcont <- genocount_all %>% filter(Country.condition ==  "HK control")
gene_HKcont_count_final <- gene_HKcont %>% group_by(ENSGid_gene) %>% summarise(HK.control.carrier = sum(n)) 

# "TW case"  
gene_TWcase <- genocount_all %>% filter(Country.condition ==  "TW case")
gene_TWcase_count_final <- gene_TWcase %>% group_by(ENSGid_gene) %>% summarise(TW.case.carrier = sum(n)) 

# "TW control"  
gene_TWcont <- genocount_all %>% filter(Country.condition ==  "TW control")
gene_TWcont_count_final <- gene_TWcont %>% group_by(ENSGid_gene) %>% summarise(TW.control.carrier = sum(n)) 

# "KR case"   
gene_KRcase <- genocount_all %>% filter(Country.condition ==  "KR case")
gene_KRcase_count_final <- gene_KRcase %>% group_by(ENSGid_gene) %>% summarise(KR.case.carrier = sum(n)) 

# "KR control"
gene_KRcont <- genocount_all %>% filter(Country.condition ==  "KR control")
gene_KRcont_count_final <- gene_KRcont %>% group_by(ENSGid_gene) %>% summarise(KR.control.carrier = sum(n)) 


# combine counts from each population and disease status #
gene_all_count_final_carriercount <- list(gene_SGMALcase_count_final,gene_SGMALcont_count_final,
                                          gene_HKcase_count_final,gene_HKcont_count_final,
                                          gene_KRcase_count_final,gene_KRcont_count_final,
                                          gene_TWcase_count_final,gene_TWcont_count_final) %>%
                                     reduce(full_join, by = 'ENSGid_gene') %>% replace(is.na(.), 0)
#write.table(gene_all_count_final_carriercount,file=paste0(setname,"_gene-level_all_carriercount.txt"), quote = F, sep = "\t", row.names= F, col.names = T) 

writeLines("\n==== Summary by population completed ====\n")



#' 
#' # 3) Association tests
writeLines("\n# 3) Association tests \n")

#' # 3.1) Prepare matrix for association test
## ----prepare matrix for association test-----------------------------------------------------------------------------------------------------------
writeLines("\n# # 3.1) Prepare matrix for association test \n")

# calculate non-carriers
genefishercounts <- gene_all_count_final_carriercount %>%
  mutate(
     ALL.case.carrier = SGMAL.case.carrier + HK.case.carrier + TW.case.carrier + KR.case.carrier,
     ALL.control.carrier = SGMAL.control.carrier + HK.control.carrier + TW.control.carrier + KR.control.carrier,
     SGMAL.case.noncarrier = SGMALcase_total - SGMAL.case.carrier,
     SGMAL.control.noncarrier = SGMALcontrol_total - SGMAL.control.carrier,
     HK.case.noncarrier = HKcase_total - HK.case.carrier,
     HK.control.noncarrier = HKcontrol_total - HK.control.carrier,
     TW.case.noncarrier = TWcase_total - TW.case.carrier,
	 TW.control.noncarrier = TWcontrol_total - TW.control.carrier,
     KR.case.noncarrier = KRcase_total - KR.case.carrier,
     KR.control.noncarrier = KRcontrol_total - KR.control.carrier,
     ALL.case.noncarrier = ALLcase_total - ALL.case.carrier,
     ALL.control.noncarrier = ALLcontrol_total - ALL.control.carrier) %>%
	 select(ENSGid_gene,
		   SGMAL.case.carrier, SGMAL.case.noncarrier,
		   SGMAL.control.carrier, SGMAL.control.noncarrier,
		   HK.case.carrier, HK.case.noncarrier,
		   HK.control.carrier, HK.control.noncarrier,
		   KR.case.carrier, KR.case.noncarrier,
		   KR.control.carrier, KR.control.noncarrier,
		   TW.case.carrier, TW.case.noncarrier,
		   TW.control.carrier, TW.control.noncarrier,
		   ALL.case.carrier, ALL.case.noncarrier,
		   ALL.control.carrier, ALL.control.noncarrier)
     
colnames(genefishercounts)
#write.table(genefishercounts, file=paste0(setname,"_genefishercounts.txt"), quote = F, sep = "\t", row.names= T, col.names = T)



#' 
#' 
#' # 3.2) Fisher's test
writeLines("\n# 3.2) Fisher's test \n")
## ----fisher test-----------------------------------------------------------------------------------------------------------

# SGMAL
SGMAL_tmp <- genefishercounts %>% 
  select(c(ENSGid_gene, contains("SGMAL."))) %>% 
  remove_rownames %>% column_to_rownames(var="ENSGid_gene") 

SGMAL_pval <- apply(SGMAL_tmp,1, 
                 function(x) {
                   tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                   fisher.test(tbl, alternative="two.sided")$p.value
                 })

SGMAL_estimate <- apply(SGMAL_tmp, 1, 
                     function(x) {
                       tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                       fisher.test(tbl, alternative="two.sided")$estimate
                     })

SGMAL_ci <- apply(SGMAL_tmp, 1, 
               function(x) {
                 tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                 fisher.test(tbl, alternative="two.sided")$conf.int
               })


# HK
HK_tmp <- genefishercounts %>% 
  select(c(ENSGid_gene, contains("HK."))) %>% 
  remove_rownames %>% column_to_rownames(var="ENSGid_gene") 

HK_pval <- apply(HK_tmp,1, 
                 function(x) {
                   tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                   fisher.test(tbl, alternative="two.sided")$p.value
                 })

HK_estimate <- apply(HK_tmp, 1, 
                     function(x) {
                       tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                       fisher.test(tbl, alternative="two.sided")$estimate
                     })

HK_ci <- apply(HK_tmp, 1, 
               function(x) {
                 tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                 fisher.test(tbl, alternative="two.sided")$conf.int
               })

# KR
KR_tmp <- genefishercounts %>% 
  select(c(ENSGid_gene, contains("KR."))) %>% 
  remove_rownames %>% column_to_rownames(var="ENSGid_gene") 

KR_pval <- apply(KR_tmp,1, 
                 function(x) {
                   tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                   fisher.test(tbl, alternative="two.sided")$p.value
                 })

KR_estimate <- apply(KR_tmp, 1, 
                     function(x) {
                       tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                       fisher.test(tbl, alternative="two.sided")$estimate
                     })

KR_ci <- apply(KR_tmp, 1, 
               function(x) {
                 tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                 fisher.test(tbl, alternative="two.sided")$conf.int
               })

# TW
TW_tmp <- genefishercounts %>% 
  select(c(ENSGid_gene, contains("TW."))) %>% 
  remove_rownames %>% column_to_rownames(var="ENSGid_gene") 

TW_pval <- apply(TW_tmp,1, 
                 function(x) {
                   tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                   fisher.test(tbl, alternative="two.sided")$p.value
                 })

TW_estimate <- apply(TW_tmp, 1, 
                     function(x) {
                       tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                       fisher.test(tbl, alternative="two.sided")$estimate
                     })

TW_ci <- apply(TW_tmp, 1, 
               function(x) {
                 tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                 fisher.test(tbl, alternative="two.sided")$conf.int
               })


genefishercounts$Fisher_SGMAL_pval <- SGMAL_pval
genefishercounts$Fisher_SGMAL_estimate <- SGMAL_estimate
genefishercounts$Fisher_SGMAL_ci1 <- t(SGMAL_ci)[,1]
genefishercounts$Fisher_SGMAL_ci2 <- t(SGMAL_ci)[,2]

genefishercounts$Fisher_HK_pval <- HK_pval
genefishercounts$Fisher_HK_estimate <- HK_estimate
genefishercounts$Fisher_HK_ci1 <- t(HK_ci)[,1]
genefishercounts$Fisher_HK_ci2 <- t(HK_ci)[,2]

genefishercounts$Fisher_KR_pval <- KR_pval
genefishercounts$Fisher_KR_estimate <- KR_estimate
genefishercounts$Fisher_KR_ci1 <- t(KR_ci)[,1]
genefishercounts$Fisher_KR_ci2 <- t(KR_ci)[,2]

genefishercounts$Fisher_TW_pval <- TW_pval
genefishercounts$Fisher_TW_estimate <- TW_estimate
genefishercounts$Fisher_TW_ci1 <- t(TW_ci)[,1]
genefishercounts$Fisher_TW_ci2 <- t(TW_ci)[,2]

head(genefishercounts)
#write.table(genefishercounts, file=paste0(setname,"_gene-based_fisher-ci.txt"), quote=F, sep="\t", row.names = FALSE)



# fixing ENSGid_Gene that has been truncated due to separation with "." and "-" earlier
genefishercounts_edit <- genefishercounts %>% 
                        separate(col= ENSGid_gene, into=c("ENSGid",NA), sep= "_", remove=FALSE) %>% 
			merge(key_ENSGid_gene , by.x="ENSGid", by.y = "GeneID", all = TRUE) %>% 
                        as_tibble() %>% 
                        select(c("ENSGid_gene.y", "SGMAL.case.carrier":"Fisher_TW_ci2")) %>%
			rename(ENSGid_gene = ENSGid_gene.y)
						
genefishercounts_edit$ENSGid_gene<-gsub("'","",genefishercounts_edit$ENSGid_gene)
#write.table(genefishercounts_edit, file=paste0(setname,"_gene-based_fisher-ci_edit.txt"), quote=F, sep="\t", row.names = FALSE)


writeLines("\n==== \"genefishercounts_edit\" RDS generated ====\n")






#' 
#' 
#' # 3.3) CMH
writeLines("\n# 3.3) CMH\n")
## ----make data for gene-based CMH-------------------------------------------------------------------------------------------------------

genefishercounts.tidy <- genefishercounts_edit %>% 
                         select(ENSGid_gene,
							    SGMAL.case.carrier,SGMAL.case.noncarrier, SGMAL.control.carrier, SGMAL.control.noncarrier,
							    HK.case.carrier, HK.case.noncarrier, HK.control.carrier, HK.control.noncarrier,
							    KR.case.carrier, KR.case.noncarrier, KR.control.carrier, KR.control.noncarrier,
							    TW.case.carrier, TW.case.noncarrier, TW.control.carrier, TW.control.noncarrier) %>%
						 gather(-ENSGid_gene, key = country.casecontrol.carrierstatus, value = counts) %>% # gather all case.control into one column
						 separate(country.casecontrol.carrierstatus, into = c("country","case.control","carrier.noncarrier")) %>%  # split into case.control, carrier.noncarrier, and strata columns
						 select(ENSGid_gene,case.control,carrier.noncarrier,country,counts)
genefishercounts.tidy[,5] <- sapply(genefishercounts.tidy[,5], as.numeric)

# create a list of dataframes, one list per gene created
genefishercounts.list <- genefishercounts.tidy %>% split(.$ENSGid_gene)

## create 2x2x3 tables
genefishercounts.tables <- map(genefishercounts.list, ~xtabs(counts ~ case.control + carrier.noncarrier + country, data = .))

# initialise an empty list to store results
no_of_genes <- nrow(genefishercounts_edit) 
gene_results <- vector(mode = "list", length = no_of_genes) # length will be the number of genes

# run CMH for each gene
writeLines("\n === Running CMH ===\n")
for(i in 1:length(genefishercounts.tables)){
  gene_results[[i]] <- mantelhaen.test(genefishercounts.tables[[i]], correct = F)
  names(gene_results) <- genefishercounts_edit$ENSGid_gene
}

## extract results and store in dataframe
gene_results_df <- data.frame(matrix(ncol = 5, nrow = no_of_genes, dimnames = list(NULL, c("CMH_SGMAL_HK_KR_TW_pvalue", "CMH_SGMAL_HK_KR_TW_OR", "CMH_SGMAL_HK_KR_TW_confint1", "CMH_SGMAL_HK_KR_TW_confint2", "CMH_SGMAL_HK_KR_TW_statistic")))) # note that nrow = number of genes

rownames(gene_results_df) <- genefishercounts_edit$ENSGid_gene # name the row names as gene names

# add in CMH_discovery (SGMAL,HK,KR,TW) results
gene_results_CMH_SGMAL_HK_KR_TW_pval <- unlist(lapply(gene_results, function(i) i$p.value)) # pull out p-value
gene_results_df$CMH_SGMAL_HK_KR_TW_pvalue <- gene_results_CMH_SGMAL_HK_KR_TW_pval[1:no_of_genes] # add pvalues to the results dataframe
gene_results_CMH_SGMAL_HK_KR_TW_OR <- unlist(lapply(gene_results, function(i) i$estimate)) # pull out odds ratio
gene_results_df$CMH_SGMAL_HK_KR_TW_OR <- gene_results_CMH_SGMAL_HK_KR_TW_OR[1:no_of_genes] # add odds ratio to results dataframe
gene_results_CMH_SGMAL_HK_KR_TW_conf <- unlist(lapply(gene_results, function(i) i$conf)) %>% split(f = c(1,2))# pull out confidence interval
gene_results_df$CMH_SGMAL_HK_KR_TW_confint1 <- gene_results_CMH_SGMAL_HK_KR_TW_conf[[1]] # add confidence interval to results dataframe
gene_results_df$CMH_SGMAL_HK_KR_TW_confint2 <- gene_results_CMH_SGMAL_HK_KR_TW_conf[[2]] # add confidence interval to results dataframe
gene_results_CMH_SGMAL_HK_KR_TW_X.squared <- unlist(lapply(gene_results, function(i) i$statistic)) # pull out Mantel-Haenszel X-squared values
gene_results_df$CMH_SGMAL_HK_KR_TW_statistic <- gene_results_CMH_SGMAL_HK_KR_TW_X.squared[1:no_of_genes] # add X-squared values to results dataframe

# write gene_results_CMH_SGMAL_HK_KR_TW
#write.table(gene_results_df, file=paste0(setname,"_gene-based_CMH.txt"), quote=F, sep="\t", row.names = T)

writeLines("\n\ === CMH completed === \n")



#' 
#' # 3.4) Find exome-wide significant gene-based CMH test results
## ---------------------------------------------------------------------------------------------------------------------------------------
writeLines("\n# 3.4) Find significant gene-based CMH test results\n")

CMH_SGMAL_HK_KR_TW_sig_CMH_2.5En6 <- gene_results_df %>% filter(CMH_SGMAL_HK_KR_TW_pvalue <= 2.5e-6) %>% arrange(CMH_SGMAL_HK_KR_TW_pvalue)
writeLines("\n OUTPUT: CMH exome-wide significant genes \n")
CMH_SGMAL_HK_KR_TW_sig_CMH_2.5En6

#write.table(CMH_SGMAL_HK_KR_TW_sig_CMH_2.5En6 , file=paste0(setname,"_gene-based_CMH_SGMAL_HK_KR_TW_p2.5e-6.txt"), quote=F, sep="\t", row.names = TRUE)
		#' Exome-wide significance P < 2.5 × 10−6 (0.05/20,000) for gene-based analysis



 
#' # 4) Merge Fisher's and CMH results
## ---------------------------------------------------------------------------------------------------------------------------------------
writeLines("\n# 4) Merge Fisher's and CMH results\n")

genemerge <- merge(genefishercounts_edit, gene_results_df , by.x = "ENSGid_gene", by.y='row.names')
write.table(genemerge, file=paste0(setname,"_gene-based_Fishers-CMH.txt"), quote=F, sep="\t", row.names=FALSE)


writeLines("\n!!! ====== SCRIPT COMPLETED ====== !!!\n")



## !!! END OF SCRIPT !!! ##