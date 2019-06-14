library(tidyverse)

#This is for analyzing the association between environmental variables and SNPs/INVs

directory <- "/media/owens/Copper/wild_gwas_2018/env_associations/for_Greg/"

snp_columns <- read_tsv(paste(directory,"snps_HA412/annuus/column.names.annuus.snps.varout.txt",sep=""))
colnames(snp_columns)[-30]
snps <- read_tsv(paste(directory,"snps_HA412/annuus/var_out_annuus_2018_HA412_maf_0.03_all_covariates_spearman.txt",sep=""),
                 col_names = colnames(snp_columns)[-30])

filtered_sites <- read_tsv("/media/owens/Copper/wild_gwas_2018/gwas/annuus/Annuus.tranche90.snp.gwas.90.bi.remappedHa412HO.jan09noinv.ldfilter.sites",
                           col_names = "snp_id")

snps %>%
  inner_join(.,filtered_sites) -> snps

snps %>%
  gather(.,env_variable,pvalue, -snp_id,-chrom,-pos) %>%
  mutate(type = "snp") -> snps_long
variables <- unique(snps_long$env_variable)

inv_columns <- read_tsv(paste(directory,"inversions_HA412/annuus/column.names.varout.annuus.txt",sep=""))
invs <- read_tsv(paste(directory,"inversions_HA412/annuus/var_out_annuus_inversions_HA412_vcf_all_covariates_spearma_pvalues.txt",sep=""),
                col_names = colnames(inv_columns))

invs %>%
  gather(.,env_variable,pvalue, -inversion_pos) %>%
  mutate(type = "inv",chr="NA",pos=1) %>%
  rename(snp_id = inversion_pos) -> invs_long
  
  
ecdf(snps_long %>% filter(env_variable == "CMD_e") %>% 
       mutate(logp = abs(log10(pvalue))) %>%
       pull(logp))
  
snps_long %>%
  filter(env_variable == "DD_0_e") %>%
  ggplot(.,aes(round(abs(log10(pvalue)),2))) + stat_ecdf(geom = "step") +
  geom_vline(xintercept = abs(log10(invs_long$pvalue[which(invs_long$env_variable == "DD_0_e")])),
             linetype="dashed") +
  geom_hline(yintercept=0.95,linetype="dotted") +
  theme_bw()

perm_pvalues <- tibble(env_variable=character(),pvalue=numeric())

for (i in 1:length(variables)){
  variable <- variables[i]
  
  
  snp_tmp <- snps_long %>%filter(env_variable == variable)
  inv_tmp <- invs_long %>%filter(env_variable == variable)
  
  

  
  #permutation test 
  n_perm <- 5000
  min_inv_pvalue <- min(inv_tmp$pvalue)
  min_pvalue <- vector() 
  for (x in 1:n_perm){
    min_pvalue[x] <- sample_n(snp_tmp, nrow(inv_tmp)) %>% pull(pvalue) %>% min()
  }
  perm_pvalue <- sum(min_pvalue < min_inv_pvalue)/(n_perm+1)
  
  tmp <- tibble(env_variable=variable,pvalue=perm_pvalue)
  perm_pvalues <- rbind(perm_pvalues, tmp)
}


as.matrix(invs[,2:27])

invs.matrix <- as.matrix(invs[,2:27])

rownames(invs.matrix) <- colnames(invs.matrix)
invs.dendro <- as.dendrogram(hclust(d = dist(x = invs.matrix)))
invs.order <- order.dendrogram(invs.dendro)
colnames(invs.matrix)[invs.order]
# Order the levels according to their position in the cluster
invs_long$env_variable <- factor(x = invs_long$env_variable,
                               levels = invs_long$env_variable[invs.order], 
                               ordered = TRUE)



perm_pvalues %>%
  ggplot(.,aes(pvalue)) + geom_histogram()

invs_long %>%
  ggplot(.,aes(pvalue)) + geom_histogram() +
  theme_bw()
invs_long %>%
  ggplot(.,aes(x=snp_id,y=env_variable,fill=abs(log10(pvalue)))) + geom_tile() +
  scale_fill_viridis_c(option="magma") +
  theme(axis.text.x=element_text(angle=60, hjust=1))


