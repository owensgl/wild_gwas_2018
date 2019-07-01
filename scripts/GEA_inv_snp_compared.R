library(tidyverse)
library(qvalue)

min_inv_freq <- 0.03 #The minimum  allele frequency for an inversion to be used in the petiolaris subspecies
species_list <- c("annuus","argophyllus","petpet","petfal")
species_list_capital  <- c("Annuus","Argophyllus","Petiolaris","Petiolaris")

gea_scores <- tibble(species=character(),threshold=character(),
                     inv_positive=numeric(),inv_count=numeric(),
                     snp_positive=numeric(),snp_count=numeric(),
                     pvalue=numeric())

for (i in 1:4){
  chosen_species <- species_list[i]
  chosen_species_capital <- species_list_capital[i]
  col_names <- read_tsv(paste0("/media/owens/Copper/wild_gwas/env_associations/for_Greg/inversions_HA412/",
                        chosen_species,"/column.names.varout.txt"))
  env_associations <- read_tsv(paste0("/media/owens/Copper/wild_gwas/env_associations/for_Greg/var_out_",chosen_species,"_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.txt"),
                               col_names = c("chr_pos","chr","pos",colnames(col_names)[2:length(colnames(col_names))])) %>%
    select(-chr,-pos) %>%
    separate(chr_pos, c("chr","pos"),"__") %>%
    mutate(pos = as.numeric(pos))
  
  
  
  inversion_regions <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt")
  chr_lengths <- read_tsv("Ha412HO.chrlengths.txt")
  
  inversions <- chr_lengths %>% 
    select(chr,end) %>%
    rename(chr_len = end) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(inversion_regions, ., by=c("chr"="chr")) %>%
    # Add a cumulative position of each SNP
    arrange(chr, start) %>%
    mutate( startcum=start+tot,endcum=end+tot) %>%
    filter(species == chosen_species_capital) 
  

  

  ###Subset the genome to remove windows in inversions.
  env_associations %>% select(chr, pos) -> env_sites
  
  env_sites %>%
    mutate(window = floor(pos/1000000),
           window_start = window*1000000,
           window_end = (window+1)*1000000) %>%
    select(chr,window,window_start,window_end) %>%
    unique() -> all_windows
  
  windows_in_inversions <- tibble(chr=character(), window=character(), 
                                  window_start=numeric(),window_end=numeric())
  for (i in 1:nrow(inversions)){
    chosen_chr = inversions$chr[i]
    chosen_start = inversions$start[i]
    chosen_end = inversions$end[i]
    all_windows %>%
      filter(chr == chosen_chr,window_end > chosen_start, window_start < chosen_start) -> left_side
    all_windows %>%
      filter(chr == chosen_chr,window_end > chosen_end, window_start < chosen_end) -> right_side 
    all_windows %>%
      filter(chr == chosen_chr,window_end < chosen_end, window_start > chosen_start) -> middle
    all_windows %>%
      filter(chr == chosen_chr,window_end > chosen_end, window_start < chosen_start) -> surrounding
    
    
    windows_in_inversions <- rbind(windows_in_inversions, left_side, right_side, middle, surrounding)
    
  }
  
  windows_outside_inversions <- all_windows %>%
    anti_join(., windows_in_inversions)
  
  env_associations %>%
    mutate(window = floor(pos/1000000),
           window_start = window*1000000,
           window_end = (window+1)*1000000) %>%
    semi_join(.,windows_outside_inversions)  %>%
    select(-window,-window_start,-window_end)-> env_associations_noninv
  
  env_associations_noninv %>%
    gather(.,env_variable,bayesfactor,-chr,-pos) %>%
    mutate(type = "snp") -> snps_long
  
  variables <- unique(snps_long$env_variable)
  
  inv_associations <- read_delim(paste0("gea/",chosen_species,"_HA412_remapped/baypass/tables/varout_",chosen_species,"_HA412_remapped_Baypass_noinv_matrix_25_vars.tab.table"),
                                 delim=" ") %>%
    rename(inversion_pos = inversion_ID)
    
  if (chosen_species_capital == "Petiolaris"){
    inv_frequencies <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.petinvfreq.txt") %>%
      filter(tolower(species) == chosen_species) %>%
      mutate(inversion_pos = paste0(chr,"_",mds))
    inner_join(inv_associations, inv_frequencies) %>%
      filter(freq > min_inv_freq,freq < (1-min_inv_freq)) %>%
      select(-freq,-species,-chr,-mds) -> inv_associations
  }
  
  inv_associations %>%
    gather(.,env_variable,bayesfactor, -inversion_pos) %>%
    mutate(type = "inv",chr="NA",pos=1) %>%
    rename(snp_id = inversion_pos) -> invs_long
  
  invs_long %>%
    mutate(sig = case_when(bayesfactor > 10 ~ "sig",
                           TRUE ~ "nonsig")) %>%
    group_by(snp_id,sig) %>%
    summarise(n=n()) %>%
    mutate(freq = n / sum(n)) %>%
    filter(sig == "nonsig") ->
    invs_long_epistasis
  
  snps_long %>%
    mutate(snp_id = paste0(chr,"_",pos)) %>%
    mutate(sig = case_when(bayesfactor > 10 ~ "sig",
                           TRUE ~ "nonsig")) %>%   
    group_by(snp_id,sig,.drop=F) %>%
    summarise(n=n()) %>% 
    mutate(freq = n / sum(n)) %>%
    filter(sig == "nonsig")-> snps_long_epistasis
  
  #########
  #Check if they have one or more positive associations
  #########
  snps_long_epistasis %>%
    group_by(freq) %>%
    summarise(n=n()) %>%
    mutate(type = "snp") %>%
    mutate(any_sig = case_when(freq < 1 ~ "Yes",
                               TRUE ~ "No")) %>%
    group_by(any_sig) %>%
    summarize(count=sum(n)) -> tmp.snp
  
  invs_long_epistasis %>%
    group_by(freq) %>%
    summarise(n=n()) %>%
    mutate(type = "inv") %>%
    mutate(any_sig = case_when(freq < 1 ~ "Yes",
                               TRUE ~ "No")) %>%
    group_by(any_sig) %>%
    summarize(count=sum(n)) -> tmp.inv
  test_result <- prop.test(c(tmp.inv$count[which(tmp.inv$any_sig == "Yes")], tmp.snp$count[which(tmp.snp$any_sig == "Yes")]),
                           c(sum(tmp.inv$count),sum(tmp.snp$count)))
  
  tmp.tibble <- tibble(species=chosen_species,threshold="1",inv_positive=tmp.inv$count[which(tmp.inv$any_sig == "Yes")],
                       inv_count=sum(tmp.inv$count),snp_positive=tmp.snp$count[which(tmp.snp$any_sig == "Yes")],
                       snp_count=sum(tmp.snp$count),pvalue=test_result$p.value)
  gea_scores <- rbind(gea_scores,tmp.tibble)
  
  #########
  #Check if they have two or more positive associations
  #########
  min_scores = (length(variables)-1)/length(variables) - 0.01
  snps_long_epistasis %>%
    group_by(freq) %>%
    summarise(n=n()) %>%
    mutate(type = "snp") %>%
    mutate(any_sig = case_when(freq < min_scores ~ "Yes",
                               TRUE ~ "No")) %>%
    group_by(any_sig) %>%
    summarize(count=sum(n)) -> tmp.snp
  
  invs_long_epistasis %>%
    group_by(freq) %>%
    summarise(n=n()) %>%
    mutate(type = "inv") %>%
    mutate(any_sig = case_when(freq < min_scores ~ "Yes",
                               TRUE ~ "No")) %>%
    group_by(any_sig) %>%
    summarize(count=sum(n)) -> tmp.inv
  test_result <- prop.test(c(tmp.inv$count[which(tmp.inv$any_sig == "Yes")], tmp.snp$count[which(tmp.snp$any_sig == "Yes")]),
                           c(sum(tmp.inv$count),sum(tmp.snp$count)))
  tmp.tibble <- tibble(species=chosen_species,threshold="2+",inv_positive=tmp.inv$count[which(tmp.inv$any_sig == "Yes")],
                       inv_count=sum(tmp.inv$count),snp_positive=tmp.snp$count[which(tmp.snp$any_sig == "Yes")],
                       snp_count=sum(tmp.snp$count),pvalue=test_result$p.value)
  gea_scores <- rbind(gea_scores,tmp.tibble)
  
}

write_tsv(gea_scores,"gea/var_out_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.LNHenrichment.txt")
gea_scores %>%
  mutate(percent_inv = inv_positive/inv_count,
         percent_snp = snp_positive/snp_count) %>%
  gather(type,value,percent_inv:percent_snp) %>%
  mutate(species = case_when(species=="annuus" ~ "H. annuus",
                             species=="argophyllus"~"H. argophyllus",
                             species=="petpet"~"H. petiolaris petiolaris",
                             species=="petfal"~"H. petiolaris fallax")) %>%
  mutate(threshold=case_when(threshold=="1"~"1+ variables",
                             threshold=="2+"~"2+ variables")) %>%
  mutate(sig_dot=case_when(pvalue < 1e-20 ~ "**",
                           pvalue < 0.05 ~ "*",
                           TRUE ~ " ")) %>%
  group_by(species,threshold) %>%
  mutate(max_prop = max(value))%>%
  ggplot(.,aes(x=species,y=value,fill=type)) + geom_bar(stat="identity",position="dodge") +
  geom_text(aes(x=species,y=max_prop+0.02,label=sig_dot)) +
  facet_wrap(~threshold) +
  ylab("Proportion associated") + xlab("Species") +
  theme_bw() + 
  scale_fill_manual(values=c("grey","black"),name="Variant Type",
                    labels=c("LNH","SNP")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("gea/var_out_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.LNHenrichment.pdf")
