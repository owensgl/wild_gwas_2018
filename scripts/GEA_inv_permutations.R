library(tidyverse)
library(qvalue)

env_associations <- read_tsv("/media/owens/Copper/wild_gwas/env_associations/for_Greg/var_out_annuus_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.txt",
                             col_names = c("chr_pos","chr","pos", "Latitude_e","Longitude_e","Elevation_e","MAT_e",
                                           "MWMT_e","MCMT_e","TD_e","MAP_e","MSP_e","AHM_e",
                                           "SHM_e","DD_0_e","DD5_e","DD_18_e","DD18_e","NFFD_e",
                                           "bFFP_e","eFFP_e","FFP_e","PAS_e","EMT_e","EXT_e","Eref_e",
                                           "CMD_e","MAR_e","RH_e")) %>%
  select(-chr,-pos) %>%
  separate(chr_pos, c("chr","pos"),"__") %>%
  mutate(pos = as.numeric(pos))

inversion_regions <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt")
chr_lengths <- read_tsv("Ha412HO.chrlengths.txt")

ann_inversions <- chr_lengths %>% 
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
  filter(species == "Annuus") 

ann_inversions_length <- ann_inversions %>%
  mutate(length = end - start) %>%
  group_by(chr, mds) %>%
  summarize(length=round(sum(length)/1000000))

ann_inversion_list <- ann_inversions %>% select(chr, mds) %>% unique()

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
for (i in 1:nrow(ann_inversions)){
  chosen_chr = ann_inversions$chr[i]
  chosen_start = ann_inversions$start[i]
  chosen_end = ann_inversions$end[i]
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

inv_columns <- read_tsv(paste("/media/owens/Copper/wild_gwas/env_associations/for_Greg/inversions_HA412/annuus/column.names.varout.annuus.txt",sep=""))
invs <- read_tsv(paste("/media/owens/Copper/wild_gwas/env_associations/for_Greg/inversions_HA412/annuus/var_out_annuus_inversions_HA412_vcf_all_covariates_BayPass_IS_STD.txt",sep=""),
                col_names = colnames(inv_columns))

invs %>%
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

snps_long_epistasis %>% filter(sig == "nonsig") -> snps_long_epistasis
snps_long_epistasis %>%
  mutate(any_sig = case_when(freq < 1 ~ "Yes",
                             TRUE ~ "No")) %>%
  group_by(any_sig) %>%
  summarise(n=n()) %>%
  mutate(freq = n / sum(n)) -> snps_long_sig_sites

snps_long_epistasis %>%
  group_by(freq) %>%
  summarise(n=n()) %>%
  mutate(type = "snp") %>%
  mutate(any_sig = case_when(freq < 1 ~ "Yes",
                             TRUE ~ "No")) %>%
  group_by(any_sig) %>%
  summarize(count=sum(n))

invs_long_epistasis %>%
  group_by(freq) %>%
  summarise(n=n()) %>%
  mutate(type = "inv") %>%
  mutate(any_sig = case_when(freq < 1 ~ "Yes",
                             TRUE ~ "No")) %>%
  group_by(any_sig) %>%
  summarize(count=sum(n))

prop.test(c(5,176744),c(11,2945413))



