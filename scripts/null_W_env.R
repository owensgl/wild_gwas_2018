library(tidyverse)

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

##Build random 10K sites 
env_associations %>%
  mutate(window = floor(pos/1000000),
         window_start = window*1000000,
         window_end = (window+1)*1000000) %>%
  semi_join(.,windows_outside_inversions) -> env_associations_noninv
env_associations %>%
  sample_n(.,10000) -> null_10k_set


#Get Z scores for inversions
for (i in 1:nrow(ann_inversion_list)){
  chosen_chr <- ann_inversion_list$chr[i]
  chosen_mds <- ann_inversion_list$mds[i]
  env_associations_inv <- env_associations %>% filter(chr == "NANA") #Create empty tibble
  for( n in 1:nrow(ann_inversions)){
    if (chosen_chr != ann_inversions$chr[n]){next}
    if (chosen_mds != ann_inversions$mds[n]){next}
    env_associations %>%
      filter(chr == chosen_chr, pos > ann_inversions$start[n], pos < ann_inversions$end[n]) -> tmp
    env_associations_inv <- rbind(env_associations_inv, tmp)
  }
  inv_scores <- env_associations_inv[,3] %>% pull()
  
  length_inv_scores <- length(inv_scores)
  
  null_10k_scores  <- null_10k_set[,3] %>% pull()
  length_null_10k_scores <- 10000
  w_test <- wilcox.test(inv_scores,null_10k_scores)
  
  Z <- (2*w_test$statistic - (length_inv_scores*length_null_10k_scores))/
    sqrt(length_inv_scores * length_null_10k_scores*(length_inv_scores + length_null_10k_scores + 1) / 3)
  w_test
  
  
}


#Get scores for nor null_window_sets
n_windows <- 29
perm_repeats <- 1000
perm_vec <- rep("", times=perm_repeats)
for (i in 1:1000){
  perm_windows <- windows_outside_inversions %>%
    sample_n(.,n_windows) 
  null_window_scores <- env_associations_noninv %>%
    semi_join(.,perm_windows) %>% .[,3] %>% pull() 
  length_null_window_scores <- length(null_window_scores)
  
  null_10k_scores  <- null_10k_set[,3] %>% pull()
  length_null_10k_scores <- 10000
  w_test <- wilcox.test(null_window_scores,null_10k_scores)
  
  Z <- (2*w_test$statistic - (length_null_window_scores*length_null_10k_scores))/
    sqrt(length_null_window_scores * length_null_10k_scores*(length_null_window_scores + length_null_10k_scores + 1) / 3)
  w_test
  perm_vec[i] <- Z
}


  hist(as.numeric(perm_vec))

