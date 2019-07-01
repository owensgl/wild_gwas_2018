library(tidyverse)
library(qvalue)

env_associations <- read_tsv("/media/owens/Copper/wild_gwas/env_associations/for_Greg/var_out_argophyllus_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.txt",
                             col_names = c("chr_pos","chr","pos", "Latitude_e","Longitude_e","Elevation_e","MAT_e",
                                           "MWMT_e","MCMT_e","TD_e","MAP_e","MSP_e","AHM_e",
                                           "SHM_e","DD_0_e","DD5_e","DD_18_e","DD18_e","NFFD_e",
                                           "bFFP_e","eFFP_e","FFP_e","EMT_e","EXT_e","Eref_e",
                                           "CMD_e","MAR_e","RH_e")) %>%
  select(-chr,-pos) %>%
  separate(chr_pos, c("chr","pos"),"__") %>%
  mutate(pos = as.numeric(pos))

inversion_regions <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt")
chr_lengths <- read_tsv("Ha412HO.chrlengths.txt")

arg_inversions <- chr_lengths %>% 
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
  filter(species == "Argophyllus") 

arg_inversions_length <- arg_inversions %>%
  mutate(length = end - start) %>%
  group_by(chr, mds) %>%
  summarize(length=round(sum(length)/1000000))

arg_inversion_list <- arg_inversions %>% select(chr, mds) %>% unique()

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
for (i in 1:nrow(arg_inversions)){
  chosen_chr = arg_inversions$chr[i]
  chosen_start = arg_inversions$start[i]
  chosen_end = arg_inversions$end[i]
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

variable_n <- 3:27 #The columns with environmental variables

inv_z_scores <- tibble(chr=character(),mds=character(),variable=character(),z=numeric())
#Get Z scores for inversions
for (i in 1:nrow(arg_inversion_list)){
  chosen_chr <- arg_inversion_list$chr[i]
  chosen_mds <- arg_inversion_list$mds[i]
  print(paste("Testing",chosen_chr,chosen_mds))
  env_associations_inv <- env_associations %>% filter(chr == "NANA") #Create empty tibble
  for( n in 1:nrow(arg_inversions)){
    if (chosen_chr != arg_inversions$chr[n]){next}
    if (chosen_mds != arg_inversions$mds[n]){next}
    env_associations %>%
      filter(chr == chosen_chr, pos > arg_inversions$start[n], pos < arg_inversions$end[n]) -> tmp
    env_associations_inv <- rbind(env_associations_inv, tmp)
  }
  for (n in variable_n){
    inv_scores <- env_associations_inv[,n] %>% pull()
    
    length_inv_scores <- length(inv_scores)
    
    null_10k_scores  <- null_10k_set[,n] %>% pull()
    length_null_10k_scores <- 10000
    w_test <- wilcox.test(inv_scores,null_10k_scores)
    
    Z <- (2*w_test$statistic - (length_inv_scores*length_null_10k_scores))/
      sqrt(length_inv_scores * length_null_10k_scores*(length_inv_scores + length_null_10k_scores + 1) / 3)
    tmp_z_scores <- tibble(chr=chosen_chr,mds=chosen_mds,variable=colnames(env_associations_inv)[n],z=Z)
    inv_z_scores <- rbind(inv_z_scores,tmp_z_scores)
  }
  
}

perm_z_scores <- tibble(chr=character(),mds=character(),variable=character(),z=numeric())

#Get scores for nor null_window_sets. Looped for each inversion to make window number match (i.e. more windows for larger inversions)
perm_repeats <- 1000

for (x in 1:nrow(arg_inversion_list)){
  chosen_chr <- arg_inversion_list$chr[x]
  chosen_mds <- arg_inversion_list$mds[x]
  chosen_windows <- arg_inversions_length %>% filter(chr == chosen_chr, mds == chosen_mds) %>% pull(length)
  print(paste("Permuting",chosen_chr,chosen_mds))
  for (i in 1:perm_repeats){
    if (i %% 100 == 1){
      print(paste("Permutation",i))
    }
    perm_windows <- windows_outside_inversions %>%
      sample_n(.,chosen_windows) 
    null_windows <- env_associations_noninv %>%
      semi_join(.,perm_windows) 
    length_null_window_scores <- nrow(null_windows)
    
    for (n in variable_n){
      
      null_window_scores <- null_windows %>% .[,n] %>% pull()
      null_10k_scores  <- null_10k_set[,n] %>% pull()
      length_null_10k_scores <- 10000
      
      w_test <- wilcox.test(null_window_scores,null_10k_scores)
      
      Z <- (2*w_test$statistic - (length_null_window_scores*length_null_10k_scores))/
        sqrt(length_null_window_scores * length_null_10k_scores*(length_null_window_scores + length_null_10k_scores + 1) / 3)
      tmp_z_scores <- tibble(chr=chosen_chr,mds=chosen_mds,variable=colnames(env_associations_inv)[n],z=Z)
      perm_z_scores <- rbind(perm_z_scores, tmp_z_scores)
    }
  }
}



#Compare permutation Z scores to inversion Z scores

inv_p_scores <- tibble(chr=character(),mds=character(),variable=character(),empirical_p=numeric())
for (i in 1:nrow(arg_inversion_list)){
  chosen_chr <- arg_inversion_list$chr[i]
  chosen_mds <- arg_inversion_list$mds[i]
  print(paste("Testing",chosen_chr,chosen_mds))
  for (n in variable_n){
    chosen_variable <- colnames(env_associations_inv)[n]
    inv_score <- inv_z_scores %>%
      filter(chr == chosen_chr, mds == chosen_mds, variable == chosen_variable) %>% pull(z)
    perm_z_scores %>% 
      filter(chr == chosen_chr, mds == chosen_mds, variable == chosen_variable) %>% 
      mutate(rank = case_when(z < inv_score ~"Lower",
                              z >= inv_score ~ "Higher")) %>%
      group_by(rank,.drop=FALSE) %>%
      summarize(count=n()) -> ranking
    if (nrow(ranking) == 1 & ranking$rank[1] == "Lower"){
      empirical_p <- 1/(ranking$count[1]+1)
    }else if (nrow(ranking) == 1){
      empirical_p <- 1
    }else{
      empirical_p <- (ranking$count[1]+1)/(ranking$count[2] +ranking$count[1] + 1)
    }
    tmp_p_scores <- tibble(chr=chosen_chr,mds=chosen_mds,variable=colnames(env_associations_inv)[n],empirical_p=empirical_p)
    inv_p_scores <- rbind(inv_p_scores,tmp_p_scores)
  }
}


inv_p_scores    
#Make qvalues for each trait
tmp_p_scores <- inv_p_scores %>% filter(chr == "NANANA") 
for (n in variable_n){
  chosen_variable <- colnames(env_associations_inv)[n]
  inv_p_scores %>% filter(variable == chosen_variable) -> tmp
  tmp$q <- qvalue(tmp$empirical_p, pi0 = 1)$qvalues
  tmp_p_scores <- rbind(tmp_p_scores, tmp)
}
inv_p_scores <- tmp_p_scores   

pdf("var_out_argophyllus_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.null_w.pdf")
inv_p_scores %>%
  mutate(qsig = case_when(q < 0.05 ~ "sig",
                          TRUE ~ "nonsig"),
         pmax = case_when(empirical_p <0.000999002 ~ "max",
                          TRUE ~ "nonmax")) %>%
  mutate(chr_mds = paste0(chr,"_",mds)) %>%
  ggplot(.) + geom_tile(aes(y=variable,x=chr_mds,fill=qsig)) +
  scale_fill_manual(values=c("white","grey")) +
  geom_point(data=inv_p_scores %>%
               mutate(qsig = case_when(q < 0.05 ~ "sig",
                                       TRUE ~ "nonsig"),
                      pmax = case_when(empirical_p <0.000999002 ~ "max",
                                       TRUE ~ "nonmax")) %>%
               mutate(chr_mds = paste0(chr,"_",mds)) %>% 
               filter(pmax == "max"), aes(y=variable, x=chr_mds))  +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

pdf("var_out_argophyllus_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.null_w.hist.pdf")
hist(inv_p_scores$empirical_p)
dev.off()

hist(inv_p_scores$empirical_p)
write_tsv(perm_z_scores, "var_out_argophyllus_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.null_w.perm.txt")
write_tsv(inv_z_scores, "var_out_argophyllus_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.null_w.inv_z.txt")
write_tsv(inv_p_scores, "var_out_argophyllus_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.null_w.txt")


