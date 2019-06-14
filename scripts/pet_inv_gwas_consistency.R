library(tidyverse)

petpet <- read_tsv(paste0("gwas/petpet/Ha412HO_inv.v3.allgwas.nostr.txt")) %>%
  mutate(species = "petpet")

petfal <- read_tsv(paste0("gwas/petfal/Ha412HO_inv.v3.allgwas.nostr.txt")) %>%
  mutate(species = "petfal")

rbind(petpet,petfal) %>%
  #filter(p < 0.5) %>%
  mutate(inv_id = paste(chr,".",direction,mds)) %>%
  group_by(inv_id, trait,trait_types) %>%
  mutate(count = sum(!is.na(estimate))) %>% filter(count == 2) %>%
  summarize(signal=case_when(estimate[1] > 0 & estimate[2] > 0 ~ "match" ,
                          estimate[1] < 0 & estimate[2] < 0 ~ "match" ,
                          TRUE ~ "opposite")) %>% 
  ggplot(., aes(x = as.factor(inv_id), fill = factor(signal))) +
  geom_bar(position="fill")
