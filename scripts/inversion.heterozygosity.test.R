#Visualizing heterozygosity across inversions
library(tidyverse)

het <- read_tsv("Annuus.tranche90.snp.chr1inv.remappedtarget.het.txt") %>%
  rename(sample = name)

info <- read_tsv("MDS_outliers/Ha412HO/annuus/Ha412HO_inv.jan09.pcasites.Ha412HOChr01.neg1.genotypes.txt")


window_size = 1000000
het %>% inner_join(.,info) %>%
  mutate(window = floor(start/window_size)*window_size) %>%
  group_by(sample,window,triangle_genotype) %>%
  summarize(hets = sum(hets), counts = sum(counts)) %>%
  filter(counts > 500) %>%
  mutate(inversion = case_when(window >= 6000000 & window < 10000000 ~ "Inverted",
                               TRUE ~ "normal")) %>%
  mutate(het = hets/counts) %>%
  ggplot(.,aes(x=window,y=het,linetype=as.factor(triangle_genotype),color=sample))+ geom_line()  
  #facet_wrap(~sample,nrow=20)


het %>% inner_join(.,info) %>%
  mutate(window = floor(start/window_size)*window_size) %>%
  group_by(sample,window,triangle_genotype) %>%
  summarize(hets = sum(hets), counts = sum(counts)) %>%
  filter(counts > 500) %>%
  mutate(inversion = case_when(window >= 6000000 & window < 10000000 ~ "Inverted",
                               TRUE ~ "normal")) %>%
  mutate(het = hets/counts) %>%
  ggplot(.,aes(x=as.factor(inversion),y=het,fill=as.factor(triangle_genotype))) + 
  geom_boxplot() 
