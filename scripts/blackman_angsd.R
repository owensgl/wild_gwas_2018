library(tidyverse)
library(ggthemes)
library(ggExtra)
library(ggpubr)

data <- read_tsv("/media/owens/Copper/wild_gwas_2018/annuus/Annuus.tranche90.snp.env.90.bi.500.mdsoutlierfst.blackmanangsd.txt")
data %>% pull(mds) %>% unique() -> mds_list

info <- read_tsv("/media/owens/Copper/wild_gwas_2018/annuus/Angds_all_annuus.beagle.info.txt") %>%
  rename(sample = Beagle_ID)

pdf("Annuus.tranche90.snp.env.90.bi.500.mdsoutlierfst.blackmanangsd.scatter.pdf")
for (mds_chosen in mds_list){
  
  data %>%
    inner_join(.,info) %>%
    filter(mds == mds_chosen) %>%
    mutate(total = (g00 + g01 + g11)*2,
           het= g01 / (g00 + g11 + g01),
           percent_1 = (g01 + (g11*2))/total) -> tmp
  print(
  ggscatterhist(
    tmp, x = "percent_1", y = "het",
    palette = "Set1",
    color = "Type", # comment out this and last line to remove the split by species
    margin.plot = "histogram", # I'd suggest removing this line to get density plots
    margin.params = list(fill = "Type", color = "black", size = 0.2),
    title = paste(mds_chosen)
    
  ) 
  )
}
dev.off()

pdf("Annuus.tranche90.snp.env.90.bi.500.mdsoutlierfst.blackmanangsd.hist.pdf",height=6,width=24)
data %>%
  inner_join(.,info) %>%
  #filter(mds == mds_chosen) %>%
  mutate(total = (g00 + g01 + g11)*2,
         het= g01 / (g00 + g11 + g01),
         percent_1 = (g01 + (g11*2))/total) %>%
  ggplot(.,aes(x=percent_1)) +
  geom_histogram(aes(fill=Type)) +
  theme_few() + 
  facet_grid(Type~mds,scales="free_y") +
  scale_fill_brewer(palette = "Set1") +
  ylab("Count") + xlab("Percent inversion 1 alleles") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

write_tsv(data %>%
            inner_join(.,info) %>%
            mutate(total = (g00 + g01 + g11)*2,
                   het= g01 / (g00 + g11 + g01),
                   percent_1 = (g01 + (g11*2))/total),"Annuus.tranche90.snp.env.90.bi.500.mdsoutlierfst.blackmanangsd.output.txt")


library(cowplot) 
# Main plot


