#Plotting fst outliers, but actual calls.

library(tidyverse)
library(ggthemes)
library(forcats)
library(ggrastr)

fst_calls <- read_tsv("/media/owens/Copper/wild_gwas_2018/annuus_sam/Annuus.tranche90.snp.fullsam.90.bi.500.mdsoutlierfst.sitecalls.txt")
samples_genotyped <- read_tsv("Annuus.tranche90.snp.env.90.bi.500.mds_cluster_genotyped.txt") %>%
  mutate(mds_coord = gsub("_","-",mds_coord))

labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T)


fst_calls %>%
  rename(name = sample, site_genotype = genotype) %>%
  inner_join(., labels) %>% 
  inner_join(., samples_genotyped) -> fst_calls


pdf("Annuus.tranche90.snp.fullsam.90.bi.500.mdsoutlierfst.sitecalls.pdf")
for (mds_chosen in sort(unique(fst_calls$mds_coord))){
  print(
    fst_calls %>%
      filter(mds_coord == mds_chosen) %>%
      ggplot(.,aes(y=fct_reorder(name, genotype),x=as.factor(pos),fill=as.factor(site_genotype))) + 
      geom_tile_rast() + 
      theme_few() + 
      scale_fill_brewer(palette = "Set1",name="Site Genotype") +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      ylab("Sample") + xlab("Position") +
      ggtitle(mds_chosen)
  )
}
dev.off()
