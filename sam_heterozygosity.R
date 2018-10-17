#SAM heterozygosity
library(tidyverse)
library(ggthemes)


het <- read_tsv("/media/owens/Copper/wild_gwas_2018/annuus_sam/Annuus.tranche90.snp.fullsam.90.bi.hetwindows.txt")


pdf("Annuus.tranche90.snp.fullsam.90.bi.hetwindows.pdf")
for (sample_name in sort(unique(het$name))){
  print(
  het %>%
    filter(name == sample_name) %>%
    filter(counts > 100) %>%
    ggplot(.,aes(x=start,y=perc_het)) +
    geom_line() + 
    facet_wrap(~chr) +
    coord_cartesian(ylim=c(0,1)) +
    theme_bw() +
    ggtitle(sample_name)
  )
}
dev.off()

het %>%
  group_by(name) %>%
  summarize(total_counts = sum(counts),total_het = sum(hets),
            perc_het = total_het/total_counts) %>%
  arrange(desc(perc_het)) %>%
  write_tsv(.,"Annuus.tranche90.snp.fullsam.90.bi.hetwindows.total.txt")



  