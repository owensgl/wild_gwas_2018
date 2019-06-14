#This prints out the summarized LD files. 

library(tidyverse)
library(viridis)


directory <- "/media/owens/Copper/wild_gwas_2018"
species <- "Niveus"
species_abr <- "niveus"
prefix <- "tranche90.snp.gwas.90.bi.remappedHa412HO"

all_data <- tibble(chr=character(),win1=numeric(),win2=numeric(),
                   n=numeric(),mean_r2=numeric(),max_r2=numeric(),
                   max_2_r2=numeric(),r2_min_0.2=numeric(),r2_min_0.8=numeric())
for (chr in sprintf("%02d", seq(1,17))){
  data <-  read_tsv(paste(directory,"/",species_abr,"/",species,".",prefix,".thin100.maf5.Chr",chr,".windows.ld.gz",sep=""))
  all_data <- rbind(all_data,data)
}

pdf(paste("LD/Ha412HO/",species_abr,"/",species,".",prefix,".thin100.maf5.ld.pdf",sep=""),height=15,width=15)
all_data %>%
  mutate(chr = str_replace(chr, "Ha412HO","")) %>%
  filter(n > 50) %>%
  ggplot(.,aes(x=win1/1000000,y=win2/1000000)) + geom_tile(aes(fill=max_2_r2)) +
  theme_bw() + scale_fill_viridis(name="2nd highest R2") +
  xlab("MB") + ylab("MB") +
  facet_wrap(~chr,scale="free")
dev.off()