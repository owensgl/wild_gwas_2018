library(tidyverse)
library(ggthemes)
directory <- "/media/owens/Copper/wild_gwas_2018/petiolaris"

pet29.30 <- read_tsv(paste(directory,"/Petiolaris.tranche90.snp.gwas.90.bi.pet29pet30.fst.txt",sep=""))
pet29.30$comparison <- "pet29.30"

pet29.31 <- read_tsv(paste(directory,"/Petiolaris.tranche90.snp.gwas.90.bi.pet29pet31.fst.txt",sep=""))
pet29.31$comparison <- "pet29.31"

pet30.31 <- read_tsv(paste(directory,"/Petiolaris.tranche90.snp.gwas.90.bi.pet30pet31.fst.txt",sep=""))
pet30.31$comparison <- "pet30.31"

pet30.50 <- read_tsv(paste(directory,"/Petiolaris.tranche90.snp.gwas.90.bi.pet30pet50.fst.txt",sep=""))
pet30.50$comparison <- "pet30.50"

fst <- rbind(pet29.30,pet29.31, pet30.31,pet30.50)

window_size <- 1000000
pdf("Pet_30_fst.v1.pdf",height=8,width=12)

fst %>%
  mutate(window = floor(pos/window_size)*window_size) %>% 
  group_by(comparison, chr, window) %>%
  summarize(win_fst = sum(FstNum)/sum(FstDenom)) %>% 
  ggplot(aes(x=window/1000000,y=win_fst)) + geom_line(aes(group=comparison,color=comparison),alpha=0.8) +
  facet_wrap(~chr,scales="free_x") +theme_few() +
  scale_color_brewer(palette = "Set1",name="Location") +
  xlab("MBase") + ylab("Fst") 
dev.off()
