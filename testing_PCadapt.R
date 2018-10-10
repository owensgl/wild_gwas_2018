install.packages("devtools")
devtools::install_github("bcm-uga/pcadapt")
library(pcadapt)


tmp.pcadapt <- read.pcadapt("/media/owens/Copper/wild_gwas_2018/argophyllus/Argophyllus.tranche90.snp.gwas.90.bi.bed", type = "bed")
sites <- read_tsv("/media/owens/Copper/wild_gwas_2018/argophyllus/Argophyllus.tranche90.snp.gwas.90.bi.bim", 
                  col_names = c("chr","x","x1","pos","x2","x3")) %>% select(chr, pos)
x <- pcadapt(input = tmp.pcadapt, K = 20) 


sites$pvalues <- -log10(x$pvalues)


sites %>%
  filter(chr == "HanXRQChr17") %>%
  ggplot(.,aes(x=pos,y=pvalues)) + geom_point()
