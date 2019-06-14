library(tidyverse)

dtf <- read_tsv("arg_dtf_tmp.txt") %>% rename(sample = individual)

geno <- read_tsv(paste("MDS_outliers/Ha412HO/argophyllus/Ha412HO_inv.jan09.pcasites.Ha412HOChr06.pos4.genotypes.txt",sep="")) %>%
  select(sample, triangle_genotype)

pdf("Arg_dtf_linkage_tmp_pos1.pdf")
geno %>%
  inner_join(dtf) %>%
  ggplot(.,aes(x=as.factor(triangle_genotype),y=DTF)) + geom_boxplot() + geom_jitter() +
  theme_bw() + xlab("Inversion Chr6-end, 64.1% PVE")
dev.off()

geno %>%
  inner_join(dtf) %>%
  lm(DTF ~ triangle_genotype,.) %>% anova(.)


geno <- read_tsv(paste("MDS_outliers/Ha412HO/argophyllus/Ha412HO_inv.jan09.pcasites.Ha412HOChr06.pos4.genotypes.txt",sep="")) %>%
  select(sample, triangle_genotype)

geno %>%
  inner_join(dtf) %>%
  lm(DTF ~ triangle_genotype,.) %>% anova(.)
pdf("Arg_dtf_linkage_tmp_pos1.pdf")
geno %>%
  inner_join(dtf) %>%
  ggplot(.,aes(x=as.factor(triangle_genotype),y=DTF)) + geom_boxplot() + geom_jitter() +
  theme_bw() + xlab("Inversion Chr6-end, 64.1% PVE")
dev.off()