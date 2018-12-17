library(tidyverse)

tmp <- read_tsv("Argophyllus.tranche90.snp.gwas.90.bi.500.mdsoutlierfst.ha412.txt",
                col_names = c("xrq_chr","xrq_pos","ha412_chr","ha412_pos","bitscore"))


tmp %>%
  filter(ha412_chr == 'Scaffold_273;HRSCAF=505') %>%
  filter(xrq_chr == 'HanXRQChr03') %>%
  ggplot(.,aes(x=xrq_pos/1000000,y=ha412_pos/1000000))+ geom_point() +
  theme_bw() +
  xlab("XRQ MB") + ylab("HA412 MB")

tmp %>%
  filter(ha412_chr == 'Scaffold_4;HRSCAF=6') %>%
  filter(xrq_chr == 'HanXRQChr10') %>% 
  ggplot(.,aes(x=xrq_pos/1000000,y=ha412_pos/1000000))+ geom_point() +
  theme_bw() +
  xlab("XRQ MB") + ylab("HA412 MB")

tmp %>%
  filter(ha412_chr == 'Scaffold_404;HRSCAF=676') %>%
  filter(xrq_chr == 'HanXRQChr06') %>% 
  ggplot(.,aes(x=xrq_pos/1000000,y=ha412_pos/1000000))+ geom_point() +
  theme_bw() +
  xlab("XRQ MB") + ylab("HA412 MB")


blast_test <- read_tsv("/media/owens/Copper/wild_gwas_2018/argophyllus/blast.chr6inv.txt")
bwa_test <- read_tsv("/media/owens/Copper/wild_gwas_2018/argophyllus/bwa.chr6inv.txt")

full_join(blast_test,bwa_test) %>%
  mutate(chr_comparison = case_when(blast_chr == bwa_chr ~"Same",
                                    TRUE ~ "Different")) %>%
  mutate(pos_comparison = blast_pos - bwa_pos) %>%
  group_by(chr_comparison) %>%
  count()
