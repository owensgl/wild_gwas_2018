library(tidyverse)




data <- read_tsv("/media/owens/Copper/wild_gwas/annuus/Annuus.tranche90.snp.env.90.bi.remappedHa412HO.Ha412HOChr01.neg1.fst.txt.gz",
            n_max=680865)

pdf("Chr1_allele_freq.pdf")
data %>%
  filter(pos >6199548, pos < 9699260) %>% 
  ggplot(.,aes(x=INV_0_freq,y=INV_2_freq)) +
  geom_point(,alpha=0.1) + 
  geom_density_2d() +
  theme_bw()

data %>%
  filter(pos >6199548, pos < 9699260) %>% 
  mutate(type = case_when(INV_0_freq > 0.99 & INV_2_freq < 0.01 ~ "fix_dif",
                          INV_0_freq > 0.99 & INV_2_freq < 0.99 ~ "inv2_seg",
                          INV_0_freq < 0.01 & INV_2_freq < 0.99 ~ "inv2_seg",
                          INV_0_freq < 0.99 & INV_2_freq < 0.01 ~ "inv0_seg",
                          INV_0_freq < 0.99 & INV_2_freq > 0.99 ~ "inv0_seg",
                          TRUE ~ "both_seg")) %>%
  ggplot(.,aes(type)) + geom_bar() +
  theme_bw()

window_size <- 100000
data %>%
  filter(pos >6199548, pos < 9699260) %>% 
  mutate(type = case_when(INV_0_freq > 0.99 & INV_2_freq < 0.01 ~ "fix_dif",
                          INV_0_freq > 0.99 & INV_2_freq < 0.99 ~ "inv2_seg",
                          INV_0_freq < 0.01 & INV_2_freq < 0.99 ~ "inv2_seg",
                          INV_0_freq < 0.99 & INV_2_freq < 0.01 ~ "inv0_seg",
                          INV_0_freq < 0.99 & INV_2_freq > 0.99 ~ "inv0_seg",
                          TRUE ~ "both_seg")) %>%
  mutate(window = floor(pos/window_size)) %>%
  group_by(window,type,.drop = FALSE) %>%
  summarize(count=n()) %>%
  group_by(window) %>%
  mutate(total = sum(count)) %>%
  ggplot(.,aes(x=window/10,y=count/total)) + geom_line(aes(color=type)) +
  theme_bw() + 
  scale_color_brewer(palette = "Set1") +
  xlab("Mb")



dev.off()


