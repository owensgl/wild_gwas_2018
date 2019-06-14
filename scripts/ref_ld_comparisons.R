#For plotting comparative LD
library(tidyverse)

window_size = 5000
XRQ <- read_tsv("/media/owens/Childs/wild_gwas/annuus/Annuus.tranche90.snp.env.90.bi.localLD.geno.ld") %>%
  rename(r2 = `R^2`) %>%
  mutate(dist_window = floor((POS2 - POS1)/window_size))

XRQ %>% 
  group_by(dist_window) %>%
  select(r2) %>%
  filter(!is.na(r2)) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n()))) -> summarized_XRQ

summarized_XRQ$ref <- "XRQ"



HA412 <- read_tsv("/media/owens/Childs/wild_gwas/annuus/Annuus.tranche90.snp.env.90.bi.remappedHa412HO.localLD.geno.ld") %>%
  rename(r2 = `R^2`) %>%
  mutate(dist_window = floor((POS2 - POS1)/window_size))



HA412 %>% 
  group_by(dist_window) %>%
  select(r2) %>%
  filter(!is.na(r2)) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n()))) -> summarized_HA412

summarized_HA412$ref <- "Ha412"

pdf("Reference_LD_comparisons.pdf")
summarized_HA412 %>%
  rbind(summarized_XRQ) %>%
  filter(dist_window < 40) %>%
  ggplot(.,aes(x=dist_window*window_size,y=mean,color=ref)) + 
  geom_point() + geom_errorbar(aes(ymin=mean-se,ymax=mean+se)) +
  theme_bw() + scale_color_brewer(palette = "Set1",name="Reference genome") +
  ylab("Mean r2") + xlab("bp")
dev.off()


