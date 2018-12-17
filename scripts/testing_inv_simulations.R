#Testing the slim inversion simulations
library(tidyverse)
library(ggthemes)
slim_test <- read_tsv("/media/owens/Copper/wild_gwas_2018/slim_example_fstcount.txt")


slim_test %>%
  mutate(total =  (count00 + count01 + count11)*2,
         perc_1 = ((count11*2) + count01)/total,
         het = (count01*2)/total) %>%
  ggplot(.,aes(x=perc_1,y=het)) +
  geom_point(aes(color=as.factor(inv_genotype))) +
  geom_segment(aes(x=0,xend=0.5,y=0,yend=1)) +
  geom_segment(aes(x=0.5,xend=1,y=1,yend=0)) +
  geom_segment(aes(x=0,xend=1,y=0,yend=0)) +
  theme_few() +
  ylab("Heterozygosity") +
  xlab("Proportion 1 allele")

slim_test %>%
  mutate(total =  (count00 + count01 + count11)*2,
         perc_1 = ((count11*2) + count01)/total,
         het = (count01*2)/total) %>%
  group_by(inv_genotype) %>%
  summarize(mean_perc_1 = mean(perc_1),count=n()) %>% head()


slim_test %>%
  mutate(total =  (count00 + count01 + count11)*2,
         perc_1 = ((count11*2) + count01)/total,
         het = (count01*2)/total) %>%
  group_by(inv_genotype) %>%
  mutate(mean_perc_1 = mean(perc_1)) %>% head()
  filter(inv_genotype == 1) %>%
  ggplot(.,aes(x=perc_1,y=het)) +
  geom_point(aes(color=as.factor(inv_genotype))) +
  geom_vline(xintercept = mean_perc_1)
  
