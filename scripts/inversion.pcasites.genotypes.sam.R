#Plotting inversion genotypes in a heatmap
library(tidyverse)


directory <- "MDS_outliers/Ha412HO/"
chosen_species <- "annuus"
chosen_date <- "jan09"
filter <- "pcasites"
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T) %>% rename(sample=name)
pop_loc <- read_tsv("pop_loc_allnum.txt")


vcf <- read_tsv(paste(directory,chosen_species,"/Ha412HO_inv.",chosen_date,".",filter,".vcf",sep=""),comment="##")

vcf %>% 
  mutate(locus = paste(`#CHROM`,ID,sep=".")) %>% 
  select(-`#CHROM`,-POS,-ID,-REF,-ALT,-QUAL,-FILTER,-INFO,-FORMAT) %>%
  gather(key = sample, value = value, -locus) %>%
  mutate(genotype = case_when(value == '0/0' ~ "0",
                              value == '0/1' ~ "1",
                              value == '1/1' ~ "2",
                              TRUE ~ "NA")) -> genotypes

genotypes <- labels %>% select(sample, is_wild, population) %>%
  inner_join(.,genotypes) 



pdf(paste(directory,chosen_species,"/Ha412HO_inv.",chosen_date,".",filter,".nonwild.genotypes.pdf",sep=""),
    width=5,height=30)
genotypes %>%
  filter(is_wild != "wild") %>%
  filter(genotype != "NA") %>%
  ggplot(.,aes(x=locus,y=sample,fill=genotype)) + geom_tile() +
  scale_fill_brewer(palette = "Set1",name="Inversion call") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=5))

dev.off()

pdf(paste(directory,chosen_species,"/Ha412HO_inv.",chosen_date,".",filter,".nonwild.genotypes.MTKH.pdf",sep=""),
    width=7,height=3)
genotypes %>% 
  filter(sample == "SAM053" | sample == "SAM254" | sample == "SAM261") %>%
  mutate(sample = case_when(sample == "SAM053" ~ "HA412HO",
                            sample == "SAM261" ~ "XRQ",
                            TRUE ~ sample)) %>%
  filter(genotype != "NA") %>%
  ggplot(.,aes(x=locus,y=sample,fill=genotype)) + geom_tile() +
  scale_fill_brewer(palette = "Set1",name="Inversion call") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=5))

dev.off()


