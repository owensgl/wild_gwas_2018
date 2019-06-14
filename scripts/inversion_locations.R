#This is for making a figure of inversion locations
library(tidyverse)
library(statebins)
chr_sizes <- read_tsv("Ha412HO.chrlengths.txt")

inv_locations <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt")

inv_locations %>%
  separate(chr, c("ref","chr_n"), "Chr") %>%
  mutate(chr_n = as.numeric(chr_n) %>%
  filter(species != "Niveus") %>%
  mutate(species_vert = case_when(species == "Annuus" ~ 0.66,
                                  species == "Argophyllus" ~ 0.33,
                                  species == "Petiolaris" ~ 0)) -> inv_locations_formatted
  
  
pdf("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.pdf",height=5,width=5)
chr_sizes %>%
  separate(chr, c("ref","chr_n"), "Chr") %>%
  mutate(chr_n = as.numeric(chr_n)) %>%
  ggplot(.) + 
  geom_rect(data=inv_locations_formatted,aes(xmin=start/1000000,xmax=end/1000000,ymin=chr_n+species_vert-0.5,ymax=chr_n+species_vert+0.33-0.5,fill=species),color=NA) +
  geom_rect(aes(xmin=start/1000000,xmax=end/1000000,ymin=chr_n-0.5,ymax=chr_n+0.5),color="black",fill=NA) +
  scale_fill_brewer(palette = "Set1",name="Species") +
  scale_y_reverse(breaks=c(1:17)) +
  xlab("MB") + ylab("Ha412HO Chromosome") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  theme_classic()
dev.off()
  



