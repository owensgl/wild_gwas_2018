#Plots of Texanus allele frequencies

library(tidyverse)

sample_info <- read_tsv("resources/sample_info_file_all_samples_2018_12_07.tsv") %>%
  rename(sample=name,pop=population)
pop_info <- read_tsv("pop_loc.txt")


inversions <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt") %>%
  filter(spe == "annuus")


all_genos <- tibble(sample=character(),triangle_genotype=numeric(),pop=character(),species=character(),
                    inversion=character())
for (i in 1:11){
  

chr_data <- read_tsv(paste0("MDS_outliers/Ha412HO/annuus/Ha412HO_inv.v3.pcasites.Ha412HOChr",
                            sprintf("%02d",inversions$chr[i]),".",
                            inversions$direction[i],inversions$mds[i],".genotypes.txt"))

usa <- map_data('state')
states <- map_data("state")
target_state <- map_data('state')
lat_range <- c(25, 50)
long_range <- c(-125,-93)
pie_size <- 0.4


chr_data %>%
  inner_join(.,sample_info) %>%
  filter(species== "Ann") %>%
  inner_join(.,pop_info) %>%
  group_by(pop,lat,long) %>%
  summarize(mean_loci=mean(triangle_genotype)) -> plotting_data

chr_data %>%
  inner_join(.,sample_info) %>%
  filter(species== "Ann") %>%
  inner_join(.,pop_info) %>%
  mutate(inversion = paste0("Chr",sprintf("%02d",inversions$chr[i]),".",
                            inversions$direction[i],inversions$mds[i])) %>%
  select(sample,triangle_genotype,pop,species,inversion) -> tmp
all_genos <- rbind(all_genos,tmp)
  

}
all_genos %>%
  inner_join(pop_info) %>%
  group_by(inversion,pop,lat,long) %>%
  summarize(mean_loci=mean(triangle_genotype)) -> plotting_data

ggplot(target_state, aes(long, lat)) +
  geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
  coord_quickmap() +
  geom_point(data=plotting_data,
             aes(x=long, y=lat, color=mean_loci/2), 
             size=3,alpha=0.8) +
  scale_color_viridis_c(name="Allele_frequency") +theme_bw() +
  xlab("Longitude") + ylab("Latitude") +
  scale_x_continuous(limits = long_range, expand = c(0, 0)) +
  scale_y_continuous(limits = lat_range, expand = c(0, 0)) +
  facet_wrap(~inversion)

texanus_populations <- read_tsv("resources/annuus.texanusfst.grouplist.txt",col_names=c("pop","texanus"))
all_genos %>%
  inner_join(.,pop_info) %>%
  inner_join(.,texanus_populations) %>%
  group_by(inversion,texanus) %>%
  summarize(mean_loci=mean(triangle_genotype)) %>%
  ggplot(.,aes(x=texanus,y=mean_loci/2)) + geom_bar(stat="identity",position="dodge") +
  facet_wrap(~inversion)
