library(tidyverse)

#Count petiolaris inversion frequencies for later filtering

sample_info <- read_tsv("resources/sample_info_file_all_samples_wildspecies.tsv") %>%
  rename(sample=name)

inversions<- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt") %>%
  filter(species == "Petiolaris") %>%
  select(chr,mds) %>% unique()

freqs <- tibble(chr=character(),mds=character(),species=character(),freq=numeric())
for (i in 1:nrow(inversions)){
  chosen_chr = inversions$chr[i]
  chosen_mds = inversions$mds[i]
  
  genotypes <- read_tsv(paste0("MDS_outliers/Ha412HO/petiolaris/Ha412HO_inv.v3.pcasites.",chosen_chr,".",chosen_mds,".genotypes.txt"))
  genotypes %>% inner_join(.,sample_info) %>%
    filter(species == "PetPet") %>%
    summarize(derived=sum(triangle_genotype),count=n()*2,
              freq= derived/count) %>%
    pull(freq) -> freq.petpet
  genotypes %>% inner_join(.,sample_info) %>%
    filter(species == "PetFal") %>%
    summarize(derived=sum(triangle_genotype),count=n()*2,
              freq= derived/count) %>%
    pull(freq) -> freq.petfal

  freq_tmp <- tibble(chr=chosen_chr,mds=chosen_mds,species="PetPet",freq=freq.petpet)
  freqs <- rbind(freqs,freq_tmp)
  freq_tmp <- tibble(chr=chosen_chr,mds=chosen_mds,species="PetFal",freq=freq.petfal)
  freqs <- rbind(freqs,freq_tmp)
  
}
write_tsv(freqs,"MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.petinvfreq.txt")
freqs %>%
  ggplot(.,aes(x=paste0(chr,"_",mds),y=freq,fill=species)) + geom_bar(stat="identity",position="dodge") +

