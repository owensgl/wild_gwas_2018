library(tidyverse)
#This is for Z-normalization the traits for GWAS
species_specific_list <- c("annuus","argophyllus","petpet","petfal")
species_general_list <- c("Annuus","Argophyllus","Petiolaris","Petiolaris")
tag_list <- c("gwas","gwas","petpet","petfal")

for (i in 1:4){
  species_specific <- species_specific_list[i]
  species_general <- species_general_list[i]
  tag <- tag_list[i]
  samplelist <- read_tsv(paste0("/media/owens/Copper/wild_gwas/gwas/",species_specific,"/",species_general,
                         ".tranche90.snp.",tag,".90.bi.samplelist.txt"),
                         col_names=c("sample")) %>%
    filter(sample != "ARG0300")
  phenotypes <- read_tsv(paste0("resources/",tolower(species_general),"_all_phenotypes_May2019.txt")) %>%
    rename(sample = FID) %>% semi_join(.,samplelist)
    
  traits <- colnames(phenotypes)[3:length(colnames(phenotypes))]
  
  for (x in 1:length(traits)){
    chosen_trait <- traits[x]
    tmp <- phenotypes %>% select(sample) %>% mutate(family = sample)
    print(paste("printing",chosen_trait))
    tmp$trait <- phenotypes %>% select(!!chosen_trait) %>% pull() %>% as.numeric() %>% scale()
    write_tsv(tmp, paste("gwas/",species_specific,"/",chosen_trait,".znorm.txt",sep=""),col_names = F)
    tmp$trait <- phenotypes %>% select(!!chosen_trait) %>% pull() %>% as.numeric() 
    write_tsv(tmp, paste("gwas/",species_specific,"/",chosen_trait,".txt",sep=""),col_names = F)
    
  }
  
  
}
