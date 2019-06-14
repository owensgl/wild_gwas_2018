library(tidyverse)
#This is for Z-normalization the traits for GWAS
species <- "petpet"
traits <- c( "TLN", "LIR", "DTF","DTF_corrected","DTF_deduced",
             "late_early","Stem_diamater_at_flowering",
             "Plant_height_at_flowering", "Primary_branches",
             "SLA_mm2_per_mg", "Leaf_total_N", "Leaf_total_C",
             "Disk_diameter", "Petal_length", "Petal_width",
             "Guides_3_petals", "Stem_colour", "Leaf_Area",
             "Leaf_Maximum_height", "Phyllaries_length", 
             "Phyllaries_width", "Seed_area", "Seed_HW_ratio")


phenotypes <- read_tsv(paste0("gwas/",species,"/",toupper(species),"_all_phenotypes_jan_2019.txt")) %>%
  rename(sample = FID)

for (i in 1:length(traits)){
  chosen_trait <- traits[i]
  tmp <- phenotypes %>% select(sample) %>% mutate(family = sample)
  if (any(names(phenotypes) == chosen_trait)){
    print(paste("printing",chosen_trait))
    tmp$trait <- phenotypes %>% select(!!chosen_trait) %>% pull() %>% scale()
    write_tsv(tmp, paste("gwas/",species,"/",chosen_trait,".znorm.txt",sep=""),col_names = F)
  }

}
