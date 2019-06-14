#Plotting fst outliers, but actual calls. for pet and neg F1 maps

library(tidyverse)
library(ggthemes)
library(forcats)
library(ggrastr)
library(gridExtra)
library(scatterpie)

#Load in all inversion genotypes
folder <- "MDS_outliers"
chosen_species <- "petiolaris"
chosen_species_file <- "petiolaris_syn"
chosen_species_abbreviation <- c("PetPet","PetFal")
filtering <- "pcasites"
prefix <- "Ha412HO_inv.v3"
base_directory <- "/media/owens/Copper/wild_gwas_2018/petmap"
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T)
pop_loc <- read_tsv("pop_loc_allnum.txt") %>% rename(population = pop)
pca_genotypes <- pca_genotypes %>% filter(species == chosen_species)
inversion_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt") %>%
  rename(chr_n = chr) %>%
  mutate(chr = paste0("Ha412HOChr",sprintf("%02d", chr_n)),
         mds = paste0(direction,mds))
inversions <- inversion_list %>% filter(spe == chosen_species) %>%
  select(chr,mds) %>% unique()

pdf("Kate_F1_maps_inversiongenotypes.pdf",height=5,width=9)

for (n in (1:nrow(inversions))){
  chosen_chr <- inversions[n,1] %>% pull()
  chosen_mds <- inversions[n,2] %>% pull()
  

  
  
  
  fst_calls <- read_tsv(paste(base_directory,"/",prefix,".",chosen_species_file,".", chosen_chr,".",chosen_mds,".",
                              filtering,".genotypes.txt.gz",sep=""),guess_max= 100000) %>%
    filter(genotype == "00" | genotype == "01" | genotype == "11")
  
  min_depth_genotyping <- 5
  min_depth_visualizing <- 2
  fst_calls %>%
    filter(depth >= min_depth_genotyping) %>%
    group_by(sample,genotype) %>%
    count() %>%
    filter(!is.na(genotype)) %>%
    group_by(genotype) %>%
    count()-> check
  if (nrow(check) < 3){next}
  
  
  fst_calls %>%
    filter(depth >= min_depth_genotyping) %>%
    group_by(sample,genotype) %>%
    count() %>%
    filter(!is.na(genotype)) %>%
    spread(genotype, n,fill=0) %>%
    mutate(total =  (`00` + `01` + `11`)*2,
           perc_1 = ((`11`*2) + `01`)/total,
           het = (`01`*2)/total,
           left_formula = ((-(2/3)*perc_1) + (2/3)),
           right_formula = (((2/3)*perc_1))) %>%
    mutate(fst_call = case_when(perc_1 < 0.5  & left_formula >= het ~ 0,
                                perc_1 < 0.5  & left_formula < het ~ 1,
                                perc_1 >= 0.5  & right_formula >= het ~ 2,
                                perc_1 >= 0.5  & right_formula < het ~ 1)) -> fst_samples
  
  

  fst_samples %>%
    separate(sample,c("pop","sample"),"_") -> fst_samples

  
  print(
    fst_samples %>%
      filter(total > 1) %>%
      ggplot(.,aes(x=perc_1,y=het)) +
      geom_segment(aes(x=0,xend=0.5,y=0,yend=1),color="grey") +
      geom_segment(aes(x=0.5,xend=1,y=1,yend=0),color="grey") +
      geom_segment(aes(x=0,xend=1,y=0,yend=0),color="grey") +
      stat_bin2d(binwidth=c(0.05, 0.05)) +
      theme_few() +
      ylab("Heterozygosity") +
      xlab("Proportion `1` allele") + 
      facet_wrap(~pop) +
      ggtitle(paste(chosen_species,chosen_chr,chosen_mds)) +
      scale_fill_viridis()
  )
  
}
dev.off()