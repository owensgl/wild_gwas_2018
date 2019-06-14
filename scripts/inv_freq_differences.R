library(tidyverse)
library(viridis)
prefix <- "Ha412HO_inv.jan09.pcasites.annuus"
chosen_species <- "annuus"
chosen_abbreviation <- "ANN"
chosen_noncapital_species <- "annuus"
chosen_capital_species <- "Annuus"
inv_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt") %>%
  filter(spe == chosen_noncapital_species)
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T) %>% rename(sample=name)
pop_info <- read_tsv("pop_loc_allnum.txt") %>% rename(population = pop)
  
  inversion_genotypes <- tibble(sample=character(),genotype=numeric(),mds=character()) 
  n_inv <- nrow(inv_list %>% filter(spe == !!chosen_species))
  for (n in 1:n_inv){
    chosen_noncapital_species <- pull(inv_list[n,1])
    chosen_chr <- paste("Ha412HOChr",sprintf("%02d",pull(inv_list[n,4])),sep="")
    chosen_mds <- paste(pull(inv_list[n,6]),pull(inv_list[n,5]),sep="")
    inversion_n <- paste("inv_",n,sep="")
    
    info <- read_tsv(paste("MDS_outliers/Ha412HO/",chosen_noncapital_species,"/Ha412HO_inv.jan09.pcasites.",chosen_chr,".",chosen_mds,".genotypes.txt",sep="")) %>%
      select(sample, triangle_genotype) %>%
      rename( genotype = triangle_genotype) %>%
      mutate(mds = paste(chosen_chr,".",chosen_mds,sep=""))
    inversion_genotypes <- rbind(inversion_genotypes,info)
  }
  
  
  usa <- map_data('state')
  states <- map_data("state")
  target_state <- map_data('state')
  lat_range <- c(25, 50)
  long_range <- c(-125,-93)
  pie_size <- 0.4
  
  inversion_genotypes %>%
    filter(mds != "Ha412HOChr01.neg1", mds != "Ha412HOChr05.pos1") %>%
    inner_join(.,labels) %>%
    group_by(mds) %>%
    mutate(ref_genotype = genotype[which(sample == "SAM053")]) %>%
    mutate(ref_distance = abs(genotype - ref_genotype)) %>%
    group_by(sample,population) %>%
    filter(species == "Ann",is_wild == "wild",!is.na(population) ) %>%
    summarize(mean_ref_distance = sum(ref_distance)) %>%
    inner_join(.,pop_info) %>%
    group_by(population,lat,long) %>%
    summarize(pop_mean_ref_distance = mean(mean_ref_distance)) -> pop_ref_distance
  
  pdf(paste(prefix,"inversion_refdist.v2.pdf" ,sep=""),height=12,width=20)
  
  inversion_genotypes %>%
    filter(mds != "Ha412HOChr01.neg1", mds != "Ha412HOChr05.pos1") %>%
    inner_join(.,labels) %>%
    group_by(mds) %>%
    mutate(ref_genotype = genotype[which(sample == "SAM053")]) %>%
    mutate(ref_distance = abs(genotype - ref_genotype)) %>%
    group_by(population,mds) %>%
    filter(species == "Ann",is_wild == "wild",!is.na(population) ) %>%
    summarize(mean_ref_distance = mean(ref_distance)/2) %>%
    group_by(population) %>%
    mutate(overall_ref_distance = sum(mean_ref_distance)) %>%
    ggplot(.,aes(x=fct_reorder(population,overall_ref_distance),y=mds,fill=mean_ref_distance)) + 
    geom_tile() +
    scale_fill_viridis(name="Allele freq\ndif from Ha412") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
    dev.off()
  
  pdf(paste(prefix,"inversion_refdist.pdf" ,sep=""))
  ggplot(target_state, aes(long, lat)) +
    geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
    coord_quickmap() +
    geom_point(data=pop_ref_distance,
                    aes(x=long, y=lat, color=pop_mean_ref_distance), alpha=.8,size=4) +
    scale_color_viridis(name="Inversion\ndistance") +theme_bw() +
    xlab("Longitude") + ylab("Latitude") +
    scale_x_continuous(limits = long_range, expand = c(0, 0)) +
    scale_y_continuous(limits = lat_range, expand = c(0, 0)) +
    ggtitle("Number of inversion differences from cultivar")
  
  ggplot(target_state, aes(long, lat)) +
    geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
    coord_quickmap() +
    geom_point(data=pop_ref_distance,
               aes(x=long, y=lat, color=pop_mean_ref_distance), alpha=.8,size=4) +
    scale_color_viridis(name="Inversion\ndistance",limits=c(1,5)) +theme_bw() +
    xlab("Longitude") + ylab("Latitude") +
    scale_x_continuous(limits = long_range, expand = c(0, 0)) +
    scale_y_continuous(limits = lat_range, expand = c(0, 0)) +
    ggtitle("Number of inversion differences from cultivar")
  dev.off()
      
  inversion_genotypes %>%
    filter(mds != "Ha412HOChr01.neg1", mds != "Ha412HOChr05.pos1") %>%
    inner_join(.,labels) %>%
    group_by(mds) %>%
    mutate(ref_genotype = genotype[which(sample == "SAM053")]) %>%
    mutate(ref_distance = abs(genotype - ref_genotype)) %>%
    group_by(sample,population) %>%
    filter(species == "Ann",is_wild == "wild",!is.na(population) ) %>%
    summarize(mean_ref_distance = sum(ref_distance)) %>%
    inner_join(.,pop_info) %>%
    group_by(population,lat,long) %>%
    summarize(pop_mean_ref_distance = mean(mean_ref_distance)) -> pop_ref_distance
  View(pop_ref_distance)
  