library(tidyverse)
library(grid)
library(gridExtra)
library(RColorBrewer)


inv_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt")
inv_ranges <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.inversions.regions.v1.txt")
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T) %>% rename(sample=name)
pop_loc <- read_tsv("pop_loc_allnum.txt")
genome_mask <- read_delim("Han412-HO.maskedbases.txt",col_names = c("chr","window","end","masked"),delim=" ") 
window_size <- 1000000
min_sites <- 500

for (n in 1:11){

  
  
  #Load all the information for the specific inversion
  chosen_species <- pull(inv_list[n,2])
  chosen_noncapital_species <- pull(inv_list[n,1])
  chosen_chr <- paste("Ha412HOChr",sprintf("%02d",pull(inv_list[n,4])),sep="")
  chosen_mds <- paste(pull(inv_list[n,6]),pull(inv_list[n,5]),sep="")
  inv_ranges %>%
    filter(species == chosen_species, chr == chosen_chr, mds == chosen_mds) %>%
    mutate(start = floor(start/window_size)*window_size,
           end = ceiling(end/window_size) * window_size) -> chosen_ranges
  
  info <- read_tsv(paste("MDS_outliers/Ha412HO/",chosen_noncapital_species,"/Ha412HO_inv.jan09.pcasites.",chosen_chr,".",chosen_mds,".genotypes.txt",sep=""))
  
  
  #Load in IBD distances for each population
  distance_files <- list.files(paste("Inv_ibd/Ha412HO/",chosen_noncapital_species,"/",sep=""),paste(chosen_chr,".",chosen_mds,".*multispecies.ibd.txt",sep=""))
  if (length(distance_files) == 0){
    next
  }
  full_dist <- tibble(chr=character(),window=character(),sample1=character(),
                      sample2=character(),ibd_dist=numeric(),sites=numeric(),
                      test_sample=character())
  for (i in 1:length(distance_files)){
    target_sample <- strsplit(distance_files[i],"\\.")[[1]][6]
    dist <- read_tsv(paste("Inv_ibd/Ha412HO/",chosen_noncapital_species,"/",distance_files[i],sep="")) %>%
      mutate(test_sample = target_sample)
    full_dist <- rbind(full_dist, dist)
    
  }
  
  
  full_dist %>%
    mutate(inversion_window = "normal") %>%
    mutate(inv_chr = chosen_chr) -> full_dist
  
  for (i in 1:nrow(chosen_ranges)){
    chosen_start <- pull(chosen_ranges[i,2])
    chosen_end <- pull(chosen_ranges[i,3])
    full_dist %>%
      mutate(inversion_window = case_when(window >= chosen_start & window < chosen_end & chr == chosen_chr ~ "inverted",
                                          TRUE ~ inversion_window)) -> full_dist
  }
  
  relatives <- unique(full_dist$test_sample)
  
  pdf(paste("Inv_ibd/Ha412HO/",chosen_noncapital_species,"/Ha412HO_inv.jan09.",chosen_chr,".",chosen_mds,".ibd_relatives.pdf",sep=""),
      height=12,width=18)
  for (i in 1:length(relatives)){
    relative <- relatives[i]
    inv_1 <- full_dist %>% filter(test_sample == relative) %>%
      filter(sample1 != relative) %>% pull(sample1) %>% unique()
    inv_2 <- full_dist %>% filter(test_sample == relative) %>%
      filter(sample2 != inv_1) %>% pull(sample2) %>% unique()
    #Figure out which sample has which inversion genotype
    inv_1_genotype <- info %>% filter(sample == inv_1) %>% pull(triangle_genotype)
    if (inv_1_genotype == 0){
      geno_0_sample <- inv_1
      geno_2_sample <- inv_2
    }else{
      geno_0_sample <- inv_2
      geno_2_sample <- inv_1
    }
    relative_species <- labels %>% filter(sample == relative) %>% pull(species)
    
    
    plot_box <- full_dist %>%
      filter(test_sample == relative) %>%
      filter(sites > min_sites) %>%
      mutate(comparison = paste(sample1,sample2,sep="-")) %>% 
      filter(sample1 == relative ) %>%
      group_by(chr, window, inversion_window) %>%
      summarize(dif = ibd_dist[which(comparison == paste(relative,"-",geno_0_sample,sep=""))] - 
                  ibd_dist[which(comparison == paste(relative,"-",geno_2_sample,sep=""))]) %>%
      inner_join(.,genome_mask) %>%
      mutate(dif_norm = dif/(window_size-masked)) %>%
      ggplot(.,aes(x=inversion_window,y=dif_norm,fill=inversion_window)) +
      geom_boxplot() + theme_bw() +
      #coord_cartesian(ylim=c(-500,500)) +
      scale_fill_brewer(palette = "Set1",name="Genomic\nregion") +
      ylab(paste(relative_species,"_Inv0 distance - ",relative_species,"_Inv2 distance",sep="")) +
      xlab("Genomic region")
    
    plot_line <- full_dist %>%
      filter(test_sample == relative) %>%
      filter(chr == chosen_chr) %>% 
      filter(sites > min_sites) %>%
      mutate(comparison = paste(sample1,sample2,sep="-")) %>% 
      filter(sample1 == relative ) %>%
      group_by(chr, window, inversion_window) %>%
      summarize(dif = ibd_dist[which(comparison == paste(relative,"-",geno_0_sample,sep=""))] - 
                  ibd_dist[which(comparison == paste(relative,"-",geno_2_sample,sep=""))]) %>%
      inner_join(.,genome_mask) %>%
      mutate(dif_norm = dif/(window_size-masked)) %>%
      ggplot(.,aes(x=window/window_size, y=dif_norm)) + geom_line() +
      geom_segment(data=chosen_ranges,
                   aes(x=start/window_size,xend=end/window_size,y=0,yend=0),color="red",alpha=0.7,size=2) +
      theme_bw() + xlab("MB") + 
      ylab(paste(relative_species,"_Inv0 distance - ",relative_species,"_Inv2 distance",sep=""))
    
    
    plot_bar <- full_dist %>%
      filter(test_sample == relative) %>%
      filter(sites > min_sites) %>%
      mutate(comparison = paste(sample1,sample2,sep="-")) %>% 
      group_by(chr, window, inversion_window) %>%
      summarize(min_comp = case_when( ((ibd_dist[which(comparison ==  paste(relative,"-",geno_0_sample,sep=""))] < ibd_dist[which(comparison ==  paste(relative,"-",geno_2_sample,sep=""))]) &
                                         (ibd_dist[which(comparison ==  paste(relative,"-",geno_0_sample,sep=""))] < ibd_dist[which(comparison ==  paste(inv_1,"-",inv_2,sep=""))])) ~ paste("Rel","_Inv0",sep=""),
                                      ((ibd_dist[which(comparison ==  paste(relative,"-",geno_2_sample,sep=""))] < ibd_dist[which(comparison ==  paste(relative,"-",geno_0_sample,sep=""))]) &
                                         (ibd_dist[which(comparison ==  paste(relative,"-",geno_2_sample,sep=""))] < ibd_dist[which(comparison == paste(inv_1,"-",inv_2,sep=""))])) ~ paste("Rel","_Inv2",sep=""),
                                      TRUE ~ "Inv0_Inv2")) %>%
      ggplot(.,aes(x=inversion_window,fill=min_comp)) +
      geom_bar(position="fill") +
      theme_bw() +
      scale_fill_brewer(palette = "Set2",name="Sister pair") +
      ylab("Proportion of windows") +
      xlab("Genomic region")
    
    print(
      grid.arrange( plot_line,plot_box,plot_bar,
                    layout_matrix = rbind(c(1, 1),
                                          c(2, 3)),
                    top = textGrob(paste(relative_species,"compared to", chosen_species,chosen_chr,chosen_mds),gp=gpar(fontsize=18,font=1))
      )
    )
    
  }
  dev.off()
}

