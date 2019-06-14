#Visualizing IBD distance across inversions
library(tidyverse)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(scatterpie)
inv_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt")
inv_ranges <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.inversions.regions.v1.txt")
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T) %>% rename(sample=name)
pop_loc <- read_tsv("pop_loc_allnum.txt")

window_size <- 1000000


for (n in 1:nrow(inv_list)){

  full_dist <- tibble(chr=character(),window=numeric(),sample1=character(),
                      sample2=character(),ibd_dist=numeric(),sites=numeric(),
                      triangle_genotype1=character(),triangle_genotype2=character(),
                      inversion_geno=character(),inversion_window=character(),
                      inv_chr=character(),
                      mds=character(),species=character(),pop=character())
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
  info %>% select(sample, triangle_genotype) %>%
    rename(sample1 = sample, triangle_genotype1 = triangle_genotype) -> info1
  
  info %>% select(sample, triangle_genotype) %>%
    rename(sample2 = sample, triangle_genotype2 = triangle_genotype) -> info2
  
  
  
  #Load in IBD distances for each population
  ibd_files <- list.files(paste("Inv_ibd/Ha412HO/",chosen_noncapital_species,"/",sep=""),paste(chosen_chr,".",chosen_mds,".*.txt",sep=""))
  if (length(ibd_files) == 0){
    next
  }
  
  
  for (x in 1:length(ibd_files)){
    chosen_pop <- str_split(ibd_files[x],"\\.")[[1]][6] #Need to update with new names
    dist <- read_tsv(paste("Inv_ibd/Ha412HO/",chosen_noncapital_species,"/",chosen_species,".tranche90.snp.",chosen_chr,".",chosen_mds,".",chosen_pop,".remappedtarget.ibd.txt",sep="")) 
    #Check for overlap of windows with inversions, kinda complicates because of multiple blocks
    
    dist %>%
      inner_join(., info1) %>%
      inner_join(., info2) %>% 
      mutate(inversion_geno = case_when(triangle_genotype1 == triangle_genotype2 ~ "Same",
                                        TRUE ~ "Different"))  %>%
      mutate(inversion_window = "normal") %>%
      mutate(inv_chr = chosen_chr) -> dist
    
    for (i in 1:nrow(chosen_ranges)){
      chosen_start <- pull(chosen_ranges[i,2])
      chosen_end <- pull(chosen_ranges[i,3])
      dist %>%
        mutate(inversion_window = case_when(window >= chosen_start & window < chosen_end & chr == chosen_chr ~ "inverted",
                                            TRUE ~ inversion_window)) -> dist
    }
    dist %>%
      mutate(pop = chosen_pop,species=chosen_species,mds=chosen_mds) -> dist
    
    full_dist <- rbind(full_dist,dist)
  }
  if (nrow(full_dist) == 0){next}
  
  line_plot <- full_dist %>%
    filter(chr == chosen_chr) %>%
    mutate(comparison = paste(sample1,sample2,sep="-")) %>%
    ggplot(.) + 
    geom_line(aes(x=window/window_size,y=ibd_dist, color=inversion_geno,group=comparison),alpha=0.5) +
    theme_bw() + scale_color_brewer(palette = "Set1",name="Inversion\ngenotype") +
    facet_wrap(~pop,nrow=10) +
    geom_segment(data=chosen_ranges,
                 aes(x=start/window_size,xend=end/window_size,y=0,yend=0),size=2) +
    xlab("MB") + ylab("IBD distance")
  
  
  colourCount = length(unique(full_dist$pop))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  

  
  pop_plot <- full_dist %>%
    mutate(comparison = paste(sample1,sample2,sep="-")) %>%
    filter(inversion_geno == "Same") %>% 
    ggplot(.,aes(x=inversion_window,y=ibd_dist,fill=as.factor(triangle_genotype1))) +
    geom_boxplot() +
    facet_wrap(~pop) +
    theme_bw() +
    scale_fill_brewer(palette = "Set1",name="Inversion\nGenotype") +
    ylab("IBD distance") + xlab("Genomic region") +
    theme(axis.text.x = element_text(angle = 60, hjust=1))
  
  
  dif_plot1 <- full_dist %>%
    mutate(comparison = paste(sample1,sample2,sep="-")) %>%
    filter(inversion_geno == "Same") %>% 
    group_by(chr, window,pop, inversion_window) %>%
    mutate(count = n()) %>%
    filter(count == 2) %>%
    summarize(dif = ibd_dist[which(triangle_genotype1 == 2)] - ibd_dist[which(triangle_genotype1 == 0)]) %>%
    ggplot(.,aes(x=inversion_window,y=dif,fill=pop),fill="white") +
    geom_boxplot() + theme_bw() +
    ylab("IBD[genotype == 2] - IBD[genotype == 0]") +
    xlab("Genomic region") +
    scale_fill_manual(values=getPalette(colourCount),name="Population") +
    theme(legend.position="none") +
    geom_hline(yintercept = 0,linetype="dotted")
  
  
  #Plotting whether it supports 0 or 2 as derived.
  full_dist %>%
    filter(inversion_geno == "Same") %>% 
    mutate(comparison = paste(sample1,sample2,sep="-")) %>% 
    group_by(chr,window,inversion_window,pop) %>%
    mutate(count = n()) %>%
    filter(count == 2) %>%
    summarize(dif = ibd_dist[which(triangle_genotype1 == 2)] - ibd_dist[which(triangle_genotype1 == 0)]) %>%
    group_by(inversion_window,pop) %>%
    summarize(mean_dif = mean(dif,na.rm=T)) %>%
    filter(inversion_window == "inverted") %>%
    mutate(sign= case_when(mean_dif <0 ~ "2_derived",
                           mean_dif > 0 ~ "0_derived")) %>%
    group_by(sign) %>%
    count() %>%
    ungroup()-> derived_signal
    
  if((derived_signal %>% filter(sign == "2_derived") %>% nrow()) == 0){
    tmp <- tibble(sign="2_derived",n=as.numeric(0))
    derived_signal <- rbind(derived_signal, tmp)
  }else if ((derived_signal %>% filter(sign == "0_derived") %>% nrow()) == 0){
    tmp <- tibble(sign="0_derived",n=as.numeric(0))
    derived_signal <- rbind(derived_signal, tmp)
  }
    
    signal_plot <- ggplot(derived_signal,aes(x=sign,y=n,fill=sign)) + geom_bar(stat="identity") +
      theme_bw() + scale_fill_brewer(palette = "Set2") +
      ylab("Population counts") +
      xlab("Signal")
    
    if (derived_signal %>% filter(sign == "2_derived") %>% pull(n) > 
        derived_signal %>% filter(sign == "0_derived") %>% pull(n)){
      derived <- 2
    }else if (derived_signal %>% filter(sign == "2_derived") %>% pull(n) < 
              derived_signal %>% filter(sign == "0_derived") %>% pull(n)){
      derived <- 0
    }else{
      derived <- "NA"
    }
  
  
  dif_plot2 <- full_dist %>%
    mutate(comparison = paste(sample1,sample2,sep="-")) %>% 
    group_by(chr,window,inversion_window,pop) %>%
    summarize(dif = mean(ibd_dist[which(inversion_geno == "Different")]) - mean(ibd_dist[which(inversion_geno == "Same")]))  %>%
    ggplot(.,aes(x=inversion_window,y=dif,fill=pop)) +
    geom_boxplot() +
    ylab("IBD distance") + xlab("Genomic region") +
    theme_bw() +
    ylab("IBD[genotype == different] - IBD[genotype == same]") +
    xlab("Genomic region") +
    scale_fill_manual(values=getPalette(colourCount),name="Population") +
    theme(legend.position="none") + 
    geom_hline(yintercept = 0,linetype="dotted")
  
  full_dist %>%
    mutate(comparison = paste(sample1,sample2,sep="-")) %>% 
    group_by(chr,window,inversion_window,pop) %>%
    summarize(dif = mean(ibd_dist[which(inversion_geno == "Different")]) - mean(ibd_dist[which(inversion_geno == "Same")])) %>%
    group_by(inversion_window,pop) %>%
    summarize(mean_dif = mean(dif,na.rm=T))
  
  
  pdf(paste("Inv_ibd/Ha412HO/",chosen_noncapital_species,"/Ha412HO_inv.jan09.",chosen_chr,".",chosen_mds,".ibd_dist.pdf",sep=""),height=15,width=15)
  
  grid.arrange( line_plot,
                layout_matrix = rbind(c(1, 1),
                                      c(1, 1)),
                top = textGrob(paste(chosen_species,chosen_chr,chosen_mds),gp=gpar(fontsize=18,font=1))
  )
  
  grid.arrange( pop_plot, dif_plot1,dif_plot2,signal_plot,
    layout_matrix = rbind(c(1, 2),
                          c(1, 3),
                          c(1, 4)),
    top = textGrob(paste(chosen_species,chosen_chr,chosen_mds),gp=gpar(fontsize=18,font=1))
  )
  
  if (chosen_noncapital_species == "annuus"){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(25, 50)
    long_range <- c(-125,-93)
    pie_size <- 0.4
  }else if(chosen_noncapital_species == "argophyllus" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- subset(states, region %in% c("texas"))
    lat_range <- c(25, 30)
    long_range <- c(-99,-96)
    pie_size <- 0.2
  }else if(chosen_noncapital_species== "petpet" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(30, 43)
    long_range <- c(-115,-95)
    pie_size <- 0.4
  }else if(chosen_noncapital_species == "petfal" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(30, 40)
    long_range <- c(-115,-100)
    pie_size <- 0.4
  }else if(chosen_noncapital_species == "petiolaris" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(30, 43)
    long_range <- c(-115,-95)
    pie_size <- 0.4
  }else if(chosen_noncapital_species == "niveus" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(30, 34)
    long_range <- c(-109,-102)
    pie_size <- 0.2
  }
  
  
  
  
  vcf <- read_tsv(paste("MDS_outliers/Ha412HO/",chosen_noncapital_species,"/Ha412HO_inv.jan09.pcasites.vcf",sep=""),comment="##")
  
  
  #Make maps of proportions to check genotyping
  
  vcf %>% 
    mutate(locus = paste(`#CHROM`,ID,sep=".")) %>% 
    select(-`#CHROM`,-POS,-ID,-REF,-ALT,-QUAL,-FILTER,-INFO,-FORMAT) %>%
    gather(key = sample, value = value, -locus) %>%
    mutate(genotype = case_when(value == '0/0' ~ "0",
                                value == '0/1' ~ "1",
                                value == '1/1' ~ "2",
                                TRUE ~ "NA")) %>%
    filter(locus == paste(chosen_chr,chosen_mds,sep="."))-> genotypes
  
  genotypes <- labels %>% select(sample, population) %>%
    inner_join(.,genotypes) 
  
  genotypes <- pop_loc %>% rename(population = pop) %>%
    inner_join(.,genotypes) 
  
  
  if (derived == 2){
    genotypes %>%
      mutate(genotype = case_when(genotype == 0 ~ "Ancestral",
                                  genotype == 1 ~ "Het",
                                  genotype == 2 ~ "Derived")) -> genotypes
  } else if (derived == 0){
    genotypes %>%
      mutate(genotype = case_when(genotype == 2 ~ "Ancestral",
                                  genotype == 1 ~ "Het",
                                  genotype == 0 ~ "Derived")) -> genotypes
  }
  if (derived == 0 | derived == 2){
    map_plot <- ggplot(target_state, aes(long, lat)) +
      geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
      coord_quickmap() +
      geom_scatterpie(data=genotypes %>% 
                        group_by(population, lat, long, locus, genotype) %>% 
                        tally() %>%
                        spread(., genotype, n,fill=0),
                      aes(x=long, y=lat, r=pie_size), 
                      cols=c("Ancestral","Het","Derived"), color=NA, alpha=.8) +
      scale_fill_manual(name="Genotype",values=c(brewer.pal(3,"Set1")[1],brewer.pal(3,"Set1")[3],brewer.pal(3,"Set1")[2]),
                        breaks=c("Ancestral", "Het", "Derived"),
                        labels=c("Ancestral", "Het", "Derived")) +theme_bw() +
      xlab("Longitude") + ylab("Latitude") +
      scale_x_continuous(limits = long_range, expand = c(0, 0)) +
      scale_y_continuous(limits = lat_range, expand = c(0, 0)) +
      facet_wrap(~locus) 
  }else{
    
    map_plot <- ggplot(target_state, aes(long, lat)) +
      geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
      coord_quickmap() +
      geom_scatterpie(data=genotypes %>% 
                        group_by(population, lat, long, locus, genotype) %>% 
                        tally() %>%
                        spread(., genotype, n,fill=0),
                      aes(x=long, y=lat, r=pie_size), 
                      cols=c("0","1","2"), color=NA, alpha=.8) +
      scale_fill_brewer(name="Genotype",palette = "Set1") +theme_bw() +
      xlab("Longitude") + ylab("Latitude") +
      scale_x_continuous(limits = long_range, expand = c(0, 0)) +
      scale_y_continuous(limits = lat_range, expand = c(0, 0)) +
      facet_wrap(~locus) 
  }
  print(map_plot)
  
  
  dev.off()
  
}












