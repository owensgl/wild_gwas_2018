#Testing for hardyweinberg for all inversions
library(tidyverse)
library(HWxtest)
library(scatterpie)


directory <- "MDS_outliers/Ha412HO/"
all_species <- c("annuus","argophyllus","petpet","petfal","niveus")
chosen_species <- "annuus"
chosen_date <- "jan09"
filter <- "pcasites"
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T) %>% rename(sample=name)
pop_loc <- read_tsv("pop_loc_allnum.txt")

for (n in 1:length(all_species)){
  chosen_species <- all_species[n]
  pdf(paste(directory,chosen_species,"/Ha412HO_inv.",chosen_date,".",filter,".HW-Utest.pdf",sep=""),
      width=15,height=15)
  
  if (chosen_species == "annuus"){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(25, 50)
    long_range <- c(-125,-93)
    pie_size <- 0.4
  }else if(chosen_species == "argophyllus" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- subset(states, region %in% c("texas"))
    lat_range <- c(25, 30)
    long_range <- c(-99,-96)
    pie_size <- 0.2
  }else if(chosen_species== "petpet" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(30, 43)
    long_range <- c(-115,-95)
    pie_size <- 0.4
  }else if(chosen_species == "petfal" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(30, 40)
    long_range <- c(-115,-100)
    pie_size <- 0.4
  }else if(chosen_species == "petiolaris" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(30, 43)
    long_range <- c(-115,-95)
    pie_size <- 0.4
  }else if(chosen_species == "niveus" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(30, 34)
    long_range <- c(-109,-102)
    pie_size <- 0.2
  }
  
  
  
  
  vcf <- read_tsv(paste(directory,chosen_species,"/Ha412HO_inv.",chosen_date,".",filter,".vcf",sep=""),comment="##")
  
  
  #Make maps of proportions to check genotyping
  
  vcf %>% 
    mutate(locus = paste(`#CHROM`,ID,sep=".")) %>% 
    select(-`#CHROM`,-POS,-ID,-REF,-ALT,-QUAL,-FILTER,-INFO,-FORMAT) %>%
    gather(key = sample, value = value, -locus) %>%
    mutate(genotype = case_when(value == '0/0' ~ "0",
                                value == '0/1' ~ "1",
                                value == '1/1' ~ "2",
                                TRUE ~ "NA")) -> genotypes
  
  genotypes <- labels %>% select(sample, population) %>%
    inner_join(.,genotypes) 
  
  genotypes <- pop_loc %>% rename(population = pop) %>%
    inner_join(.,genotypes) 
  
  print(
    ggplot(target_state, aes(long, lat)) +
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
  )
  print(
  genotypes %>%
    group_by(sample, genotype) %>%
    tally() %>%
    spread(., genotype, n,fill=0) %>%
    ggplot(.,aes(x=`1`)) + geom_bar(aes(y = (..count..)/sum(..count..))) +
    ylab("Proportion") +
    xlab("Heterozygous inversions") +
    theme_bw()
  )
  
  vcf %>% 
    mutate(locus = paste(`#CHROM`,ID,sep=".")) %>% 
    select(-`#CHROM`,-POS,-ID,-REF,-ALT,-QUAL,-FILTER,-INFO,-FORMAT) %>%
    gather(key = sample, value = value, -locus) %>% 
    spread(key = locus, value = value) -> vcf
  
  labels %>% 
    select(sample,population) %>% 
    inner_join(.,vcf) -> vcf
  
  
  rownames(vcf) <- vcf$sample
  
  vcf %>% select(-sample) %>% rename(pop = population) -> vcf
  
  hw_result <- hwx.test(vcf, statName="U")
  
  hw_utest <- hwdf(hw_result) 
  
  hw_utest$tmp <- rownames(hw_utest)
  
  hw_utest <- hw_utest %>% separate(tmp, c("pop","chr","mds"),"\\.") %>%
    inner_join(.,pop_loc)  %>%mutate(locus = paste(chr,mds,sep="."))
  
  
  
  test_result <- t.test(hw_utest$`obs-U`)
  t_test_result <- tibble(locus="All",t=test_result$statistic,
                          df = test_result$parameter, p.value= test_result$p.value,
                          mean=test_result$estimate,min=test_result$conf.int[1],
                          max=test_result$conf.int[2])
  for (locus_n in 1:length(unique(hw_utest$locus))){
    n_pops <- hw_utest %>% filter(locus == unique(hw_utest$locus)[locus_n]) %>% pull(`obs-U`) %>% length()
    
    if(n_pops <= 1){
      tmp <- tibble(locus=unique(hw_utest$locus)[locus_n],t=NA,
                    df = NA, p.value= NA,
                    mean=NA,min=NA,
                    max=NA)
    }else{
      test_result <- t.test(hw_utest %>% filter(locus == unique(hw_utest$locus)[locus_n]) %>% pull(`obs-U`))
      tmp <- tibble(locus=unique(hw_utest$locus)[locus_n],t=test_result$statistic,
                    df = test_result$parameter, p.value= test_result$p.value,
                    mean=test_result$estimate,min=test_result$conf.int[1],
                    max=test_result$conf.int[2])
    }
    
    t_test_result <- rbind(t_test_result, tmp)
    
  }
  write_tsv(t_test_result,paste(directory,chosen_species,"/Ha412HO_inv.",chosen_date,".",filter,".HW-Utest.txt",sep=""))
  print(
    ggplot(t_test_result,aes(x=locus,y=mean)) +
      geom_point() + geom_errorbar(aes(ymin=min,ymax=max),width=0.2) +
      theme_bw() +
      geom_hline(yintercept = 0,linetype="dotted") +
      theme(axis.text.x = element_text(angle = 45, hjust=1)) +
      ylab("Mean U score") +
      xlab("Locus")
  )
  print(
    ggplot(target_state, aes(long, lat)) +
      geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
      coord_quickmap() +
      geom_point(data=hw_utest,
                 aes(x=long, y=lat,color=`obs-U`), 
                 cols=c("0","1","2"), size=4, alpha=.9) +
      scale_color_distiller(name="U score",palette = "Spectral") +theme_bw() +
      xlab("Longitude") + ylab("Latitude") +
      scale_x_continuous(limits = long_range, expand = c(0, 0)) +
      scale_y_continuous(limits = lat_range, expand = c(0, 0)) +
      facet_wrap(~locus) +
      ggtitle("U score")
  )
  dev.off()
  
}
  