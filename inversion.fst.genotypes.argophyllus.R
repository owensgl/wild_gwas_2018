#Plotting fst outliers, but actual calls.

library(tidyverse)
library(ggthemes)
library(forcats)
library(ggrastr)
library(gridExtra)

prefix <- "Argophyllus.tranche90.snp.gwas.90.bi.500"
input_directory <- "/media/owens/Copper/wild_gwas_2018/argophyllus"

fst_calls <- read_tsv(paste(input_directory, "/",prefix, ".mdsoutlierfst.sitecalls.txt",sep=""))

labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T)


min_depth = 5 #Minimum number of reads required

fst_calls %>%
  rename(name = sample, site_genotype = genotype) %>%
  inner_join(., labels)   -> fst_calls

fst_calls %>%
  filter(depth >= min_depth) %>%
  group_by(name, mds_coord,  site_genotype) %>%
  count() %>%
  filter(!is.na(site_genotype)) %>%
  spread(site_genotype, n,fill=0) %>%
  mutate(total =  (`00` + `01` + `11`)*2,
         perc_1 = ((`11`*2) + `01`)/total,
         het = (`01`*2)/total) %>%
  mutate(fst_call = case_when(perc_1 < 0.1 ~ 0,
                              perc_1 > 0.9 ~ 2,
                              TRUE ~ 1)) -> fst_samples




pdf(paste(prefix, ".mdsoutlierfst.calls.pdf",sep=""),width=8,height=6)
for (mds_chosen in sort(unique(fst_samples$mds_coord))){
  if(fst_samples %>%filter(total > 40, mds_coord == mds_chosen) %>% nrow() == 0){next}
  fst_samples %>%
    inner_join(.,labels) %>%
    filter(total > 40, mds_coord == mds_chosen) %>%
    ggplot(.,aes(x=perc_1,y=het)) +
    geom_segment(aes(x=0,xend=0.5,y=0,yend=1),color="grey") +
    geom_segment(aes(x=0.5,xend=1,y=1,yend=0),color="grey") +
    geom_segment(aes(x=0,xend=1,y=0,yend=0),color="grey") +
    geom_point() +
    theme_few() +
    ylab("Heterozygosity") +
    xlab("Proportion `1` allele") +
    facet_wrap(~species)-> p1
  
  fst_samples %>%
    inner_join(.,labels) %>%
    filter(total > 40, mds_coord == mds_chosen) %>%
    ggplot(.,aes(perc_1)) +
    geom_histogram() +
    theme_few() +
    ylab("Count") +
    xlab("Proportion `1` allele") +
    facet_wrap(~species,scales="free_y")  ->p2
  print(  
    grid.arrange(p1, p2, nrow = 2,top=mds_chosen)
  )
}
dev.off()




pdf(paste(prefix, ".mdsoutlierfst.sitecalls.pdf",sep=""),width=5,height=10)
for (mds_chosen in sort(unique(fst_calls$mds_coord))){
  if (fst_calls %>%
      filter(depth >= min_depth) %>%
      filter(mds_coord == mds_chosen,!is.na(chr)) %>% nrow() == 0){next()}
  print(
    fst_calls %>%
      filter(depth >= min_depth) %>%
      filter(mds_coord == mds_chosen) %>%
      inner_join(., fst_samples) %>% 
      ggplot(.,aes(y=fct_reorder(name, perc_1),x=as.factor(pos),fill=as.factor(site_genotype))) + 
      geom_tile_rast() + 
      theme_few() + 
      scale_fill_brewer(palette = "Set1",name="Site Genotype") +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      ylab("Sample") + xlab("Position") +
      ggtitle(paste(mds_chosen, "min depth =", min_depth)) +
      facet_wrap(~species,nrow=3,scales="free_y")
  )
}
dev.off()


