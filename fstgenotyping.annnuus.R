library(tidyverse)
library(gridExtra)
#This is for calling genotypes by high Fst inversion sites

directory <- "/media/owens/Copper/wild_gwas_2018/annuus_sam"
fst_genotyped_file <- "Annuus.tranche90.snp.fullsam.90.bi.500.mdsoutlierfst.calls.txt"

fst_genotyped <- read_tsv(paste(directory,"/",fst_genotyped_file,sep=""))
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T)

head(fst_genotyped)
pdf("Annuus.tranche90.snp.env.90.bi.500.mdsoutlierfst.calls.pdf",width=8,height=6)
for (mds in unique(fst_genotyped$mds_coord)){
    callable_samples <- fst_genotyped %>%
      mutate(total =  (`00` + `01` + `11`)*2,
             perc_1 = ((`11`*2) + `01`)/total,
             het = (`01`*2)/total) %>%
      filter(total > 40) %>% 
      filter(mds_coord == mds) %>% nrow()
    if(callable_samples == 0){next}
    fst_genotyped %>%
      mutate(total =  (`00` + `01` + `11`)*2,
             perc_1 = ((`11`*2) + `01`)/total,
             het = (`01`*2)/total) %>%
      filter(total > 40) %>%
      filter(mds_coord == mds) %>%
      inner_join(.,labels) %>% 
      filter(species == "Ann") %>%
      ggplot(.,aes(x=perc_1,y=het)) +
      geom_point() +
      geom_segment(aes(x=0,xend=0.5,y=0,yend=1)) +
      geom_segment(aes(x=0.5,xend=1,y=1,yend=0)) +
      geom_segment(aes(x=0,xend=1,y=0,yend=0)) +
      theme_few() +
      ylab("Heterozygosity") +
      xlab("Proportion `1` allele") +
      facet_wrap(~is_wild)-> p1
    
    fst_genotyped %>%
      mutate(total =  (`00` + `01` + `11`)*2,
             perc_1 = ((`11`*2) + `01`)/total,
             het = (`01`*2)/total) %>%
      filter(total > 40) %>%
      filter(mds_coord == mds) %>%
      inner_join(.,labels) %>% 
      filter(species == "Ann") %>% 
      mutate(reference_count = perc_1[which(name == "SAM261")]) %>% 
      ggplot(.,aes(perc_1)) +
      geom_histogram() +
      theme_few() +
      ylab("Count") +
      xlab("Proportion `1` allele") +
      facet_wrap(~is_wild,scales="free_y") +
      geom_vline(data = fst_genotyped %>%
                   filter(name == "SAM261") %>%
                   mutate(total =  (`00` + `01` + `11`)*2,
                          perc_1 = ((`11`*2) + `01`)/total,
                          het = (`01`*2)/total) %>%
                   filter(total > 40) %>%
                   filter(mds_coord == mds) %>% head(),
                   aes(xintercept = perc_1),color="red") ->p2
  print(  
    grid.arrange(p1, p2, nrow = 2,top=mds)
  )
}
dev.off()
