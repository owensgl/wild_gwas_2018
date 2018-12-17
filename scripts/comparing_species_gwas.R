library(tidyverse)


labels <- read.table(info_name,header=T)

gene.gwas <- data.frame(Chr=as.numeric(),StartPos=as.numeric(),
                        EndPos=as.numeric(),Gene=as.character(),
                        TotalSites=as.numeric(),lowestp=as.numeric(),
                        secondp=as.numeric(),
                        phenotype=as.character(),
                        taxon=character())

for (spe in c("ANN","PET")){
  for (trait_choice in c("DTF","SLA")){
    for (i in 1:length(chrlist)){
      print(paste("Loading data for ",trait_choice," ",spe,sep=""))
      filename <- paste("/media/owens/Copper/wild_gwas_2018/prelim_gwas/",spe,"/",trait_choice,".HanXRQChr",chrlist[i],".lrt0.gene.txt",sep="")
      
      gene.tmp <- read.table(filename,header=T)
      gene.tmp$phenotype <- trait_choice
      gene.tmp$taxon <- spe
      gene.gwas <- rbind(gene.gwas,gene.tmp)
    }
  }
}
gene.gwas %>% 
  filter(!is.na(second_lrt)) %>%
  group_by(phenotype, taxon) %>%
  mutate(rank = row_number(desc(second_lrt)),
         percentile = rank(second_lrt)/length(second_lrt)) %>%
   ungroup()-> gene.gwas

gene.gwas %>% filter(Chr == i) %>%
  filter(phenotype == trait_choice, taxon == "ANN") %>%
  ggplot(.,aes(x=(StartPos/1000000),y=second_lrt)) + 
  geom_point(size=1) +
  geom_point(data=fulldata %>% filter(LRT > min_LRT, Chr == i) %>%
               filter(trait == trait_choice, taxon == "PET") %>%
               mutate(capped_LRT = ifelse(LRT > 100,100,LRT)),
             aes(x=(Position/1000000),y=-capped_LRT))


pdf("wild_gwas.annuus.MAF10.genewise.comparison.jan5.pdf")
for (trait_choice in c("DTF","SLA")){
  print(
  ggplot(merge(gene.gwas %>%
                 filter(phenotype == trait_choice, taxon == "ANN") %>%
                 mutate(ann_lrt = second_lrt),
               gene.gwas %>%
                 filter(phenotype == trait_choice, taxon == "PET") %>%
                 mutate(pet_lrt = second_lrt) %>%
                 select(StartPos,pet_lrt))
         ,aes(x=ann_lrt,y=pet_lrt),color="black") +
    geom_point() +
    geom_point(data=merge(gene.gwas %>%
                            filter(phenotype == trait_choice, taxon == "ANN") %>%
                            mutate(ann_lrt = second_lrt,ann_percentile = percentile),
                          gene.gwas %>%
                            filter(phenotype == trait_choice, taxon == "PET") %>%
                            mutate(pet_lrt = second_lrt,pet_percentile = percentile) %>%
                            select(StartPos,pet_lrt,pet_percentile)) %>%
                 mutate(highlight = ifelse((pet_percentile > 0.95 & ann_percentile > 0.95),1,0 )) %>%
                 filter(highlight == 1),color="red") +
    theme_bw() +
    xlab("H.annuus LRT") + 
    ylab("H.petiolaris LRT") +
    ggtitle(paste("Preliminary GWAS for ",trait_choice,"\nHightlights overlapping top 5%",sep=""))
  )
}
dev.off()

for (trait_choice in c("DTF","SLA")){
  write.table(
merge(gene.gwas %>%
        filter(phenotype == trait_choice, taxon == "ANN") %>%
        mutate(ann_lrt = second_lrt,ann_percentile = percentile),
      gene.gwas %>%
        filter(phenotype == trait_choice, taxon == "PET") %>%
        mutate(pet_lrt = second_lrt,pet_percentile = percentile) %>%
        select(StartPos,pet_lrt,pet_percentile)) %>%
  mutate(highlight = ifelse((pet_percentile > 0.95 & ann_percentile > 0.95),1,0 )) %>%
  filter(highlight == 1) %>% select("Chr","StartPos","EndPos","Gene",
                                    "phenotype","pet_lrt","ann_lrt",
                                    "pet_percentile","ann_percentile"),
  file=paste("wild_gwas.annuus.MAF10.genewise.",trait_choice,".tophits.jan5.txt"),
  row.names=F,quote=F
  )
}
  



merge(gene.gwas %>%
        filter(phenotype == trait_choice, taxon == "ANN") %>%
        mutate(ann_lrt = second_lrt,ann_percentile = percentile),
      gene.gwas %>%
        filter(phenotype == trait_choice, taxon == "PET") %>%
        mutate(pet_lrt = second_lrt,pet_percentile = percentile) %>%
        select(StartPos,pet_lrt,pet_percentile)) %>%
  mutate(pet_outlier = ifelse((pet_percentile > 0.95),1,0 ),
         ann_outlier = ifelse((ann_percentile > 0.95),1,0)) %>%
  group_by(pet_outlier, ann_outlier) %>%
  summarize(counts=n()) %>% ungroup() %>%select(counts) %>% .$counts -> counts
  chisq.test(matrix(counts,nrow=2))
