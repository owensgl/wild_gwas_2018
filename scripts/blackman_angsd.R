library(tidyverse)
library(ggthemes)
library(ggExtra)
library(ggpubr)

directory <- "/media/owens/Copper/wild_gwas/annuus"
info <- read_tsv("Angsd_most_aDNA.sampleinfo.txt") %>%
  rename(sample = Beagle_ID)
all_calls <- tibble(sample=character(),mds=character(),g00=numeric(),g01=numeric(),g11=numeric(),
                    Path=character(),Other_name=character(),Type=character(),Use=character(),
                    Lat=numeric(),Long=numeric(),total=numeric(),het=numeric(),
                    percent_1=numeric(),left_formula=numeric(),right_formula=numeric(),
                    sample_call=numeric())
mds_list <- read_tsv(paste(directory,"/Jan09.invlist.txt",sep=""),col_names = c("mds"))
bar_order <- c("Wild", "Archaeological", "Ethnographic","Landrace","Cultivar")
for (n in 1:nrow(mds_list)){
  mds_chosen <- pull(mds_list[n,1])
  
  data <- read_tsv(paste(directory, "/Ha412HO_inv.jan09.annuus.",mds_chosen,".blackman.mostaDNA.txt.gz",sep=""))

  #pdf(paste("Ha412HO_inv.jan09.annuus.",mds_chosen,".blackman.mostaDNA.pdf",sep=""),height=6,width=9)
  print(
  data %>%
    inner_join(.,info) %>%
    #filter(mds == mds_chosen) %>%
    mutate(total = (g00 + g01 + g11)*2,
           het= g01 / (g00 + g11 + g01),
           percent_1 = (g01 + (g11*2))/total) %>%
    ggplot(.,aes(x=percent_1)) +
    geom_histogram(aes(fill=Type)) +
    theme_few() + 
    facet_grid(Type~mds,scales="free_y") +
    scale_fill_brewer(palette = "Set1") +
    ylab("Count") + xlab("Percent inversion 1 alleles") +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  )
  
  data %>%
    inner_join(.,info) %>%
    #filter(mds == mds_chosen) %>%
    mutate(total = (g00 + g01 + g11)*2,
           het= g01 / (g00 + g11 + g01),
           percent_1 = (g01 + (g11*2))/total,
           left_formula = ((-(2/3)*percent_1) + (2/3)),
           right_formula = (((2/3)*percent_1))) %>%
    mutate(sample_call = case_when(percent_1 < 0.5  & left_formula >= het ~ 0,
                                   percent_1 < 0.5  & left_formula < het ~ 1,
                                   percent_1 >= 0.5  & right_formula >= het ~ 2,
                                   percent_1 >= 0.5  & right_formula < het ~ 1)) -> data_calls
  
  print(
  data_calls %>%
   # filter(total > 10) %>%
    ggplot(.,aes(x=percent_1,y=het)) +
    geom_segment(aes(x=0,xend=0.5,y=0,yend=1),color="grey") +
    geom_segment(aes(x=0.5,xend=1,y=1,yend=0),color="grey") +
    geom_segment(aes(x=0,xend=1,y=0,yend=0),color="grey") +
    geom_point(alpha=0.4,aes(color=as.factor(sample_call)),size=4) +
    geom_segment(aes(x=0.25,xend=0.5,y=0.5,yend=(1/3)),linetype="dotted") +
    geom_segment(aes(x=0.5,xend=0.75,y=(1/3),yend=0.5),linetype="dotted") +
    geom_segment(aes(x=0.5,xend=0.5,y=0,yend=(1/3)),linetype="dotted") +
    scale_color_brewer(palette = "Set1",name="Genotype") +
    theme_few() +
    ylab("Heterozygosity") +
    xlab("Proportion `1` allele") +
    ggtitle(paste(mds_chosen)) +
    facet_wrap(~Type)
  )
    
  #For less heterozygosity
  data %>%
    inner_join(.,info) %>%
    #filter(mds == mds_chosen) %>%
    mutate(total = (g00 + g01 + g11)*2,
           het= g01 / (g00 + g11 + g01),
           percent_1 = (g01 + (g11*2))/total,
           left_formula = ((-(1/3)*percent_1) + (1/3)),
           right_formula = (((1/3)*percent_1))) %>%
    mutate(sample_call = case_when(percent_1 < 0.5  & left_formula >= het ~ 0,
                                   percent_1 < 0.5  & left_formula < het ~ 1,
                                   percent_1 >= 0.5  & right_formula >= het ~ 2,
                                   percent_1 >= 0.5  & right_formula < het ~ 1)) -> data_calls

  print(
  data_calls %>%
   # filter(total > 10) %>%
    ggplot(.,aes(x=percent_1,y=het)) +
    geom_segment(aes(x=0,xend=0.5,y=0,yend=1),color="grey") +
    geom_segment(aes(x=0.5,xend=1,y=1,yend=0),color="grey") +
    geom_segment(aes(x=0,xend=1,y=0,yend=0),color="grey") +
    geom_point(alpha=0.4,aes(color=as.factor(sample_call)),size=4) +
    geom_segment(aes(x=0.12,xend=0.5,y=0.25,yend=(1/6)),linetype="dotted") +
    geom_segment(aes(x=0.5,xend=0.88,y=(1/6),yend=0.25),linetype="dotted") +
    geom_segment(aes(x=0.5,xend=0.5,y=0,yend=(1/6)),linetype="dotted") +
    scale_color_brewer(palette = "Set1",name="Genotype") +
    theme_few() +
    ylab("Heterozygosity") +
    xlab("Proportion `1` allele") +
    ggtitle(paste(mds_chosen, "LOWER HET THRESHOLD")) +
    facet_wrap(~Type)
  )
  print(
  data_calls %>%
    #filter(total > 10) %>%
    ggplot(.,aes(x=Type,fill=as.factor(sample_call),group=as.factor(sample_call))) + 
    geom_bar(position="fill",stat="count") +
    scale_fill_brewer(palette = "Set1",name="Genotype") +
    theme_few() +
    ggtitle(paste(mds_chosen, "LOWER HET THRESHOLD")) +
    scale_x_discrete(limits = bar_order)
  
  )

  #dev.off()
  all_calls <- rbind(all_calls,data_calls)
}

pdf(paste("Ha412HO_inv.jan09.annuus.blackman.mostaDNA.pdf",sep=""),height=6,width=9)
all_calls %>%
  #filter(total > 10) %>%
  ggplot(.,aes(x=Type,fill=as.factor(sample_call),group=as.factor(sample_call))) + 
  geom_bar(position="fill",stat="count") +
  scale_fill_brewer(palette = "Set1",name="Genotype") +
  theme_few() +
  facet_wrap(~mds) +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  scale_x_discrete(limits = bar_order)
dev.off()

pdf(paste("Ha412HO_inv.jan09.annuus.blackman.mostaDNA.june2019filtered.pdf",sep=""),height=6,width=9)
chosen_mds = c( "Ha412HOChr05.pos1", "Ha412HOChr11.pos1", "Ha412HOChr13.pos1", "Ha412HOChr14.neg1", "Ha412HOChr15.pos1", "Ha412HOChr16.neg1")
all_calls %>%
  filter(mds %in% chosen_mds ) %>%
  mutate(new_call = case_when(sample_call == 0 & total < 3 ~ "REMOVED",
                                 sample_call == 2 & total < 3 ~ "REMOVED",
                                 sample_call == 1 & total < 6 ~ "REMOVED",
                                 TRUE ~ as.character(sample_call))) %>%  View()
  ggplot(.,aes(x=Type,fill=as.factor(new_call),group=as.factor(new_call))) + 
  geom_bar(position="fill",stat="count") +
  scale_fill_brewer(palette = "Set1",name="Genotype") +
  theme_few() +
  facet_wrap(~mds) +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  scale_x_discrete(limits = bar_order)
all_calls %>%
  filter(mds %in% chosen_mds ) %>%
  mutate(new_call = case_when(sample_call == 0 & total < 3 ~ "REMOVED",
                              sample_call == 2 & total < 3 ~ "REMOVED",
                              sample_call == 1 & total < 6 ~ "REMOVED",
                              TRUE ~ as.character(sample_call))) %>%  
  mutate(new_call = na_if(new_call, "REMOVED")) %>%
  ggplot(.,aes(x=Type,fill=as.factor(new_call),group=as.factor(new_call))) + 
  geom_bar(position="fill",stat="count") +
  scale_fill_brewer(palette = "Set1",name="Genotype") +
  theme_few() +
  facet_wrap(~mds) +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  scale_x_discrete(limits = bar_order)
all_calls %>%
  filter(mds %in% chosen_mds ) %>%
  mutate(new_call = case_when(sample_call == 0 & total < 3 ~ "REMOVED",
                              sample_call == 2 & total < 3 ~ "REMOVED",
                              sample_call == 1 & total < 6 ~ "REMOVED",
                              TRUE ~ as.character(sample_call))) %>%  
  mutate(new_call = na_if(new_call, "REMOVED")) %>%
  filter(!is.na(new_call)) %>%
  ggplot(.,aes(x=Type,fill=as.factor(new_call),group=as.factor(new_call))) + 
  geom_bar(position="fill",stat="count") +
  scale_fill_brewer(palette = "Set1",name="Genotype") +
  theme_few() +
  facet_wrap(~mds) +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  scale_x_discrete(limits = bar_order)
dev.off()



pdf(paste("Ha412HO_inv.jan09.annuus.blackman.mostaDNA.markercounts.pdf",sep=""),height=6,width=9)
all_calls %>%
  #filter(total > 10) %>%
  ggplot(.,aes(x=Type,y=total)) + geom_boxplot() +
  facet_wrap(~mds,scales = "free_y") +
  scale_x_discrete(limits = bar_order) +
  theme_few() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  ggtitle("Number of markers used") 

all_calls %>%
 # filter(total > 10) %>%
  filter(Type == "Archaeological") %>%
  ggplot(.,aes(x=mds,y=total)) + geom_boxplot() +
  #facet_wrap(~mds,scales = "free") +
  theme_few() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  ggtitle("Number of markers used Archaeological") + geom_jitter()
dev.off()
  
pdf(paste("Ha412HO_inv.jan09.annuus.blackman.mostaDNA.markerplots.pdf",sep=""),height=12,width=12)

#Plotting site by site
for (n in 1:nrow(mds_list)){
  mds_chosen <- pull(mds_list[n,1])
  
  chosen_types <- c("Wild","Archaeological","Ethnographic","Landrace","Cultivar")
  data <- read_tsv(paste(directory, "/Ha412HO_inv.jan09.annuus.",mds_chosen,".blackman.mostaDNA.sites.txt.gz",sep=""))
  for (x in 1:length(chosen_types)){
    chosen_type <- chosen_types[x]
    data %>%
      inner_join(info) %>%
      filter(Type == chosen_type) %>%
      select(pos) %>% unique() %>% pull() -> positions
    n_labels <- ceiling(length(positions)/40)
    x_labels <- sort(positions)[seq(0, length(positions), by= n_labels)]
    
    print(
    data %>%
      inner_join(info) %>%
      filter(Type == chosen_type) %>%
      inner_join(all_calls) %>%
      ggplot(aes(x=as.factor(pos),y=fct_reorder(sample, percent_1))) + geom_tile(aes(fill=as.factor(genotype))) + 
      scale_fill_brewer(palette = "Set1",name="Genotype") +
      theme(axis.text.x=element_text(angle=60, hjust=1)) +
      theme_bw() +
      scale_x_discrete(breaks=x_labels) +
      theme(axis.text.x=element_text(angle=60, hjust=1)) +
      ylab("Sample") + xlab("Position") +
      ggtitle(paste(chosen_type,mds_chosen))
    )
  }

}
dev.off()
  
write_tsv(all_calls %>% select(-left_formula,-right_formula)
          ,paste("Ha412HO_inv.jan09.annuus.blackman.mostaDNA.summary.txt",sep=""))

pdf(paste("Ha412HO_inv.jan09.annuus.blackman.mostaDNA.depthtomarkers.pdf",sep=""),height=12,width=12)

all_calls %>% 
  mutate(sites_genotyped = total/2) %>%
  filter(!is.na(Nuclear_coverage)) %>%
  mutate(n_sites_gt_5 = case_when(sites_genotyped >= 5 ~"TRUE",
                                  TRUE ~ "FALSE")) %>%
  ggplot(.,aes(x=Nuclear_coverage,y=sites_genotyped)) +
  geom_point(aes(color=n_sites_gt_5)) +
  facet_wrap(~mds,scales="free_y") +
  theme_bw() +
  scale_color_brewer(palette = "Set1")
dev.off()


