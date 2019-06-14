library(tidyverse)
#This is for taking the ML distance matrixes from inversion trees and then trying to see if its introgression or ancient inversion

directory <- "/media/owens/Copper/wild_gwas_2018/inv_gen_dist"
info1 <- read_tsv("sample_info_apr_2018.tsv") %>%
  rename(sample1 = name,
         species1 = species) %>% select(sample1, species1)
info2 <- read_tsv("sample_info_apr_2018.tsv") %>%
  rename(sample2 = name,
         species2 = species) %>% select(sample2, species2)

inversion_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt")

all_distances <- tibble(sample1=character(),sample2=character(),distance=numeric(),
                        inv_species=character(),species1=character(),species2=character())

for (i in 1:nrow(inversion_list)){
  
  
  chosen_species <- pull(inversion_list[i,1])
  chosen_inversion <-  paste("Ha412HOChr",sprintf("%02d",pull(inversion_list[i,4])),".", pull(inversion_list[i,6]),pull(inversion_list[i,5]),
                             sep="")
  if (!file_test("-f", paste(directory,"/",chosen_species,"/",chosen_inversion,"/",chosen_inversion,".mldist",sep=""))){
    next
  }
  dist_matrix <- read_delim(paste(directory,"/",chosen_species,"/",chosen_inversion,"/",chosen_inversion,".mldist",sep=""),
                            delim=" ",skip=1,col_names = F)
  
  colnames(dist_matrix) <- c("sample",pull(dist_matrix[,1]))
  
  dist_tidy <- gather(dist_matrix, sample1, distance, -sample) %>% 
    rename(sample2 = sample1, sample1 = sample) %>%
    mutate(distance = str_trim(distance)) %>%
    mutate(distance = as.numeric(distance)) %>%
    mutate(inv_species = chosen_species) %>%
    mutate(inversion = chosen_inversion) %>%
    inner_join(.,info1) %>% inner_join(.,info2) 
  all_distances <- rbind(all_distances, dist_tidy)
  
}
  

#For Annuus inversion testing
  

annuus_distances <- tibble(sample1=character(),sample2=character(),distance=numeric(),
                           inv_species=character(),species1=character(),species2=character(),
                           type=character())
all_distances %>%
  filter(inv_species != "annuus", inv_species != "argophyllus") %>%
  filter(species1 == "Ann",species2 == "Arg") %>%
  mutate(type = "baseline",inversion="baseline") -> ann_arg_baseline

ann_arg_baseline %>% pull(distance) %>% quantile(0.90) -> ann_arg_top
ann_arg_baseline %>% pull(distance) %>% quantile(0.10) -> ann_arg_bottom

for (i in 1:nrow(inversion_list)){
  chosen_species <- pull(inversion_list[i,1])
  chosen_inversion <-  paste("Ha412HOChr",sprintf("%02d",pull(inversion_list[i,4])),".", pull(inversion_list[i,6]),pull(inversion_list[i,5]),
                             sep="")
  if (chosen_species != "annuus"){
    next
  }
  geno <- read_tsv(paste("MDS_outliers/Ha412HO/annuus/Ha412HO_inv.jan09.pcasites.",chosen_inversion,".genotypes.txt",sep=""))  %>%
    rename(sample1 = "sample", type = triangle_genotype) %>%
    select(sample1, type)
  tmp <- all_distances %>%
    filter(inv_species == "annuus") %>%
    filter(inversion == chosen_inversion) %>%
    filter(species1 == "Ann",species2 == "Arg") %>%
    inner_join(.,geno) 
  annuus_distances <- rbind(annuus_distances,tmp)
}
annuus_distances %>%
  group_by(inversion, type) %>%
  summarize(mean_dist = mean(distance)) %>%
  group_by(inversion) %>%
  mutate(zero_score = mean_dist[[1]],two_score = mean_dist[[2]]) %>%
  mutate(order = case_when(zero_score > two_score ~ "Reverse",
         TRUE ~ "Normal")) %>% select(inversion, order) %>% unique() -> ordering
  
pdf("annuus.inv_gen_dist.v0.pdf",height=4,width=6)
annuus_distances %>%
  inner_join(.,ordering) %>%
  mutate(new_type = case_when(order == "Reverse" & type == 0 ~ 2,
                              order == "Reverse" & type == 2 ~ 0,
                              TRUE ~ type)) %>% 
  mutate(inversion = fct_relevel(inversion,
                                 "Ha412HOChr13.pos1","Ha412HOChr15.pos1",
                                 "Ha412HOChr16.pos2")) %>%
  ggplot(.,) + geom_boxplot(aes(x=inversion,y=distance,color=as.factor(new_type),fill=inversion)) +
  geom_hline(yintercept=arg_ann_top,linetype="dotted") + 
  geom_hline(yintercept=arg_ann_bottom,linetype="dotted") +
  #scale_fill_viridis_d(option = "magma") +
  scale_color_manual(values=rep("#000000",10)) +
  geom_rect(aes(ymin=arg_ann_bottom,ymax=arg_ann_top,xmin=0,xmax=length(unique(annuus_distances$inversion))+1
  ),alpha=0.1,fill="grey") +
  geom_boxplot(aes(x=inversion,y=distance,color=as.factor(new_type),fill=inversion)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme(legend.position="none") +
  ylab("Genetic distance") + xlab("Inversion") +
  coord_cartesian(xlim=c(1,length(unique(annuus_distances$inversion))))
dev.off()
#For Arg inversion testing


argophyllus_distances <- tibble(sample1=character(),sample2=character(),distance=numeric(),
                           inv_species=character(),species1=character(),species2=character(),
                           type=character())
all_distances %>%
  filter(inv_species != "annuus", inv_species != "argophyllus") %>%
  filter(species1 == "Arg",species2 == "Ann") %>%
  mutate(type = "baseline",inversion="baseline") -> arg_ann_baseline

arg_ann_baseline %>% pull(distance) %>% quantile(0.9) -> arg_ann_top
arg_ann_baseline %>% pull(distance) %>% quantile(0.1) -> arg_ann_bottom

for (i in 1:nrow(inversion_list)){
  chosen_species <- pull(inversion_list[i,1])
  chosen_inversion <-  paste("Ha412HOChr",sprintf("%02d",pull(inversion_list[i,4])),".", pull(inversion_list[i,6]),pull(inversion_list[i,5]),
                             sep="")
  if (chosen_inversion == "Ha412HOChr14.pos1"){
    next
  }
  if (chosen_species != "argophyllus"){
    next
  }
  geno <- read_tsv(paste("MDS_outliers/Ha412HO/argophyllus/Ha412HO_inv.jan09.pcasites.",chosen_inversion,".genotypes.txt",sep=""))  %>%
    rename(sample1 = "sample", type = triangle_genotype) %>%
    select(sample1, type)
  tmp <- all_distances %>%
    filter(inv_species == "argophyllus") %>%
    filter(inversion == chosen_inversion) %>%
    filter(species1 == "Arg",species2 == "Ann") %>%
    inner_join(.,geno) 
  argophyllus_distances <- rbind(argophyllus_distances,tmp)
}

argophyllus_distances %>%
  group_by(inversion, type) %>%
  summarize(mean_dist = mean(distance)) %>%
  group_by(inversion) %>%
  mutate(zero_score = mean_dist[[1]],two_score = mean_dist[[2]]) %>%
  mutate(order = case_when(zero_score > two_score ~ "Reverse",
                           TRUE ~ "Normal")) %>% select(inversion, order) %>% unique() -> ordering
pdf("argophyllus.inv_gen_dist.v0.pdf",height=4,width=6)
argophyllus_distances %>%
  inner_join(.,ordering) %>%
  mutate(new_type = case_when(order == "Reverse" & type == 0 ~ 2,
                              order == "Reverse" & type == 2 ~ 0,
                              TRUE ~ type)) %>% 
  mutate(inversion = fct_relevel(inversion,
                                 "Ha412HOChr10.neg1")) %>%
  ggplot(.,) + geom_boxplot(aes(x=inversion,y=distance,color=as.factor(new_type),fill=inversion)) +
  geom_hline(yintercept=arg_ann_top,linetype="dotted") + 
  geom_hline(yintercept=arg_ann_bottom,linetype="dotted") +
  #scale_fill_viridis_d(option = "magma") +
  scale_color_manual(values=rep("#000000",10)) +
  geom_rect(aes(ymin=arg_ann_bottom,ymax=arg_ann_top,xmin=0,xmax=length(unique(argophyllus_distances$inversion))+1
            ),alpha=0.1,fill="grey") +
  geom_boxplot(aes(x=inversion,y=distance,color=as.factor(new_type),fill=inversion)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme(legend.position="none") +
  ylab("Genetic distance") + xlab("Inversion") +
  coord_cartesian(xlim=c(1,length(unique(argophyllus_distances$inversion))))
dev.off()

#For ARG Chr6 inversion specifically.
pdf("argophyllus.inv_gen_dist.chr6.v0.pdf",width=6,height=4)
argophyllus_distances %>%
  inner_join(.,ordering) %>%
  mutate(new_type = case_when(order == "Reverse" & type == 0 ~ 2,
                              order == "Reverse" & type == 2 ~ 0,
                              TRUE ~ type)) %>% 
  filter(inversion == "Ha412HOChr06.pos1") %>% select(distance,new_type) %>%
  rbind(.,arg_ann_baseline %>% select(distance,type) %>% rename(new_type = type)) %>% 
  mutate(new_type = as.factor(new_type)) %>%
  mutate(new_type = fct_relevel(new_type, "baseline","0","2")) %>%
  group_by(new_type) %>%
  summarize(sd = sd(distance),mean=mean(distance)) %>%
  ggplot(.,) + geom_bar(aes(x=new_type,y=mean),stat="identity",fill="grey",color="black") +
  geom_errorbar(aes(x=new_type, ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  theme_bw() +
  coord_cartesian(ylim=c(0,0.025))
dev.off()




#For PetPet inversion testing


petpet_distances <- tibble(sample1=character(),sample2=character(),distance=numeric(),
                                inv_species=character(),species1=character(),species2=character(),
                                type=character())
all_distances %>%
  filter(inv_species != "petpet", inv_species != "petfal") %>%
  filter(species1 == "PetPet",species2 == "PetCan") %>%
  mutate(type = "baseline",inversion="baseline") -> petpet_petcan_baseline

petpet_petcan_baseline %>% pull(distance) %>% quantile(0.9) -> petpet_petcan_top
petpet_petcan_baseline %>% pull(distance) %>% quantile(0.1) -> petpet_petcan_bottom

for (i in 1:nrow(inversion_list)){
  chosen_species <- pull(inversion_list[i,1])
  chosen_inversion <-  paste("Ha412HOChr",sprintf("%02d",pull(inversion_list[i,4])),".", pull(inversion_list[i,6]),pull(inversion_list[i,5]),
                             sep="")
  if (chosen_species != "petpet"){
    next
  }
  geno <- read_tsv(paste("MDS_outliers/Ha412HO/petpet/Ha412HO_inv.jan09.pcasites.",chosen_inversion,".genotypes.txt",sep=""))  %>%
    rename(sample1 = "sample", type = triangle_genotype) %>%
    select(sample1, type)
  tmp <- all_distances %>%
    filter(inv_species == "petpet") %>%
    filter(inversion == chosen_inversion) %>%
    filter(species1 == "PetPet",species2 == "PetCan") %>%
    inner_join(.,geno) 
  petpet_distances <- rbind(petpet_distances,tmp)
}

petpet_distances %>%
  group_by(inversion, type) %>%
  summarize(mean_dist = mean(distance)) %>%
  group_by(inversion) %>%
  mutate(zero_score = mean_dist[[1]],two_score = mean_dist[[2]]) %>%
  mutate(order = case_when(zero_score > two_score ~ "Reverse",
                           TRUE ~ "Normal")) %>% select(inversion, order) %>% unique() -> ordering
petpet_distances %>%
  #rbind(.,arg_ann_baseline) %>%
  ggplot(.,aes(x=inversion,y=distance,fill=as.factor(type))) + geom_boxplot() +
  geom_hline(yintercept=petpet_petcan_top,linetype="dotted") + 
  geom_hline(yintercept=petpet_petcan_bottom,linetype="dotted") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1")

  
