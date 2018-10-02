pdf("Example_MDS_annuus.v0.pdf",height=4,width=18)
win.regions %>%
  gather(., mds, value, colnames(win.regions)[5:(ncol(win.regions)-1)]) %>% 
  filter(mds == "mds01") %>%
  mutate(chrom = gsub("HanXRQChr","",chrom)) %>%
  ggplot(aes(x=mid/1000000,y=value)) + geom_point(alpha=0.5) +
  facet_grid(~chrom,scales="free") + theme_few() +
  ylab("MDS coord 1") + xlab("MB")

win.regions %>%
  gather(., mds, value, colnames(win.regions)[5:(ncol(win.regions)-1)]) %>% 
  filter(mds == "mds05") %>%
  mutate(chrom = gsub("HanXRQChr","",chrom)) %>%
  ggplot(aes(x=mid/1000000,y=value)) + geom_point(alpha=0.5) +
  facet_grid(~chrom,scales="free") + theme_few() +
  ylab("MDS coord 5") + xlab("MB")

win.regions %>%
  gather(., mds, value, colnames(win.regions)[5:(ncol(win.regions)-1)]) %>% 
  filter(mds == "mds11") %>%
  mutate(chrom = gsub("HanXRQChr","",chrom)) %>%
  ggplot(aes(x=mid/1000000,y=value)) + geom_point(alpha=0.5) +
  facet_grid(~chrom,scales="free") + theme_few() +
  ylab("MDS coord 11") + xlab("MB")
dev.off()



####Example Fst plot

fst <- read_tsv("/media/owens/Copper/wild_gwas_2018/annuus/Annuus.tranche90.snp.env.90.bi.mds05_neg.fst.txt.gz")
w_size <- 1000000
pdf("Example_fst_annuus.v0.pdf",height=4,width=18)

fst %>%
  mutate(window = floor(pos/w_size)*w_size) %>%
  group_by(chr,window) %>%
  summarize(w_fst = sum(FstNum,na.rm=T)/sum(FstDenom,na.rm=T)) %>%
  ggplot(.,aes(x=window/1000000,y=w_fst)) +
  geom_line()+
  facet_wrap(~chr,scales = "free_x",nrow=1) +
  theme_few() +
  xlab("MB") +
  ylab(paste("Fst Window = ",w_size/1000000, " MB",sep="")) 
dev.off()


####All inversions at once.


annuus_stats <- read_tsv("/home/owens/bin/wild_gwas_2018/Annuus.tranche90.snp.env.90.bi.500.inversion_stats.txt")
annuus_windows <- read_tsv("/home/owens/bin/wild_gwas_2018/Annuus.tranche90.snp.env.90.bi.500.mds_cluster_windows.txt")

inner_join(annuus_windows, annuus_stats)  %>% mutate(species = "Annuus")-> annuus

petpet_stats <- read_tsv("/home/owens/bin/wild_gwas_2018/Petiolaris.tranche90.snp.petpet.90.bi.500.inversion_stats.txt")
petpet_windows <- read_tsv("/home/owens/bin/wild_gwas_2018/Petiolaris.tranche90.snp.petpet.90.bi.500.mds_cluster_windows.txt")

inner_join(petpet_windows, petpet_stats) %>% mutate(species = "PetiolarisPet")-> petpet

petfal_stats <- read_tsv("/home/owens/bin/wild_gwas_2018/Petiolaris.tranche90.snp.petfal.90.bi.500.inversion_stats.txt")
petfal_windows <- read_tsv("/home/owens/bin/wild_gwas_2018/Petiolaris.tranche90.snp.petfal.90.bi.500.mds_cluster_windows.txt")

inner_join(petfal_windows, petfal_stats) %>% mutate(species = "PetiolarisFal")-> petfal

argophyllus_stats <- read_tsv("/home/owens/bin/wild_gwas_2018/Argophyllus.tranche90.snp.gwas.90.bi.500.inversion_stats.txt")
argophyllus_windows <- read_tsv("/home/owens/bin/wild_gwas_2018/Argophyllus.tranche90.snp.gwas.90.bi.500.mds_cluster_windows.txt")

inner_join(argophyllus_windows, argophyllus_stats)  %>% mutate(species = "Argophyllus") -> argophyllus

anomalus_stats <- read_tsv("/home/owens/bin/wild_gwas_2018/Anomalus.tranche90.snp.gwas.90.bi.500.inversion_stats.txt")
anomalus_windows <- read_tsv("/home/owens/bin/wild_gwas_2018/Anomalus.tranche90.snp.gwas.90.bi.500.mds_cluster_windows.txt")

inner_join(anomalus_windows, anomalus_stats)  %>% mutate(species = "Anomalus") -> anomalus

all_species <- rbind(annuus,petpet, petfal,argophyllus,anomalus)

chrlengths <- read_tsv("/home/owens/ref/HanXRQr1.0-20151230.chrlength.txt")
chrlengths %>% rename(chrom = chr) -> chrlengths
pdf("Inversion_locations.v0.pdf")
all_species %>%
  filter(betweenSS > 0.90) %>%
  filter(!is.na(het_pvalue)) %>%
  filter(het_pvalue < 0.05) %>%
  filter(window_cluster < 10) %>%
  mutate(mds_species = paste(mds_coord, ":",species,sep="")) %>%
  #filter(chrom == "HanXRQChr01") %>%
  ggplot(.,aes(x=mid/1000000,y=mds_species)) + geom_point(aes(color=species)) +
  facet_wrap(~chrom,scales="free") + theme_few() +
  geom_segment(data=chrlengths, aes(x=0,xend=length/1000000,y=0,yend=0)) +
  xlab("MB") + scale_color_brewer(palette = "Set1") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

pdf("Inversion_counts.v0.pdf")
all_species %>%
  filter(betweenSS > 0.90) %>%
  filter(!is.na(het_pvalue)) %>%
  filter(het_pvalue < 0.05) %>%
  filter(window_cluster < 10) %>%
  group_by(species, mds_coord) %>%
  summarize(n = "tmp") %>%
  group_by(species) %>%
  summarize(count = n()) %>%
  ggplot(.,aes(x=species,y=count,fill=species)) + geom_bar(stat="Identity") +
  theme_few() + 
  scale_fill_brewer(palette = "Set1") +
  ggtitle("Credible Inversions")
dev.off()



