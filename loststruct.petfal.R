library(tidyverse)
library(lostruct)
library(Matrix)
library(colorspace)
library(RColorBrewer)
library(ggmap)
library(scatterpie)
library(SNPRelate)
library(gridExtra)
library(ggExtra)

bcf.file <- "/media/owens/Copper/wild_gwas_2018/petfal/Petiolaris.tranche90.snp.petfal.90.bi.bcf"
samples <- read_tsv("/media/owens/Copper/wild_gwas_2018/petfal/Petiolaris.tranche90.snp.petfal.90.bi.samplelist.txt",col_names = F)
colnames(samples) <- c("sample")
sites <- vcf_positions(bcf.file)
win.fn.snp <- vcf_windower(bcf.file, size=1000, type="snp", sites=sites) 
system.time( snp.pca <- eigen_windows(win.fn.snp,k=2) )
#Took 19677 time to run
system.time( pcdist <- pc_dist( snp.pca ) )

na.inds <- is.na( snp.pca[1,] )
mds <- cmdscale( pcdist, eig=TRUE, k=8 )
#mds.coords <- mds$points[ ifelse( na.inds, NA, cumsum(!na.inds) ), ]
mds.coords <- mds$points
colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))
win.regions <- region(win.fn.snp)()
win.regions %>% mutate(mid = (start + end) / 2) ->  win.regions


mincirc <- lostruct:::enclosing_circle( mds.coords[,1:2] )
mds.corners <- corners( mds.coords[,1:2], prop=.05 )
corner.cols <- c("red","blue","purple")
ccols <- rep("black",nrow(mds.coords))
for (k in 1:ncol(mds.corners)) {
  ccols[ mds.corners[,k] ] <- corner.cols[k]
}
win.regions$corner <-"0"
for (i in 1:3){
  win.regions$corner[mds.corners[,i]] <- i
}
win.regions$mds1 <- "NA"
win.regions$mds2 <- "NA"
win.regions$mds3 <- "NA"
win.regions$mds4 <- "NA"
win.regions$mds5 <- "NA"
win.regions$mds6 <- "NA"
win.regions$mds7 <- "NA"
win.regions$mds8 <- "NA"

for (i in 1:8){
  j = i + 5
  win.regions[,j] <- mds.coords[,i]
}
head(win.regions)

win.regions %>% ggplot(.,aes(x=mds1,y=mds2,color=corner)) + geom_point() + 
  scale_color_brewer(palette = "Set1") + theme_bw()

win.colors <- c("black",brewer.pal(3,"Set1"))
pdf("Petiolaris.tranche90.snp.petfal.90.bi.lostruct.1000.v0.pdf",height=16,width=25)
win.regions %>%
  gather(., mds, value, mds1, mds2, mds3,mds4, mds5, mds6, mds7, mds8) %>%
  ggplot(aes(x=mid,y=value,color=corner)) + geom_point(alpha=0.5) +
  facet_grid(mds~chrom,scales="free") + theme_bw() +
  scale_color_manual(values=win.colors,name="Corner")
dev.off()

all_pcs <- tibble(PC=character(),name=character(),value=numeric(),
                  chrom=character(),start=numeric(),end=numeric(),mid=numeric())
for (i in 1:nrow(win.regions)){
  tmp <- as.data.frame(snp.pca[i,4:length(snp.pca[i,])])
  tibble::rownames_to_column(tmp) -> tmp
  colnames(tmp) <- c("x","value")
  tmp %>% separate(x, into = c('PC', 'name'), sep = 5) -> tmp
  tmp$PC <- str_sub(tmp$PC, 1, 4)
  tmp$chrom <- win.regions$chrom[i]
  tmp$start <- win.regions$start[i]
  tmp$end <- win.regions$end[i]
  tmp$mid <- win.regions$mid[i]
  all_pcs <- rbind(all_pcs,tmp)
  print(paste("Loaded",i))
  
}
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T)
pop_loc <- read_tsv("pop_loc_allnum.txt")
pop_loc %>% rename(population = pop) %>% inner_join(.,labels) -> labels

all_pcs %>% inner_join(., labels) -> all_pcs

all_pcs %>% inner_join(., win.regions) -> all_pcs

saveRDS(all_pcs, file = "Petiolaris.tranche90.snp.petfal.90.bi.lostruct.1000.all_pcs.rds")



all_pcs %>% filter(mds3 < -0.2) %>% distinct(chrom, start, end)

all_pcs %>% filter(start == 107169083) %>% 
  spread(., PC, value) %>% 
  ggplot(.,aes(x=PC_1,y=PC_2,color=long)) + geom_point() +
  scale_color_viridis() + theme_bw()

pdf("Petiolaris.tranche90.snp.petfal.90.bi.lostruct.1000.closeup.v0.pdf")


all_pcs %>%
  filter(chrom == "HanXRQChr05") %>%
  select(chrom, start, end, mid, mds2, mds5, mds8) %>%
  distinct(chrom, start, end, mid, mds2, mds5, mds8) %>% 
  gather(., mds, value,  mds2, mds5, mds8) %>%
  ggplot(.,aes(x=mid/1000000,y=value)) + geom_point() +
  theme_bw() + xlab("MB") + ylab("MDS") +
  ggtitle("MDS 5b") +
  facet_wrap(~mds,nrow=3)


all_pcs %>%
  filter(chrom == "HanXRQChr09") %>%
  select(chrom, start, end, mid, mds1, mds2) %>%
  distinct(chrom, start, mid, end, mid, mds1, mds2) %>% 
  gather(., mds, value, mds1, mds2) %>%
  ggplot(.,aes(x=mid/1000000,y=value)) + geom_point() +
  theme_bw() + xlab("MB") + ylab("MDS") +
  ggtitle("MDS 9a") +
  facet_wrap(~mds,nrow=5)



all_pcs %>%
  filter(chrom == "HanXRQChr11") %>%
  select(chrom, start, end, mid, mds1, mds2, mds5 ) %>%
  distinct(chrom, start, end, mid, mds1, mds2, mds5) %>% 
  gather(., mds, value, mds1, mds2, mds5) %>%
  ggplot(.,aes(x=mid/1000000,y=value)) + geom_point() +
  theme_bw() + xlab("MB") + ylab("MDS") +
  ggtitle("MDS 11b") +
  facet_wrap(~mds,nrow=5)

all_pcs %>%
  filter(chrom == "HanXRQChr12") %>%
  select(chrom, start, end, mid,mds1, mds2 ) %>%
  distinct(chrom, start, end, mid,mds1, mds2) %>% 
  gather(., mds, value, mds1, mds2) %>%
  ggplot(.,aes(x=mid/1000000,y=value)) + geom_point() +
  theme_bw() + xlab("MB") + ylab("MDS") +
  ggtitle("MDS 12a") +
  facet_wrap(~mds,nrow=5)

all_pcs %>%
  filter(chrom == "HanXRQChr14") %>%
  select(chrom, start, end, mid, mds1, mds2,mds3, mds4, mds5 ) %>%
  distinct(chrom, start, end, mid, mds1, mds2,mds3, mds4, mds5) %>% 
  gather(., mds, value, mds1, mds2,mds3, mds4, mds5) %>% 
  ggplot(.,aes(x=mid/1000000,y=value)) + geom_point() +
  theme_bw() + xlab("MB") + ylab("MDS") +
  ggtitle("MDS 14b") +
  facet_wrap(~mds,nrow=5)



dev.off()
###Try comparing to structure
inversions = c("5b", "9a","11b","12a","14b")


pdf("Petiolaris.tranche90.snp.petfal.90.bi.inversions.pdf")
all_calls <- tibble(name=character(),haplotype=numeric(),region=character())
for (inversion in inversions){
  struc = read_delim(paste("/media/owens/Copper/wild_gwas_2018/petfal/Petiolaris.tranche90.snp.petfal.90.bi.",inversion,".2.meanQ",sep=""),col_names = F,
                     delim=" ")
  
  het = read_tsv(paste("/media/owens/Copper/wild_gwas_2018/petfal/Petiolaris.tranche90.snp.petfal.90.bi.",inversion,".het.txt",sep=""))
  colnames(struc) <- c("k1","nothing","k2")
  struc %>% select(-nothing) -> struc
  cbind(struc, samples) %>% rename(name = sample)-> struc
  
  
  folder <- "/media/owens/Copper/wild_gwas_2018/petfal/"
  gds.file <- paste("Petiolaris.tranche90.snp.petfal.90.bi.",inversion,".gds",sep="")
  vcf.file <- paste("Petiolaris.tranche90.snp.petfal.90.bi.",inversion,".vcf",sep="")
  snpgdsVCF2GDS(paste(folder,vcf.file,sep="/"), paste(folder,gds.file,sep="/"), method="biallelic.only", ignore.chr.prefix = "HanXRQChr")
  
  genofile <- snpgdsOpen(paste(folder,gds.file,sep="/"))
  set.seed(1000)
  
  pca <- snpgdsPCA(genofile, num.thread=10, 
                   eigen.cnt = 0)
  pc.percent <- pca$varprop*100
  snpgdsClose(genofile)
  
  tab <- data.frame(name = pca$sample.id,
                    EV1 = pca$eigenvect[,1],  
                    EV2 = pca$eigenvect[,2],    
                    stringsAsFactors = FALSE)
  inner_join(tab, labels) -> tab
  
  lat <- ggplot(tab) + geom_point(aes(x=EV1,y=EV2,col=lat),size=3, alpha=0.5) +
    scale_color_viridis(name="Latitude") + theme_bw() +
    ylab(paste("PC",2," (",round(pc.percent[2],3)," PVE)",sep="")) +
    xlab(paste("PC",1," (",round(pc.percent[1],3)," PVE)",sep="")) +
    ggtitle(paste("Region",inversion))
  long <- ggplot(tab) + geom_point(aes(x=EV1,y=EV2,col=long),size=3, alpha=0.5) +
    scale_color_viridis(name="Longitude") + theme_bw() +
    ylab(paste("PC",2," (",round(pc.percent[2],3)," PVE)",sep="")) +
    xlab(paste("PC",1," (",round(pc.percent[1],3)," PVE)",sep="")) +
    ggtitle(paste("Region",inversion))
  print(
    grid.arrange(lat, long, nrow= 2)
  )
  frame()
  tab %>%
    inner_join(., struc) %>% 
    ggplot(.,aes(x=EV1,y=k1)) + geom_point() +
    theme_bw() + xlab("PC1") + ylab("Structure assignment") +
    ggtitle(paste("Region",inversion)) -> p
  print(
    ggExtra::ggMarginal(p, type = "histogram")
  )
  print( 
    tab %>%
      inner_join(., struc) %>% 
      mutate(haplotype = case_when(k1 < 0.25 ~ 0,
                                   k1 > 0.75 ~ 2,
                                   TRUE ~ 1)) %>% 
      inner_join(., het) %>% 
      ggplot(.,aes(y=percent_het,x=haplotype,group=as.character(haplotype))) +  geom_boxplot()  +
      theme_bw() + xlab("Haplotype") + ylab("Heterozygosity") +
      ggtitle(paste("Region",inversion))
  )
  usa <- map_data('state')
  states <- map_data("state")
  texas <- subset(states, region %in% c("texas"))
  
  print(
    ggplot(usa, aes(long, lat)) +
      geom_map(map=usa, aes(map_id=region), fill=NA, color="black") +
      coord_quickmap() +
      geom_scatterpie(data=tab %>%
                        inner_join(., struc) %>% 
                        mutate(haplotype = case_when(k1 < 0.25 ~ 0,
                                                     k1 > 0.75 ~ 2,
                                                     TRUE ~ 1)) %>% 
                        group_by(population, lat, long, haplotype) %>% 
                        tally() %>%
                        spread(., haplotype, n,fill=0),
                      aes(x=long, y=lat, r=1.0), 
                      cols=c("0","1","2"), color=NA, alpha=.8) +
      scale_fill_brewer(name="Haplotype",palette = "Set1") +theme_bw() +
      xlab("Longitude") + ylab("Latitude") + 
      ggtitle(paste("Region",inversion)) +
      coord_map(xlim = c(-115,-96),ylim = c(30, 42))
  )
  
  tab %>%
    inner_join(., struc) %>% 
    mutate(haplotype = case_when(k1 < 0.25 ~ 0,
                                 k1 > 0.75 ~ 2,
                                 TRUE ~ 1)) %>% 
    select(name, haplotype)  -> calls
  calls$region <- inversion
  all_calls <- rbind(all_calls,calls)
  dune_info <- read_tsv("dune_sample_info.txt")
  print(
  tab %>%
    inner_join(.,dune_info) %>%
    inner_join(., struc) %>% 
    mutate(haplotype = case_when(k1 < 0.25 ~ 0,
                                 k1 > 0.75 ~ 2,
                                 TRUE ~ 1)) %>%
    ggplot(.,aes(ecotype)) + geom_bar(aes(fill=as.factor(haplotype))) + 
    facet_wrap(~spe) +
    theme_bw() + scale_fill_brewer(palette = "Set1",name="Haplotype") +
    ggtitle(paste("Region",inversion))
  )
  
}



dev.off()

pdf("Petiolaris.tranche90.snp.petfal.90.bi.inversions.genotyped.pdf")

all_calls %>% 
  inner_join(labels) %>%
  ggplot(.,aes(x=region,y=fct_reorder(name, long))) + geom_tile(aes(fill=as.factor(haplotype))) +
  scale_fill_brewer(palette = "Set4",name="Haplotype") +
  ylab("Sample") + xlab("Region") +
  theme(axis.text.y = element_blank()) +ggtitle("Arranged by Longitude")

dev.off()

all_calls %>% 
  inner_join(labels) %>%
  write_tsv(.,"Petiolaris.tranche90.snp.petfal.90.bi.inversions.genotyped.txt")

read_tsv("Petiolaris.tranche90.snp.petfal.90.bi.inversions.genotyped.txt") %>%
  mutate(color = case_when(haplotype == 0 ~ brewer.pal(6,"Set1")[4],
                           haplotype == 1 ~ brewer.pal(6,"Set1")[5],
                           haplotype == 2 ~ brewer.pal(6,"Set1")[6])) %>%
  select(name,region,color) %>%
  write_tsv(., "Petiolaris.tranche90.snp.petfal.90.bi.inversions.phylocolors.txt")

pdf("Petiolaris.tranche90.snp.petfal.90.bi.inversions.genotyped.linkage.pdf")
read_tsv("Petiolaris.tranche90.snp.petfal.90.bi.inversions.genotyped.txt") %>%
 filter(region == "9a" | region == "11b" | region == "14b") %>% 
  select(name, haplotype, region) %>%
  spread(., region, haplotype) %>%
  mutate(merged_haplotype = paste(`9a`,`11b`,sep="-")) %>%
  ggplot(.,aes(merged_haplotype)) + geom_bar() + 
  theme_bw() +
  ggtitle("9a + 11b")

read_tsv("Petiolaris.tranche90.snp.petfal.90.bi.inversions.genotyped.txt") %>%
  filter(region == "9a" | region == "11b" | region == "14b") %>% 
  select(name, haplotype, region) %>%
  spread(., region, haplotype) %>%
  mutate(merged_haplotype = paste(`9a`,`14b`,sep="-")) %>%
  ggplot(.,aes(merged_haplotype)) + geom_bar() + 
  theme_bw() +
  ggtitle("9a + 14b")

read_tsv("Petiolaris.tranche90.snp.petfal.90.bi.inversions.genotyped.txt") %>%
  filter(region == "9a" | region == "11b" | region == "14b") %>% 
  select(name, haplotype, region) %>%
  spread(., region, haplotype) %>%
  mutate(merged_haplotype = paste(`11b`,`14b`,sep="-")) %>%
  ggplot(.,aes(merged_haplotype)) + geom_bar() + 
  theme_bw() +
  ggtitle("11b + 14b")
dev.off()
####Plotting fst across inversions
inversions = c("5b", "9a","11b","12a","14b")

pdf("Petiolaris.tranche90.snp.petfal.90.bi.inversions.genotyped.fst.pdf",height=17,width=20)
for (inversion in inversions){
  
  fst <- read_tsv(paste("/media/owens/Copper/wild_gwas_2018/petfal/Petiolaris.tranche90.snp.petfal.90.bi.",inversion,".fst.txt",sep=""))
  
  window_size <- 10000
  print(
    fst %>%
      mutate(window = floor(pos/window_size) * window_size) %>%
      group_by(chr, window) %>%
      summarize(mean_fst = sum(FstNum)/sum(FstDenom)) %>% 
      ggplot(.,aes(x=window/1000000,y=mean_fst)) + geom_line() +
      coord_cartesian(ylim=c(0,1)) + facet_wrap(~chr,nrow=17) +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.text.y=element_blank()
      ) +
      ylab("Fst, 10KB window") + xlab("MB") +
      ggtitle(paste("Inversion", inversion))
  )
  window_size <- 100000
  print(
    fst %>%
      mutate(window = floor(pos/window_size) * window_size) %>%
      group_by(chr, window) %>%
      summarize(mean_fst = sum(FstNum)/sum(FstDenom)) %>% 
      ggplot(.,aes(x=window/1000000,y=mean_fst)) + geom_line() +
      coord_cartesian(ylim=c(0,1)) + facet_wrap(~chr,nrow=17) +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.text.y=element_blank()
      ) +
      ylab("Fst, 100KB window") + xlab("MB") +
      ggtitle(paste("Inversion", inversion))
  )
}
dev.off()
