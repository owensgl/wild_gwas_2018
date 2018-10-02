library(tidyverse)
library(SNPRelate)

folder <- "/media/owens/Copper/wild_gwas_2018/petiolaris/"
gds.file <- paste("Petiolaris.tranche90.snp.gwas.90.bi.gds",sep="")
vcf.file <- paste("Petiolaris.tranche90.snp.gwas.90.bi.vcf",sep="")
snpgdsVCF2GDS(paste(folder,vcf.file,sep="/"), paste(folder,gds.file,sep="/"), method="biallelic.only", ignore.chr.prefix = "HanXRQChr")

genofile <- snpgdsOpen(paste(folder,gds.file,sep="/"))
set.seed(1000)

snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,method="r", num.thread=10)
snpset.id <- unlist(snpset)

pca <- snpgdsPCA(genofile, num.thread=10, snp.id=snpset.id,
                 eigen.cnt = 40)
pc.percent <- pca$varprop*100
snpgdsClose(genofile)

labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T)
pop_loc <- read_tsv("pop_loc_allnum.txt")
pop_loc %>% rename(population = pop) %>% inner_join(.,labels) -> labels
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

ggplot(tab) + geom_point(aes(x=EV1,y=EV2,col=lat),size=3, alpha=0.5) +
  scale_color_viridis(name="Latitude") + theme_bw() +
  ylab(paste("PC",2," (",round(pc.percent[2],3)," PVE)",sep="")) +
  xlab(paste("PC",1," (",round(pc.percent[1],3)," PVE)",sep="")) +
  ggtitle(paste("Region",inversion)) +
  facet_wrap(~population)


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
