
library(lostruct)
library(Matrix)
library(colorspace)
library(RColorBrewer)
library(ggmap)
library(scatterpie)
library(SNPRelate)
library(gridExtra)
library(ggExtra)

###Try comparing to structure
inversions = c("13a", "13b", "5a", "11a", "14a")

labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T)
samples <- read_tsv("/media/owens/Copper/wild_gwas_2018/annuus_sam/Annuus.tranche90.snp.all.90.bi.samplelist.txt",col_names = F)
colnames(samples) <- c("sample")


pdf("Annuus.tranche90.snp.all.90.bi.inversions.pdf")
all_calls <- tibble(name=character(),haplotype=numeric(),region=character())
for (inversion in inversions){
  struc = read_delim(paste("/media/owens/Copper/wild_gwas_2018/annuus_sam/Annuus.tranche90.snp.all.90.bi.",inversion,".2.meanQ",sep=""),col_names = F,
                     delim=" ")
  
  het = read_tsv(paste("/media/owens/Copper/wild_gwas_2018/annuus_sam/Annuus.tranche90.snp.all.90.bi.",inversion,".het.txt",sep=""))
  colnames(struc) <- c("k1","nothing","k2")
  struc %>% select(-nothing) -> struc
  cbind(struc, samples) %>% rename(name = sample)-> struc
  
  
  folder <- "/media/owens/Copper/wild_gwas_2018/annuus_sam/"
  gds.file <- paste("Annuus.tranche90.snp.all.90.bi.",inversion,".gds",sep="")
  vcf.file <- paste("Annuus.tranche90.snp.all.90.bi.",inversion,".vcf",sep="")
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
  
  print(
  ggplot(tab) + geom_point(aes(x=EV1,y=EV2,col=is_wild),size=3, alpha=0.5) +
    scale_color_brewer(palette = "Set1", name="Type") + theme_bw() +
    ylab(paste("PC",2," (",round(pc.percent[2],3)," PVE)",sep="")) +
    xlab(paste("PC",1," (",round(pc.percent[1],3)," PVE)",sep="")) +
    ggtitle(paste("Region",inversion))
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
   print(
   tab %>%
     inner_join(., struc) %>% 
     mutate(haplotype = case_when(k1 < 0.25 ~ 0,
                                  k1 > 0.75 ~ 2,
                                  TRUE ~ 1)) %>% 
     inner_join(., het) %>% 
     ggplot(.,aes(as.factor(haplotype),fill=is_wild)) +  geom_bar()  +
     theme_bw() + xlab("Haplotype") + ylab("Heterozygosity") +
     ggtitle(paste("Region",inversion)) +
     scale_color_brewer(palette = "Set1",name="Type")
   )
  
  tab %>%
    inner_join(., struc) %>% 
    mutate(haplotype = case_when(k1 < 0.25 ~ 0,
                                 k1 > 0.75 ~ 2,
                                 TRUE ~ 1)) %>% 
    select(name, haplotype)  -> calls
  calls$region <- inversion
  all_calls <- rbind(all_calls,calls)
  
}
dev.off()

pdf("Annuus.tranche90.snp.all.90.bi.inversions.genotyped.pdf")

all_calls %>% 
  inner_join(labels) %>%
  filter(is_wild == "landrace") %>%
  ggplot(.,aes(x=region,y=name)) + geom_tile(aes(fill=as.factor(haplotype))) +
  scale_fill_brewer(palette = "Set4",name="Haplotype") +
  ylab("Sample") + xlab("Region") +
  theme(axis.text.y = element_blank()) +ggtitle("Landraces")
all_calls %>% 
  inner_join(labels) %>%
  filter(is_wild == "cultivar") %>%
  ggplot(.,aes(x=region,y=name)) + geom_tile(aes(fill=as.factor(haplotype))) +
  scale_fill_brewer(palette = "Set4",name="Haplotype") +
  ylab("Sample") + xlab("Region") +
  theme(axis.text.y = element_blank()) +ggtitle("Cultivars")

dev.off()

all_calls %>% 
  inner_join(labels) %>% 
  write_tsv(.,"Annuus.tranche90.snp.all.90.bi.inversions.genotyped.txt")

all_calls %>% 
  inner_join(labels) %>%
  mutate(color = case_when(haplotype == 0 ~ brewer.pal(6,"Set1")[4],
                           haplotype == 1 ~ brewer.pal(6,"Set1")[5],
                           haplotype == 2 ~ brewer.pal(6,"Set1")[6])) %>%
  select(name,region,color) %>%
  write_tsv(., "Annuus.tranche90.snp.all.90.bi.inversions.phylocolors.txt")


###Plotting the Fst of inversions



  