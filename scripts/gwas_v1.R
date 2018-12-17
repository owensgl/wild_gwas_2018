library(qvalue)
library(tidyverse)
library(qqman)


chrlist <- seq(17) %>% formatC(., width = 2, format = "d", flag = "0")
fulldata <- data.frame(species=factor(),trait=factor(),Chromosome=factor(),Position=numeric(),Frequency=numeric(),
                       N=numeric(),LRT=numeric(),pvalue=numeric(),taxon=factor(),
                       trait=factor())
for (spe in c("ANN","PET")){
  for (trait in c("DTF","SLA")){
    for (i in 1:length(chrlist)){
      filename <- paste("/media/owens/Copper/wild_gwas_2018/prelim_gwas/",spe,"/",trait,".HanXRQChr",chrlist[i],".lrt0.gz",sep="")
      data <- read.table(gzfile(filename),
                         colClasses=c("factor","numeric","NULL","NULL","numeric","numeric","numeric","NULL"),header=T)
      data %>% filter(Frequency >= 0.1, Frequency <= 0.90) %>% filter(LRT >= 0) -> data
      data %>% mutate(pvalue = 1- pchisq(LRT,1)) -> data
      data$taxon <- spe
      data$trait <- trait
      fulldata <- rbind(fulldata, data)
      print(paste(spe,"-",trait, " Chromosome ",chrlist[i]," loaded!",sep = ""))
    }
  }
}


q<- qvalue(fulldata$pvalue, fdr.level=0.05)
fulldata$qvalue.fdr05 <- q$qvalues
max_LRT <- max(fulldata$LRT)
#max_LRT <- 85
min_LRT <- 5
fulldata$Chr <- as.numeric(substr(fulldata$Chromosome, 10, 11))
pdf("wild_gwas.annuus.SLA.MAF10.noPCA.pdf",width=9,height=5)

pdf("wild_gwas.annuus.MAF10.QQ.jan5.pdf",width=9,height=5)
for (spe in c("ANN","PET")){
  for (trait_choice in c("DTF","SLA")){
    print(
      qq(fulldata %>% filter(trait == trait_choice, taxon == spe) %>%
        sample_n(1000000) %>% .$pvalue, main = paste(spe,"-",trait_choice, " Q-Q plot of GWAS p-values"))
    )
  }
}
dev.off()

pdf("wild_gwas.annuus.MAF10.ANN.jan5.pdf",width=9,height=5)
for (spe in c("ANN")){
  for (trait_choice in c("DTF","SLA")){
    for (i in 1:17){
      print(
        fulldata %>% filter(LRT > min_LRT, Chr == i) %>%
          filter(trait == trait_choice, taxon == spe) %>%
          ggplot(.,aes(x=(Position/1000000),y=LRT)) + 
          geom_point(size=1) +
          theme_minimal() + 
          ggtitle(paste("Chromosome",chrlist[i],"\nPreliminary ",spe,"-",trait_choice," GWAS, MAF>0.1",sep="")) +
          xlab("MBase") + ylab("Likelihood Ratio Test") 
      )
    }
  }
}
dev.off()


genes <- read.table("/home/owens/ref/HanXRQr1.0-20151230-EGN-r1.1.gff3",header=F,na.strings="#",sep="\t",quote="")
colnames(genes) <- c("chr","program","type","start","end","score","strand","phase","info")

pdf("wild_gwas.annuus.MAF10.ANN.SLA.HanXRQChr06g0179621.jan5.pdf")
focal_chr <- "HanXRQChr06"
regions <- c(57000000, 57500000)
genelist <- genes %>% filter(chr ==  focal_chr) %>% filter(start > regions[1]) %>% filter( end < regions[2])
p <- fulldata %>% filter(Chromosome == focal_chr, Position > regions[1], Position < regions[2]) %>%
  filter(taxon == "ANN", trait == "SLA") %>%
  mutate(max_LRT = max(LRT)) %>%
  ggplot(data=.) + geom_point(size=2,alpha=0.2,aes(x=Position, y=LRT)) + 
  coord_cartesian(ylim = c(-5,60)) +
  ggtitle(paste("HanXRQChr06g0179621\n","ANN-SLA",sep="")) +
  geom_rect(data=genes %>% filter(chr ==  focal_chr) %>% filter(start > regions[1]) %>% filter( end < regions[2]),
            aes(xmin=start,xmax=end,ymin=-5,ymax=-1),fill="dark green") +
  theme_bw()
  
p
dev.off()




#######Plot top hits

fulldata %>% filter(LRT > 22) -> tophits

genes <- read.table("/home/owens/ref/HanXRQr1.0-20151230-EGN-r1.1.gff3",header=F,na.strings="#",sep="\t",quote="")
colnames(genes) <- c("chr","program","type","start","end","score","strand","phase","info")

candidate_genes <- data.frame(chr= as.character(),program=as.character(),type=as.character(),start=as.character(),
                              end=as.character(),strand=as.character(),phase=as.character(),info=as.character())
pdf("")
windows <- data.frame(chr = as.character(),start=as.numeric(),end=as.numeric())
for (i in 1:nrow(tophits)){
  focal_chr <- as.character(tophits$Chromosome[i])
  focal<- tophits$Position[i]
  surround <- 1000000
  
  in_window <- windows %>% filter(chr ==  focal_chr, start < focal,end>focal)
  if (nrow(in_window) > 0){
    print(paste("Skipping",focal_chr,focal,"because it's already been printed",sep=" "))
    next
  }
  current_window <- data.frame(chr = focal_chr,start=(focal-surround), end=(focal+surround))
  windows <- rbind(windows,current_window)
  print(paste("Printing", focal_chr,current_window$start, "to",current_window$end,sep=" "))
  genelist <- genes %>% filter(chr ==  focal_chr) %>% filter(start > (focal-surround)) %>% filter( end < (focal+surround))
  if (nrow(genelist) > 0){
    candidate_genes <- rbind(candidate_genes, genelist)
    p <- fulldata %>% filter(Chromosome == focal_chr, Position > (focal-surround), Position < (focal+surround)) %>%
      ggplot(data=.) + geom_point(size=2,aes(x=Position, y=LRT)) + 
      ggtitle(paste(focal_chr," ",focal)) +
      geom_rect(data=genes %>% filter(chr ==  focal_chr) %>% filter(start > (focal-surround)) %>% filter( end < (focal+surround)),
                aes(xmin=start,xmax=end,ymin=-5,ymax=-1),fill="dark green") +
      coord_cartesian(ylim = c(-5, max_LRT), xlim = c((focal-surround), (focal+surround))) +
      geom_line(data=hit_windows %>% filter(Chromosome == focal_chr),aes(x=w_pos1+5000,y=(10*percent_hits)),color="blue",size=1)
    
    
    
    print(
      
      p
      
    )
  }else{
    p<- fulldata %>% filter(Chromosome == focal_chr, Position > (focal-surround), Position < (focal+surround)) %>%
      ggplot(data=.) + geom_point(size=2,aes(x=Position, y=LRT)) + scale_colour_gradient(low = "red", high = "black") +
      ggtitle(paste(focal_chr," ",focal)) +
      coord_cartesian(ylim = c(-5, max_LRT),xlim = c((focal-surround), (focal+surround))) +
      geom_line(data=hit_windows %>% filter(Chromosome == focal_chr),aes(x=w_pos1+5000,y=(10*percent_hits)),color="blue",size=1)
    
    
    print(
      
      p
      
    )
  }
}
dev.off()

candidate_genes <- candidate_genes %>% filter(type == "mRNA") %>% unique(.)
write.table(candidate_genes,"annuus_herbicide.asso.nov1.MAF10.tophits.genelists.txt",
            quote=F,col.names=T,row.names=F)

##PRINT TOP REGIONS
focal_chr <- "HanXRQChr13"
regions <- c(21500000, 23500000)
genelist <- genes %>% filter(chr ==  focal_chr) %>% filter(start > regions[1]) %>% filter( end < regions[2])
p <- fulldata %>% filter(Chromosome == focal_chr, Position > regions[1], Position < regions[2]) %>%
  ggplot(data=.) + geom_point(size=2,aes(x=Position, y=LRT)) + 
  ggtitle(paste(focal_chr,":",regions[1],"to",regions[2])) +
  geom_rect(data=genes %>% filter(chr ==  focal_chr) %>% filter(start > regions[1]) %>% filter( end < regions[2]),
            aes(xmin=start,xmax=end,ymin=-5,ymax=-1),fill="dark green") +
  coord_cartesian(ylim = c(-5, max_LRT))
p
pdf("")
p
dev.off()

focal_chr <- "HanXRQChr01"
regions <- c(85459631, 87459631)
genelist <- genes %>% filter(chr ==  focal_chr) %>% filter(start > regions[1]) %>% filter( end < regions[2])
p <- fulldata %>% filter(Chromosome == focal_chr, Position > regions[1], Position < regions[2]) %>%
  ggplot(data=.) + geom_point(size=2,aes(x=Position, y=LRT)) + 
  ggtitle(paste(focal_chr,":",regions[1],"to",regions[2])) +
  geom_rect(data=genes %>% filter(chr ==  focal_chr) %>% filter(start > regions[1]) %>% filter( end < regions[2]),
            aes(xmin=start,xmax=end,ymin=-5,ymax=-1),fill="dark green") +
  coord_cartesian(ylim = c(-5, max_LRT))
p
pdf("annuus_herbicide.asso.nov1.MAF10.chr01.85459631.87459631.pdf")
p
dev.off()

focal_chr <- "HanXRQChr04"
regions <- c(155107025, 157107025)
genelist <- genes %>% filter(chr ==  focal_chr) %>% filter(start > regions[1]) %>% filter( end < regions[2])
p <- fulldata %>% filter(Chromosome == focal_chr, Position > regions[1], Position < regions[2]) %>%
  ggplot(data=.) + geom_point(size=2,aes(x=Position, y=LRT)) + 
  ggtitle(paste(focal_chr,":",regions[1],"to",regions[2])) +
  geom_rect(data=genes %>% filter(chr ==  focal_chr) %>% filter(start > regions[1]) %>% filter( end < regions[2]),
            aes(xmin=start,xmax=end,ymin=-5,ymax=-1),fill="dark green") +
  coord_cartesian(ylim = c(-5, max_LRT))
p
pdf("annuus_herbicide.asso.nov1.MAF10.chr01.85459631.87459631.pdf")
p
dev.off()

#######Try with higher MAF filter
fulldata %>% filter(Frequency >= 0.2) -> fulldata.MAF20
q<- qvalue(fulldata.MAF20$pvalue, fdr.level=0.1)
fulldata.MAF20$qvalue.fdr10 <- q$qvalues

fulldata %>% filter(LRT < 0) %>% nrow()

q<- qvalue(fulldata$pvalue)


####Trying to make a sliding window

fulldata$w_pos2 <- (((fulldata$Position / 10000) %>% floor) + 1)*10000
fulldata$w_pos1 <- fulldata$w_pos2 - 9999

fulldata %>% mutate(hit = ifelse(LRT >=10,1,0)) -> fulldata
percent_hits <- sum(fulldata$hit)/nrow(fulldata)
fulldata %>% ungroup %>%
  group_by(Chromosome, Chr, w_pos1,w_pos2) %>%
  summarise(percent_hits = mean(hit),n_sites = n(), total_hits = sum(hit)) %>% 
  filter(n_sites > 20) -> hit_windows

ggplot(hit_windows,aes(x=w_pos1,y=percent_hits)) + geom_line() + facet_grid(Chr~.)

hit_windows %>% group_by(Chromosome, Chr, w_pos1,w_pos2) %>%
  do(binomial = binom.test(total_hits,n_sites,percent_hits,alternative="greater"))

bifunc  <- function(x,n){
  bi <- binom.test(x,n,percent_hits,alternative="greater")$p.value
  data_frame(pvalue = bi)}
hit_windows %>%
  group_by(Chromosome, Chr, w_pos1,w_pos2) %>%
  do(bifunc(.$total_hits,.$n_sites)) -> window_pvalues

q<- qvalue(window_pvalues$pvalue)
window_pvalues$qvalue.fdr10 <- q$qvalues

ggplot(window_pvalues,aes(x=w_pos1,y=-log10(qvalue.fdr10))) + geom_line() + facet_grid(Chr~.)
