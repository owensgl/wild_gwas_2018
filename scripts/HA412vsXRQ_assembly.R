#Comparing reference genomes
library(tidyverse)

pdf("Compare_references_v0.pdf",height=6,width=10)
chr <- list.files("/home/owens/working/reference_comparison/",pattern="allcomp.rev.copmem.txt")
chr <- gsub(".allcomp.rev.copmem.txt","",chr)
for (chr_chosen in chr){
  data.for <- read_tsv(paste("/home/owens/working/reference_comparison/",chr_chosen,".for.copmem.txt",sep=""),comment=">",
                     col_names=c("chr", "XRQ_pos","HA412_pos","length"))
  data.for$direction <- "forward"
  data.rev <- read_tsv(paste("/home/owens/working/reference_comparison/",chr_chosen,".rev.copmem.txt",sep=""),comment=">",
                       col_names=c("chr", "XRQ_pos","HA412_pos","length"))
  data.rev$direction <- "reverse"
  
  data <- rbind(data.for,data.rev)
  print(
  data %>% filter(length > 900) %>%
    ggplot(.,aes(x=XRQ_pos,y=HA412_pos)) + geom_point(aes(color=direction)) +
    facet_wrap(~direction) +
    ggtitle(paste(chr_chosen, "\ncopmem, 900bp match",sep="")) +
    scale_color_brewer(palette = "Set1") + 
    theme_bw()
  )
}
copmem.data %>% filter(length > 300) %>%
  ggplot(.,aes(x=XRQ_pos,y=HA412_pos)) + geom_point(aes(color=direction)) +
  facet_wrap(~direction) +
  ggtitle(paste(chr_chosen, "\n150-200MB copmem, 300bp match",sep="")) +
  scale_color_brewer(palette = "Set1") + 
  theme_bw()


dev.off()

#copmem of inverted region

copmem.for.data <- read_tsv("/home/owens/working/reference_comparison/chr5_invertedregion.for.copmem.txt",comment=">",col_names=c("chr", "XRQ_pos","HA412_pos","length"))
copmem.for.data$direction <- "forward"
copmem.rev.data <- read_tsv("/home/owens/working/reference_comparison/chr5_invertedregion.rev.copmem.txt",comment=">",col_names=c("chr", "XRQ_pos","HA412_pos","length"))
copmem.rev.data$direction <- "reverse"

copmem.data <- rbind(copmem.for.data,copmem.rev.data)

copmem.data %>% filter(length > 330) %>%
  ggplot(.,aes(x=XRQ_pos,y=HA412_pos)) + geom_point(aes(color=direction)) +
  facet_wrap(~direction) +
  ggtitle(paste(chr_chosen, "\n150-200MB copmem, 300bp match",sep="")) +
  scale_color_brewer(palette = "Set1") + 
  theme_bw()


#blastn of inverted region

blast.data <- read_tsv("/home/owens/working/reference_comparison/chr5_invertedregion.blast.txt",
                       col_names=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

blast.data %>% mutate(qlength = qend - qstart, slength = send - sstart) %>% 
  mutate(direction= case_when(slength > 0 ~"forward", 
                      slength < 0 ~ "reverse")) %>% 
  filter(length >8000) %>%
#  filter(direction == "forward") %>%
#  filter(pident < 99) %>% 
  ggplot(.,aes(x=sstart,y=qstart)) + geom_point(aes(color=direction)) + facet_wrap(~direction)
                      

#copmem using masked genome.

copmem.for.data <- read_tsv("/home/owens/working/reference_comparison/ha412_sunflower_26Apr2018_CHus7_top17scaf.Scaffold_909.HRSCAF.1463.hardmsk.for.copmem.txt",comment=">",col_names=c("chr", "XRQ_pos","HA412_pos","length"))
copmem.for.data$direction <- "forward"
copmem.rev.data <- read_tsv("/home/owens/working/reference_comparison/ha412_sunflower_26Apr2018_CHus7_top17scaf.Scaffold_909.HRSCAF.1463.hardmsk.rev.copmem.txt",comment=">",col_names=c("chr", "XRQ_pos","HA412_pos","length"))
copmem.rev.data$direction <- "reverse"

copmem.rev.data$HA412_pos <- 190537180 - copmem.rev.data$HA412_pos


copmem.data <- rbind(copmem.for.data,copmem.rev.data)

pdf("Chr5_inversion_comparison.v0.pdf",height=15,width=15)
copmem.data %>% filter(length > 300) %>% View()
  ggplot(.,aes(x=XRQ_pos,y=HA412_pos)) + geom_point(aes(color=direction)) +
  ggtitle(paste("ha412_sunflower_26Apr2018_CHus7_top17scaf.Scaffold_909.HRSCAF.1463.hardmsk", "\n copmem, 300bp match",sep="")) +
  scale_color_brewer(palette = "Set1") + 
  theme_bw()
dev.off()

#Trying copmem of all masked ha412 chromosomes against all XRQ chromosomes
chr_lengths <- read_delim("/home/owens/working/reference_comparison/ha412_scaffold_lengths.txt",col_names = c("scaf","length"),delim=" ")

chr <- list.files("/home/owens/working/reference_comparison/",pattern="allcomp.rev.copmem.txt")
chr <- gsub(".allcomp.rev.copmem.txt","",chr)

pdf("XRQv2_vs_ha412_RedMasked_copmem300.v0.pdf",height=15,width=15)
for (chosen_chr in chr){
  tmp.for <- read_tsv(paste("/home/owens/working/reference_comparison/",chosen_chr,".allcomp.for.copmem.txt",sep=""),comment=">",col_names=c("chr", "XRQ_pos","HA412_pos","length"))
  tmp.for$direction <- "forward"
  tmp.rev <- read_tsv(paste("/home/owens/working/reference_comparison/",chosen_chr,".allcomp.rev.copmem.txt",sep=""),comment=">",col_names=c("chr", "XRQ_pos","HA412_pos","length"))
  tmp.rev$direction <- "reverse"
  tmp.length <- chr_lengths %>% filter(scaf == chosen_chr) %>% select(length) %>% .$length
  tmp.rev$HA412_pos <- tmp.length - tmp.rev$HA412_pos
  
  
  tmp <- rbind(tmp.for, tmp.rev)
  print(
    tmp %>% filter(length > 300) %>%
      filter(!grepl("00",chr)) %>% filter(!grepl("CP",chr)) %>% filter(!grepl("MT",chr)) %>%
      ggplot(.,aes(x=XRQ_pos,y=HA412_pos)) + geom_point(aes(color=direction)) +
      facet_wrap(~chr) +
      ggtitle(paste(chosen_chr, "\ncopmem, 300bp match, masked",sep="")) +
      scale_color_brewer(palette = "Set1") + 
      theme_bw()
  )
  matched_chr <-  tmp %>% mutate(XRQ_window = paste(chr, round(XRQ_pos / 10000))) %>% distinct(chr, XRQ_window) %>%
    group_by(chr) %>% summarize(count = n()) %>%  ungroup() %>% arrange(desc(count)) %>% head(n =1) %>% .$chr
  
  print(
    tmp %>% filter(length > 300) %>%
      filter(chr == matched_chr) %>%
      ggplot(.,aes(x=XRQ_pos,y=HA412_pos)) + geom_point(aes(color=direction)) +
      ggtitle(paste(matched_chr,"\n",chosen_chr, "\ncopmem, 300bp match, masked",sep="")) +
      scale_color_brewer(palette = "Set1") + 
      theme_bw()
    
  )
}
dev.off()



chosen_chr <- "ha412_sunflower_26Apr2018_CHus7_top17scaf.Scaffold_917.HRSCAF.1509"


