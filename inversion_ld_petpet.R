library(tidyverse)

intrachr_r <- read_tsv("/media/owens/Copper/wild_gwas_2018/petpet/Petiolaris.tranche90.snp.petpet.90.bi.reduced.geno.ld")

pdf("Petiolaris.tranche90.snp.petpet.90.bi.reduced.geno.ld.v1.pdf",height=10,width=10)
full_chrs <- sprintf("HanXRQChr%02d", seq(17))
for (i in full_chrs){
  print(
    intrachr_r %>% rename(r2 = `R^2`,chr = CHR, pos1 = POS1, pos2 = POS2) %>% 
      filter(chr == i) %>% 
      mutate(win1 = floor(pos1/1000000), win2 = floor(pos2/1000000)) %>% 
      group_by(chr, win1, win2) %>%
      arrange(desc(r2)) %>%
      slice(2) %>% 
      ggplot(.,aes(x=win1,y=win2)) + geom_tile(aes(fill=r2)) +
      theme_bw() + scale_fill_viridis(name="2nd highest R2") +
      xlab("MB") + ylab("MB") +
      ggtitle(i)
  )
  
}
dev.off() 


###Interchromosome
interchr_r <- read_tsv("/media/owens/Copper/wild_gwas_2018/niveus/Petiolaris.tranche90.snp.gwas.90.bi.reduced.interchrom.geno.ld")
interchr_r %>% rename(r2 = `R^2`,chr1 = CHR1, pos1 = POS1, chr2 = CHR2, pos2 = POS2) -> interchr_r

chrs <- c("HanXRQChr03","HanXRQChr06","HanXRQChr10")

pdf("Petiolaris.tmp.pet.snp.bi.90.filtered.reduced.geno.interld.v1.pdf",height=10,width=10)
for (i in 1:(length(chrs)-1)){
  for (j in (i+1):length(chrs)){
    print(
      interchr_r %>%
        filter(chr1 == chrs[i],chr2 == chrs[j]) %>%
        mutate(win1 = floor(pos1/1000000), win2 = floor(pos2/1000000)) %>% 
        group_by(chr1,chr2, win1, win2) %>%
        arrange(desc(r2)) %>%
        slice(2) %>%
        ggplot(.,aes(x=win1,y=win2)) + geom_tile(aes(fill=r2)) + 
        theme_bw() +
        xlab(paste(chrs[i],"MB")) + ylab(paste(chrs[j],"MB")) +
        ggtitle(paste(chrs[i],"vs",chrs[j]))+
        scale_fill_viridis(name="2nd highest R2")
    )
  }
}
dev.off()

