library(tidyverse)
library(ggthemes)

pdf("Annuus.tranche90.snp.env.90.bi.remappedHa412.maf10.thin5k.500kbwin.ld.pdf",height=10,width=10)
for (i in sprintf("%02d", seq(17))){
  
  ld <- read_tsv(paste("/media/owens/Copper/wild_gwas_2018/annuus/Annuus.tranche90.snp.env.90.bi.remappedHa412.maf10.thin5k.Chr",i,".geno.ld",sep=""),
                 col_names = c("chr","pos1","pos2","n","r2"),skip=1)
  window_size <- 500000
  print(
  ld %>%
    mutate(win1 = floor(pos1/window_size)*window_size, win2 = floor(pos2/window_size)*window_size) %>% 
    group_by(win1, win2) %>%
    arrange(desc(r2)) %>%
    slice(2) %>% 
    ggplot(.,aes(x=win1/1000000,y=win2/1000000)) + geom_tile(aes(fill=r2)) +
    theme_bw() + scale_fill_viridis_c(name="2nd highest R2") +
    xlab("MB") + ylab("MB") +
    ggtitle(paste("HA412TMPChr",i,sep=""))
  )
  
}
dev.off()

pdf("Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412.maf5.thin1k.500kbwin.ld.pdf",height=10,width=10)
for (i in sprintf("%02d", seq(17))){
  
  ld <- read_tsv(paste("/media/owens/Copper/wild_gwas_2018/argophyllus/Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412.maf5.thin1k.Chr",i,".geno.ld",sep=""),
                 col_names = c("chr","pos1","pos2","n","r2"),skip=1)
  window_size <- 500000
  print(
    ld %>%
      mutate(win1 = floor(pos1/window_size)*window_size, win2 = floor(pos2/window_size)*window_size) %>% 
      group_by(win1, win2) %>%
      arrange(desc(r2)) %>%
      slice(2) %>% 
      ggplot(.,aes(x=win1/1000000,y=win2/1000000)) + geom_tile(aes(fill=r2)) +
      theme_bw() + scale_fill_viridis_c(name="2nd highest R2") +
      xlab("MB") + ylab("MB") +
      ggtitle(paste("HA412TMPChr",i,sep=""))
  )
  
}
dev.off()

ld_files <- list.files("/media/owens/Copper/wild_gwas_2018/annuus/",pattern="window.ld") %>% .[1:153]

ld_data <- ld_files %>%
  map_dfr(function(x) { 
    read_tsv(file.path("/media/owens/Copper/wild_gwas_2018/annuus/", x)) %>% mutate(file=x) } ) %>%
  separate(.,file,into=c("species","tranche","snp","env","n","bi","beagle","remap","pop","chrom","window","ld"),sep="\\.") %>%
  select(chr,win1,win2,mean_r2, r2_min_0.2, pop)
  


  
pdf("Annuus.tranche90.snp.env.90.bi.beagle.remappedHa412.tmp.ld.pdf",width=15,height=6)
for (n in 1:17){
  n_chr <- paste("Ha412TMPChr", str_pad(n, 2, pad = "0"),sep="")
  print(
  ld_data %>%
    filter(chr == n_chr) %>% 
    ggplot(.,aes(x=win1/1000000,y=win2/1000000,fill=mean_r2)) +
    geom_tile() +
    scale_fill_viridis_c(limits=c(0,0.5)) + theme_few() +
    facet_grid(chr~pop)
  )
}
dev.off()
