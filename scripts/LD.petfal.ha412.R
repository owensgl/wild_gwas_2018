#This prints out the summarized LD files. 

library(tidyverse)
library(viridis)


directory <- "/media/owens/Copper/wild_gwas_2018"
species <- "Petiolaris"
species_abr <- "petfal"
broad_species_abr <- "petiolaris"
prefix <- "tranche90.snp.petfal.90.bi.remappedHa412HO"

inv_locations <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt") %>% filter(spe == broad_species_abr)

inv_windows <- tibble(chr=character(),win1=numeric(),win2=numeric())
prev_mds <- "N"
prev_chr <- "N"
saved_windows <- ""
for (i in 1:nrow(inv_locations)){
  chosen_windows <- (floor(inv_locations$start[i]/500000):ceiling(inv_locations$end[i]/500000))*500000
  current_mds <- paste(inv_locations$chr[i],inv_locations$mds[i])
  chosen_chr <- inv_locations$chr[i]
  if (current_mds == prev_mds){
    saved_windows <- c(saved_windows, chosen_windows)
  }else{
    prev_mds = current_mds
    if (prev_chr != "N"){
      for (j in 1:length(saved_windows)){
        for (k in 1:length(saved_windows)){
          tmp <- tibble(chr=prev_chr,win1=saved_windows[j],win2=saved_windows[k])
          inv_windows <- rbind(inv_windows, tmp)
        }
      }
    }
    saved_windows <- chosen_windows
    prev_chr <- chosen_chr
  }
}
for (j in 1:length(saved_windows)){
  for (k in 1:length(saved_windows)){
    tmp <- tibble(chr=prev_chr,win1=saved_windows[j],win2=saved_windows[k])
    inv_windows <- rbind(inv_windows, tmp)
  }
}


all_data <- tibble(chr=character(),win1=numeric(),win2=numeric(),
                   n=numeric(),mean_r2=numeric(),max_r2=numeric(),
                   max_2_r2=numeric(),r2_min_0.2=numeric(),r2_min_0.8=numeric())
for (chr in sprintf("%02d", seq(1,17))){
  data <-  read_tsv(paste(directory,"/",species_abr,"/",species,".",prefix,".thin100.maf5.Chr",chr,".windows.ld.gz",sep=""))
  all_data <- rbind(all_data,data)
}

pdf(paste("LD/Ha412HO/",species_abr,"/",species,".",prefix,".thin100.maf5.ld.pdf",sep=""),height=15,width=15)
all_data %>%
  mutate(chr = str_replace(chr, "Ha412HO","")) %>%
  filter(n > 50) %>%
  ggplot(.,aes(x=win1/1000000,y=win2/1000000)) + geom_tile(aes(fill=max_2_r2)) +
  theme_bw() + scale_fill_viridis(name="2nd highest R2") +
  xlab("MB") + ylab("MB") +
  geom_tile(data=inv_windows %>% 
              mutate(chr = str_replace(chr, "Ha412HO",""))%>%
              filter(win1 > win2),
            aes(x=win1/1000000,y=win2/1000000),fill="red") +
  facet_wrap(~chr,scale="free") 
dev.off()

