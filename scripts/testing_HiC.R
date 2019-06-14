library(tidyverse)
library(viridis)

chrs <- c("01","05")
pdf("HiC/HA412_2lane.filtered.HiCtest.v1.pdf")
for (chr in chrs){
  
  HiC <- read_tsv(paste0("HiC/HA412_2lane.filtered.chr",chr,".1000k.raw.txt.gz"))
  
  print(
    HiC %>% dplyr::select(-Regions) %>%
      gather(.,window2, value, 2:ncol(.)) %>%
      rename(window1 = `HiCMatrix (directory=HA412_2lane_filtered)`) %>%
      mutate(window1 = as.numeric(str_remove(window1, paste0("HanPSC8Chr",chr,"-"))),
             window2 = as.numeric(str_remove(window2, paste0("HanPSC8Chr",chr,"-")))) %>%
      filter(window1 != window2) %>%
      ggplot(.,aes(x=window1/1000000,y=window2/1000000,fill=value)) + geom_tile() +
      scale_fill_viridis(option="viridis",limits=c(),name="Interactions") +
      theme_bw() + ylab("MB") + xlab("MB") +
      ggtitle(paste0("HaPSC8Chr",chr," - Ha412HO HiC"))
  )
  
  
  HiC %>% dplyr::select(-Regions) %>%
    gather(.,window2, value, 2:ncol(.)) %>%
    rename(window1 = `HiCMatrix (directory=HA412_2lane_filtered)`) %>%
    mutate(window1 = as.numeric(str_remove(window1, paste0("HanPSC8Chr",chr,"-"))),
           window2 = as.numeric(str_remove(window2, paste0("HanPSC8Chr",chr,"-")))) %>%
    filter(window1 != window2) -> HiC_processed
  
  if (chr == "01"){
    HiC_processed %>%
      filter(window1 < 12000000, window2 < 12000000, window1 > 000000, window2 > 000000) -> HiC_processed
  }else if (chr == "05"){
    HiC_processed %>%
      filter(window1 < 190000000, window2 < 190000000, window1 > 140000000, window2 > 140000000)-> HiC_processed
  }
  print(
    HiC_processed %>%
      ggplot(.,aes(x=window1/1000000,y=window2/1000000,fill=value)) + geom_tile() +
      scale_fill_viridis(option="viridis",limits=c(),name="Interactions") +
      theme_bw() + ylab("MB") + xlab("MB") +
      ggtitle(paste0("HaPSC8Chr",chr," - Ha412HO HiC"))
  )
  print(
    HiC_processed %>%
      ggplot(.,aes(x=window2/1000000,y=value)) + geom_point() + 
      facet_wrap(~as.factor(paste(window1/1000000,"MB")),scales="free_y") +
      geom_vline(aes(xintercept = window1/1000000)) +
      ylab("Interacations") +
      xlab("MB") + 
      ggtitle(paste0("HaPSC8Chr",chr," - Ha412HO HiC"))
  )
  HiC %>% dplyr::select(-Regions) %>%
    gather(.,window2, value, 2:ncol(.)) %>%
    rename(window1 = `HiCMatrix (directory=HA412_2lane_filtered)`) %>%
    mutate(window1 = as.numeric(str_remove(window1, paste0("HanPSC8Chr",chr,"-"))),
           window2 = as.numeric(str_remove(window2, paste0("HanPSC8Chr",chr,"-")))) %>%
    filter(window1 != window2) %>% pull(window1) %>% max() -> chr_end
  HiC %>% dplyr::select(-Regions) %>%
    gather(.,window2, value, 2:ncol(.)) %>%
    rename(window1 = `HiCMatrix (directory=HA412_2lane_filtered)`) %>%
    mutate(window1 = as.numeric(str_remove(window1, paste0("HanPSC8Chr",chr,"-"))),
           window2 = as.numeric(str_remove(window2, paste0("HanPSC8Chr",chr,"-")))) %>%
    filter(window1 != window2) %>% pull(window1) %>% min() -> chr_start
  
  
  print(
    HiC %>% dplyr::select(-Regions) %>%
      gather(.,window2, value, 2:ncol(.)) %>%
      rename(window1 = `HiCMatrix (directory=HA412_2lane_filtered)`) %>%
      mutate(window1 = as.numeric(str_remove(window1, paste0("HanPSC8Chr",chr,"-"))),
             window2 = as.numeric(str_remove(window2, paste0("HanPSC8Chr",chr,"-")))) %>%
      filter(window1 != window2) %>%
      mutate(distance = abs(window1 - window2),
             neighbours = case_when(distance <= 2000000 ~ "close",
                                    TRUE ~ "far")) %>%
      group_by(window1) %>%
      mutate(max_value = max(value[which(neighbours == "close")])) %>%
      filter(neighbours == "far") %>%
      filter(value > max_value,value > 10) %>%
      ggplot(.,aes(x=window1/1000000,xend=window2/1000000,y=1,yend=1)) + 
      geom_curve(curvature=0.5) +
      geom_segment(aes(x=chr_start/1000000,xend=chr_end/1000000,y=1,yend=1),color="grey",alpha=0.1) +
      theme_bw() +
      ylab("HaPSC8") + xlab("MB") + 
      ggtitle(paste0("HaPSC8Chr",chr," - Ha412HO HiC\nLong-distance interactions"))
  )
}
dev.off()

