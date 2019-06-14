library(tidyverse)
library(zoo)
library(SNPRelate)
library(grid)
library(gridExtra)
library(plotly)
library(scatterpie)
library(ggpmisc)
library(viridis)
#paper_colors <- c("#5BC0EB","#9BC53D","#FDE74C","#E55934","#FA7921")
paper_colors <- c("#004C6C","#7FB800","#FFE3A2","#F6511D","#0D2C54")

#Make final figures for inversion paper
#Part 1: The MDS process.
inversions <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt")
data_directory <- "/media/owens/Copper/wild_gwas"
window_size <- "100"
chr_prefix <- "Ha412HOChr"
het_script <- "/home/owens/bin/pop_gen/vcf2het.pl"
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T) 
pop_loc <- read_tsv("pop_loc_allnum.txt")
pop_loc %>% rename(population = pop) %>% inner_join(.,labels) %>% rename(sample = name) -> labels
info <- read_tsv("sample_info_apr_2018.tsv") %>% rename(sample = name)


#Get inversion locations

directory <- "/media/owens/Copper/wild_gwas"
species <- "Annuus"
species_abr <- "annuus"
prefix <- "tranche90.snp.env.90.bi.remappedHa412HO"

inv_locations <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt") %>% filter(spe == species_abr)

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





n <- 2

genotypes <- tibble(sample=character(),EV1=character(),
                    EV2=character(),cluster=character(),
                    species=character(),chr=character(),
                    mds=character())
regions <- tibble(group=character(),start=character(),
                  end=character(),chr=character(),
                  spe=character(), species=character(),
                  dataset=character(),mds=character())

for (n in 1:nrow(inversions)){
  
  win.regions <- readRDS(paste(
    "MDS_plots/Ha412HO/",
    pull(inversions[n,1]),
    "/",
    pull(inversions[n,2]),
    ".tranche90.snp.",
    pull(inversions[n,3]),
    ".90.bi.remappedHa412HO.Chr",
    sprintf("%02d",pull(inversions[n,4])),
    ".",
    window_size,
    ".windows.rds",sep=""))
  
  win.regions %>%
    select(chrom,start,end,n,mid,paste("mds",sprintf("%02d",pull(inversions[n,5])),sep="")) %>%
    rename(mds_chosen = paste("mds",sprintf("%02d",pull(inversions[n,5])),sep=""))-> win.regions
  threshold <- pull(inversions[n,7])
  if (pull(inversions[n,6]) == "neg"){
    win.regions <- win.regions %>%
      mutate(outlier = case_when(mds_chosen < threshold ~ "1",
                                 TRUE ~ "0")) 
  }else{
    win.regions <- win.regions %>%
      mutate(outlier = case_when(mds_chosen > threshold ~ "1",
                                 TRUE ~ "0")) 
  }
  
  #Select windows surrounded by other outliers.
  max_distance_surround <- 20
  win.regions <- win.regions %>% mutate(lead_max = rollmax(x = outlier, max_distance_surround, align = "left", fill = NA),
                                        lag_max = rollmax(x = outlier, max_distance_surround, align = "right", fill = NA)) %>%
    mutate(filled_outlier = case_when(outlier == 1 ~ 1,
                                      lead_max == 1 & lag_max == 1 ~ 1,
                                      TRUE ~ 0)) 
  #Plot mds scores
  plot_mds <- win.regions %>%
    ggplot(.,aes(x=mid/1000000,y=mds_chosen,color=as.factor(filled_outlier))) + geom_point(size=0.7) +
    theme_bw() + scale_color_manual(values=c("black",paper_colors[4]),name="Selected") +
    ylab(paste("mds",sprintf("%02d",pull(inversions[n,5])),sep="")) +
    xlab(paste0("Chr",sprintf("%02d",pull(inversions[n,4]))," (Mb)")) +
    theme(legend.position = "none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) +
    labs(tag="A") +
    scale_y_continuous(breaks=c(0,0.3,0.6)) +
    scale_x_continuous(expand = c(0.01, 0.01))

  
  win.regions$ungapped_n <- 1:nrow(win.regions)
  win.regions <- win.regions %>% 
    filter(filled_outlier == 1) %>% mutate(gap = ungapped_n - lag(ungapped_n))
  
  win.regions$group <- "1"
  current_group = 1
  for (i in 2:nrow(win.regions)){
    if (win.regions$gap[i] > 1){
      current_group = current_group + 1
    }
    win.regions$group[i] <- current_group
  }
  win.regions %>% 
    group_by(group) %>%
    summarize(start = min(start),end=max(end)) %>%
    mutate(chr=paste(chr_prefix,sprintf("%02d",pull(inversions[n,4])),sep="")) -> region_selected
  #Save regions to output file
  region_output <- region_selected
  region_output$spe <- pull(inversions[n,1])
  region_output$species <- pull(inversions[n,2])
  region_output$dataset <-  pull(inversions[n,3])
  region_output$mds <- paste(pull(inversions[n,6]),pull(inversions[n,5]),sep="")
  regions <- rbind(regions, region_output)
  
  vcf_in <- paste(
    data_directory,
    "/",
    pull(inversions[n,1]),
    "/",
    pull(inversions[n,2]),
    ".tranche90.snp.",
    pull(inversions[n,3]),
    ".90.bi.remappedHa412HO.vcf.gz",
    sep="")
  vcf_out <- paste(
    data_directory,
    "/",
    pull(inversions[n,1]),
    "/",
    pull(inversions[n,2]),
    ".tranche90.snp.",
    pull(inversions[n,3]),
    ".90.bi.remappedHa412HO.tmp.vcf.gz",
    sep="")
  gds_out <- paste(
    data_directory,
    "/",
    pull(inversions[n,1]),
    "/",
    pull(inversions[n,2]),
    ".tranche90.snp.",
    pull(inversions[n,3]),
    ".90.bi.remappedHa412HO.tmp.gds",
    sep="")
  region_string <- paste(region_selected$chr[1],":",region_selected$start[1],"-",region_selected$end[1],sep="")
  if (nrow(region_selected)> 1){
    for (i in 2:nrow(region_selected)){
      region_string <- paste(region_string,",",region_selected$chr[i],":",region_selected$start[i],"-",region_selected$end[i],sep="")
    }
  }
  
  
  system(paste("bcftools view -O z -r",region_string,vcf_in,">",vcf_out))
  system(paste("tabix -p vcf",vcf_out))
  
  system(paste("zcat ",vcf_out," | perl ",het_script," > tmp.het.txt",sep=""))
  vcf_het <- read_tsv("tmp.het.txt")
  snpgdsVCF2GDS(vcf_out,gds_out, method="biallelic.only",ignore.chr.prefix="Ha412HOChr")
  genofile <- snpgdsOpen(gds_out)
  pca <- snpgdsPCA(genofile, num.thread=2)
  tab <- data.frame(sample = pca$sample.id,
                    EV1 = pca$eigenvect[,1],    # the first eigenvector
                    EV2 = pca$eigenvect[,2],    # the second eigenvector
                    stringsAsFactors = FALSE)
  snpgdsClose(genofile)
  samplelist <- as.vector(tab$sample)
  
  try_3_clusters <-try(kmeans(tab[,2], 3, centers=c(min(tab[,2]),(min(tab[,2])+max(tab[,2]))/2,max(tab[,2]))))
  
  
  if("try-error" %in% class(try_3_clusters)){
    kmeans_cluster <-kmeans(tab[,2], 2, centers=c(min(tab[,2]),max(tab[,2])))
    
  }else{
    
    kmeans_cluster <-kmeans(tab[,2], 3, centers=c(min(tab[,2]),(min(tab[,2])+max(tab[,2]))/2,max(tab[,2])))
    
  }
  
  tab$cluster <- kmeans_cluster$cluster - 1 
  
  betweenss <- kmeans_cluster$betweenss
  
  max_het <- round(max(vcf_het$percent_het),2)
  min_het <- round(min(vcf_het$percent_het),2)
  
  het_plot <- vcf_het %>%
    inner_join(tab) %>%
    ggplot(.,aes(x=as.factor(cluster),y=percent_het)) + geom_boxplot(aes(fill=as.factor(cluster))) +
    scale_fill_manual(values=paper_colors,name="Kmeans\ncluster") +
    theme_bw() +
    theme(legend.position = "none") +
    ylab("Heterozygosity") +
    labs(tag = "C")
  
  if (n == 1){
  data.tb <- tibble(x = 0.1, y = 0.18, 
                    plot = list(het_plot +
                                  theme_bw(10) +
                                  theme(legend.position = "none",
                                        panel.grid.major = element_blank(), 
                                        panel.grid.minor = element_blank(),
                                        axis.title.x=element_blank(),
                                        axis.text.x=element_blank(),
                                        axis.ticks.x=element_blank())))
  }else if (n == 2){
    data.tb <- tibble(x = 0.06, y = -0.2, 
                      plot = list(het_plot +
                                    theme_bw(10) +
                                    theme(legend.position = "none",
                                          panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank(),
                                          axis.title.x=element_blank(),
                                          axis.text.x=element_blank(),
                                          axis.ticks.x=element_blank())))
  }
  

  plot_pca <- 
    tab %>% 
    ggplot(.,aes(x=EV1,y=EV2)) + geom_point(aes(color=as.factor(cluster)),size=2) +
    geom_plot(data = data.tb, aes(x, y, label = plot,vp.width=0.7)) +
      
    theme_bw() + 
    #ggtitle(paste(" BetweenSS = ",betweenss, sep="")) + 
    scale_color_manual(values=paper_colors,name="Kmeans\ncluster") +
    theme(legend.position = "none") +
    ylab("PC2") + xlab("PC1") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(tag = "B")
    
  
  # tab %>% 
  #   inner_join(.,info) %>%
  #   ggplot(.,aes(x=EV1,y=EV2)) + geom_point(aes(color=species,shape=as.factor(cluster)),size=2) +
  #   theme_bw() + ggtitle(paste(" BetweenSS = ",betweenss, sep="")) + 
  #   scale_color_brewer(palette = "Set1") +
  #   scale_shape(name="Kmeans\ncluster") +
  #   theme(legend.position="bottom")
  
  
  #Define the map for each type
  if (pull(inversions[n,1]) == "annuus"){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(25, 50)
    long_range <- c(-125,-93)
    pie_size <- 0.4
  }else if(pull(inversions[n,1]) == "argophyllus" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- subset(states, region %in% c("texas"))
    lat_range <- c(25, 30)
    long_range <- c(-99,-96)
    pie_size = 0.1
  }else if(pull(inversions[n,1]) == "petpet" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(30, 43)
    long_range <- c(-115,-95)
    pie_size <- 0.4
  }else if(pull(inversions[n,1]) == "petfal" ){
    usa <- map_data('state')
    states <- map_data("state")
    target_state <- map_data('state')
    lat_range <- c(30, 40)
    long_range <- c(-115,-100)
    pie_size = 0.4
  }
  
  
  if(length(unique(kmeans_cluster$cluster)) == 3){
    plot_map <- ggplot(target_state, aes(long, lat)) +
      geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
      coord_quickmap() +
      geom_scatterpie(data=inner_join(tab, labels) %>% 
                        group_by(population, lat, long, cluster) %>% 
                        tally() %>%
                        spread(., cluster, n,fill=0),
                      aes(x=long, y=lat, r=pie_size), 
                      cols=c("0","1","2"), color=NA, alpha=.9) +
      scale_fill_manual(name="Cluster",values=paper_colors) +theme_bw() +
      xlab("Longitude") + ylab("Latitude") +
      scale_x_continuous(limits = long_range, expand = c(0, 0)) +
      scale_y_continuous(limits = lat_range, expand = c(0, 0)) +
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank()) +
      labs(tag = "D")
  }else{
    plot_map <- ggplot(target_state, aes(long, lat)) +
      geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
      coord_quickmap() +
      geom_scatterpie(data=inner_join(tab, labels) %>% 
                        group_by(population, lat, long, cluster) %>% 
                        tally() %>%
                        spread(., cluster, n,fill=0),
                      aes(x=long, y=lat, r=0.4), 
                      cols=c("0","1"), color=NA, alpha=.8) +
      scale_fill_brewer(name="Cluster",palette = "Set1") +theme_bw() +
      xlab("Longitude") + ylab("Latitude") +
      scale_x_continuous(limits = long_range, expand = c(0, 0)) +
      scale_y_continuous(limits = lat_range, expand = c(0, 0))
  }
  
  if (n == 2){
    coords_columns <- c("region1_start","region1_end","region2_start","region2_end","length1","length2","identity","direction1","direction2","ref1","ref2")
    
    coords_directory <- "/home/owens/working/reference_comparison"
    
    coords <- read_tsv(paste0(coords_directory,"/Ha412HO-HaPSC8/Ha412HO-HaPSC8.Chr05.b1000.c200.g500.1coords"),
                       col_names = coords_columns) 
    x_plot <- coords %>%
      filter(length1 > 1000,region1_start > 130000000,region2_start > 130000000) %>%
      mutate(ID = 1:n()) %>%
      mutate(direction = case_when( region2_start > region2_end ~ "reverse",
                                    TRUE ~ "forward")) %>%
      select(region1_start,region1_end,region2_start,region2_end,ID,direction) %>%
      gather( location, xpos, 1:4) %>%
      mutate(ypos = case_when( (location == "region1_start") ~ "Ha412HO",
                               (location == "region1_end") ~ "Ha412HO",
                               (location == "region2_start") ~ "HaPSC8",
                               (location == "region2_end") ~ "HaPSC8")) %>%
      mutate(location = fct_relevel(location, "region1_start","region1_end","region2_end","region2_start")) %>%
      arrange(location) %>%
      filter(direction == "forward") %>%
      ggplot(.,aes(x=xpos/1000000,y=fct_relevel(ypos,"HaPSC8","Ha412HO"),group=ID)) + 
      geom_polygon(aes(),fill=paper_colors[5],alpha=0.1) +
      geom_polygon(data=coords %>%
                     filter(length1 > 1000,region1_start > 130000000,region2_start > 130000000) %>%
                     mutate(ID = 1:n()) %>%
                     mutate(direction = case_when( region2_start > region2_end ~ "reverse",
                                                   TRUE ~ "forward")) %>%
                     select(region1_start,region1_end,region2_start,region2_end,ID,direction) %>%
                     gather( location, xpos, 1:4) %>%
                     mutate(ypos = case_when( (location == "region1_start") ~ "Ha412HO",
                                              (location == "region1_end") ~ "Ha412HO",
                                              (location == "region2_start") ~ "HaPSC8",
                                              (location == "region2_end") ~ "HaPSC8")) %>%
                     mutate(location = fct_relevel(location, "region1_start","region1_end","region2_end","region2_start")) %>%
                     arrange(location) %>%
                     filter(direction == "reverse"),
                     aes(),fill=paper_colors[4]) +
      theme_bw() +
      xlab("Mb") +
      coord_cartesian(ylim=c(1.5,1.5)) +
      scale_x_continuous(breaks=c(140,160,180)) +
      theme(legend.position="none",
            axis.text.y = element_text(size = 1)) +
      labs(tag = "F")
    

    
    ld_data <-read_tsv(paste(directory,"/","annuus","/Annuus.tranche90.snp.env.90.bi.remappedHa412HO.thin100.maf5.Chr05.windows.ld.gz",sep=""))

    ld_plot <- ld_data %>%
      filter(chr =="Ha412HOChr05") %>%
      mutate(chr = str_replace(chr, "Ha412HO","")) %>%
      filter(n > 50) %>%
      ggplot(.,aes(x=win1/1000000,y=win2/1000000)) + geom_tile(aes(fill=max_2_r2)) +
      theme_bw() + scale_fill_viridis(name="2nd highest R2") +
      xlab("Mb") + ylab("Mb") +
      geom_tile(data=inv_windows %>% 
                  filter(chr =="Ha412HOChr05") %>%
                  mutate(chr = str_replace(chr, "Ha412HO",""))%>%
                  filter(win1 > win2),
                aes(x=win1/1000000,y=win2/1000000),fill=paper_colors[4]) +
      theme(legend.position = "none",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.margin=unit(c(0,0,0,0),"cm")) +
      labs(tag = "G")
    
    
    insert_plots.tb <- tibble(x = c(0,200), y = c(200,0), 
                        plot = c(list(x_plot +
                                      theme_bw(10) +
                                      theme(legend.position = "none",
                                            panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank(),
                                            panel.background = element_blank(),
                                            plot.margin=unit(c(0,0,0,0),"cm"),
                                            axis.text.y = element_text(size = 5),
                                            axis.title.y=element_blank())),
                                  list(ld_plot +
                                        theme_bw(10) +
                                        theme(legend.position = "none",
                                              panel.grid.major = element_blank(), 
                                              panel.grid.minor = element_blank(),
                                              panel.background = element_blank(),
                                              plot.margin=unit(c(0,0,0,0),"cm")))))
                        
    
    
    plot_genome <- coords %>%
      filter(length1 > 1000) %>%
      mutate(direction = case_when( region2_start < region2_end ~ "reverse",
                                    TRUE ~ "forward")) %>% 
      ggplot(.,aes(x=region1_start/1000000,xend=region1_end/1000000,
                    y=region2_start/1000000,yend=region2_end/1000000,
                    color=as.factor(direction))) +
      geom_segment(size=1.4) +
      geom_plot(data = insert_plots.tb, aes(x, y, label = plot, vp.width=0.4,vp.height=0.4)) +
      scale_color_manual(values=c(paper_colors[4],paper_colors[5]),labels=c("Reverse","Forward"),
                         name="Direction") +
      theme_bw() +
      ylab("HaPSC8 Mb") + xlab("Ha412HO Mb") +
      theme(legend.position="bottom",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
      labs(tag = "E")
      
    
    
    
  }  
    pdf("Figure_01v01",width=12,height=6)
  print(
    grid.arrange( plot_mds,plot_pca,plot_map,plot_genome,
                  #nrow=3,
                  #top = textGrob(paste(pull(inversions[n,1])," ",region_string,sep=""),gp=gpar(fontsize=10,font=2)),
                  layout_matrix = rbind(c(1,1,1),
                                        c(2,3,4),
                                        c(2,3,4),
                                        c(2,3,4),
                                        c(2,3,4)))
  )
  dev.off()

}
