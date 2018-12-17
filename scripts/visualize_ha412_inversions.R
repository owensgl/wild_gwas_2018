library(tidyverse)
library(zoo)
library(SNPRelate)
library(grid)
library(gridExtra)
#Visualize candidate inversions from lostruct

inversions <- read_tsv("Ha412HO_inv.dec11.txt")
data_directory <- "/media/owens/Copper/wild_gwas_2018"
window_size <- "100"
chr_prefix <- "Ha412HOChr"
het_script <- "/home/owens/bin/pop_gen/vcf2het.pl"
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T) 
pop_loc <- read_tsv("pop_loc_allnum.txt")
pop_loc %>% rename(population = pop) %>% inner_join(.,labels) %>% rename(sample = name) -> labels


n <- 1

genotypes <- tibble(sample=character(),EV1=character(),
                    EV2=character(),cluster=character(),
                    species=character(),chr=character(),
                    mds=character())
regions <- tibble(group=character(),start=character(),
                  end=character(),chr=character(),
                  spe=character(), species=character(),
                  dataset=character(),mds=character())
pdf("Ha412HO_inv.dec11.inversions.v1.pdf",height=10,width=15)
for (n in 1:nrow(inversions)){
  
  win.regions <- readRDS(paste(
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
    ggplot(.,aes(x=mid/1000000,y=mds_chosen,color=as.factor(filled_outlier))) + geom_point() +
    theme_bw() + scale_color_brewer(palette = "Set1",name="Selected") +
    ylab(paste("mds",sprintf("%02d",pull(inversions[n,5])),sep="")) +
    xlab("Mb")
  
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

  plot_pca <- tab %>% 
    inner_join(.,vcf_het) %>%
    ggplot(.,aes(x=EV1,y=EV2)) + geom_point(aes(color=percent_het,shape=as.factor(cluster)),size=2) +
    theme_bw() + ggtitle(paste(" BetweenSS = ",betweenss, sep="")) + 
    scale_color_viridis_c(name="Heterozygosity",
                          limits=c(min_het,
                                   max_het),
                          breaks=c(min_het,
                                   max_het),
                          labels=c(min_het,
                                   max_het)) +
    scale_shape(name="Kmeans\ncluster") +
    theme(legend.position="bottom")
  
  
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
                      cols=c("0","1","2"), color=NA, alpha=.8) +
      scale_fill_brewer(name="Cluster",palette = "Set1") +theme_bw() +
      xlab("Longitude") + ylab("Latitude") +
      scale_x_continuous(limits = long_range, expand = c(0, 0)) +
      scale_y_continuous(limits = lat_range, expand = c(0, 0))
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
  
  print(
    grid.arrange( plot_mds,plot_pca,plot_map,
                  #nrow=3,
                  top = textGrob(paste(pull(inversions[n,1])," ",region_string,sep=""),gp=gpar(fontsize=10,font=2)),
                  layout_matrix = rbind(c(1,1),
                                        c(2,3),
                                        c(2,3),
                                        c(2,3),
                                        c(2,3)))
  )
  tab$species <- pull(inversions[n,1])
  tab$chr <- paste(chr_prefix,sprintf("%02d",pull(inversions[n,4])),sep="")
  tab$mds <- paste(pull(inversions[n,6]),pull(inversions[n,5]),sep="")
  genotypes <- rbind(genotypes,tab)
}
dev.off()
write_tsv(genotypes,"Ha412HO_inv.dec11.inversions.genotypes.v1.txt")
write_tsv(regions,"Ha412HO_inv.dec11.inversions.regions.v1.txt")

chrlengths <- read_tsv("Ha412HO.chrlengths.txt")

pdf("Ha412HO_inv.dec11.inversions.regions.v1.pdf",height=2,width=20)
regions %>%
  filter(spe != "argophyllus" | chr !="Ha412HOChr14" | mds != "pos1") %>%
  filter(spe != "petpet" | chr !="Ha412HOChr16" | mds != "neg1") %>%
  filter(spe != "petpet" | chr !="Ha412HOChr07" | mds != "neg1") %>%
  filter(spe != "petpet" | chr !="Ha412HOChr04" | mds != "pos1") %>%
  mutate(subspecies = case_when(spe == "petpet" ~ "H. petiolaris petiolaris",
                                spe == "petfal" ~ "H. petiolaris fallax",
                                spe == "annuus" ~ "H. annuus",
                                spe == "argophyllus" ~ "H. argophyllus")) %>%
  mutate(chr = gsub("Ha412HO","",chr)) %>%
  ggplot(.) + 
  geom_segment(aes(x=start/1000000,xend=end/1000000,y=as.factor(subspecies),yend=as.factor(subspecies),color=subspecies),size=3) +
  geom_segment(data=chrlengths %>% mutate(chr = gsub("Ha412HO","",chr)),aes(x=start/1000000,xend=end/1000000,y=1,yend=1),alpha=0) +
  facet_wrap(~chr,scales="free_x",nrow=1) +theme_bw() +
  scale_color_brewer(palette = "Set1",name="Subspecies") + 
  theme(legend.position="bottom",axis.title.y=element_blank()) +
  xlab("Mb") 
dev.off()

#Plotting Fst peaks. Collected using 
#perl /home/owens/bin/pop_gen/run_vcf2inversionfst.pl /home/owens/bin/wild_gwas_2018/Ha412HO_inv.dec11.inversions.genotypes.v1.txt
library(changepoint)
#/media/owens/Copper/wild_gwas_2018/argophyllus/Ha412HO_inv.dec11.argophyllus.Ha412HOChr03.pos2.fst.txt.gz

fst_regions <- tibble(group=character(),start=character(),
                      end=character(),chr=character(),
                      spe=character(),species=character(),
                      dataset=character(),mds=character())
pdf("Ha412HO_inv.dec11.inversions.fstregions.v1.pdf",width=15,height=10)
for (n in 1:nrow(inversions)){
  n_clusters <- inversions$n_blocks[n]
  fst <- read_tsv(paste(
    "/media/owens/Copper/wild_gwas_2018/",
    pull(inversions[n,1]),
    "/Ha412HO_inv.dec11.",
    pull(inversions[n,1]),
    ".Ha412HOChr",
    sprintf("%02d",pull(inversions[n,4])),
    ".",
    pull(inversions[n,6]),
    pull(inversions[n,5]),
    ".fst.txt.gz",sep="")) 
  if (nrow(fst) == 0){
    next
  }
  fst <- fst %>%
    mutate(Fst = case_when(Fst < 0 ~ 0,TRUE ~ Fst))
  
  
  #Run changepoint analysis to pull out elevated regions
  fst_cpt<- cpt.mean(fst$Fst,method='BinSeg',Q=(2 * n_clusters))
  
  
  changepoints <- sort(cpts.full(fst_cpt)[(2 * n_clusters),])
  fst$num <- 1:nrow(fst)
  
  fst$cluster <- 0
  for (i in 1:n_clusters){
    fst %>%
      mutate(cluster = case_when(num >= changepoints[(2 * i) - 1] &
                                   num < changepoints[(2 * i)] ~ 1,
                                 TRUE ~ cluster)) -> fst
  }
  
  
  window_size = 1000000
  full_plot <- fst %>%
    mutate(window = (floor(pos/window_size)*window_size)+(window_size/2)) %>%
    group_by(window) %>%
    summarize(window_fst = sum(FstNum)/sum(FstDenom)) %>%
    mutate(window_fst = case_when(window_fst < 0 ~ 0, TRUE ~ window_fst)) %>%
    ggplot(.,aes(x=window/1000000,y=window_fst)) + 
    geom_point(data=fst %>% filter(Fst > 0.05),
               aes(x=pos/1000000, y=Fst,color=as.factor(cluster)),alpha=0.05) +
    geom_line(color="black") +
    scale_color_brewer(palette = "Set1",name="Cluster") +theme_bw() +
    xlab("Mb") + ylab("Fst") +
    guides(color=FALSE) +
    coord_cartesian(ylim=c(0,1))
  
  edge_plots <- list()
  for (i in 1:n_clusters){
    
    left_edge <- fst %>% filter(num == changepoints[(2 * i) - 1]) %>% pull(pos)
    right_edge <- fst %>% filter(num == (changepoints[(2 * i)]-1)) %>% pull(pos)
    
    pulled_regions <- tibble(group="1",start=left_edge,
                             end=right_edge,chr=paste("Ha412HOChr",sprintf("%02d",pull(inversions[n,4])),sep=""),
                             spe=pull(inversions[n,1]),species=pull(inversions[n,2]),
                             dataset=pull(inversions[n,3]),
                             mds=paste(pull(inversions[n,6]),pull(inversions[n,5]),sep=""))
    fst_regions <- rbind(fst_regions,pulled_regions)
    
    
    window_size = 10000
    print_surround <- 500000
    plot1 <- fst %>%  
      mutate(window = (floor(pos/window_size)*window_size)+(window_size/2)) %>%
      mutate(Fst = case_when(Fst < 0 ~ 0,
                             TRUE ~ Fst)) %>%
      group_by(window) %>%
      summarize(window_fst = sum(FstNum)/sum(FstDenom)) %>%
      ggplot(.,aes(x=window,y=window_fst)) + 
      geom_point(data=fst,
                 aes(x=pos, y=Fst,color=as.factor(cluster)),alpha=1) +
      geom_line(color="black") +
      coord_cartesian(xlim=c(left_edge-print_surround,left_edge+print_surround),
                      ylim=c(0,1)) +
      scale_color_brewer(palette = "Set1",name="Cluster") +theme_bw() +
      xlab("Basepair") + ylab("Fst") +
      guides(color=FALSE) +
      ggtitle(paste("Left edge = ",left_edge,sep=""))
    
    edge_plots[[(2 * i) - 1]] <- plot1
    
    plot2 <- fst %>%  
      mutate(window = (floor(pos/window_size)*window_size)+(window_size/2)) %>%
      mutate(Fst = case_when(Fst < 0 ~ 0,
                             TRUE ~ Fst)) %>%
      group_by(window) %>%
      summarize(window_fst = sum(FstNum)/sum(FstDenom)) %>%
      ggplot(.,aes(x=window,y=window_fst)) + 
      geom_point(data=fst,
                 aes(x=pos, y=Fst,color=as.factor(cluster)),alpha=1) +
      geom_line(color="black") +
      coord_cartesian(xlim=c(right_edge-print_surround,right_edge+print_surround),
                      ylim=c(0,1)) +
      scale_color_brewer(palette = "Set1",name="Cluster") +theme_bw() +
      xlab("Basepair") + ylab("Fst") +
      guides(color=FALSE) +
      ggtitle(paste("Right edge = ",right_edge,sep=""))
    
    edge_plots[[(2 * i)]] <- plot2
  }
  
  title <- paste(pull(inversions[n,1])," ",paste("Ha412HOChr",sprintf("%02d",pull(inversions[n,4])),sep="")," ",pull(inversions[n,6]),pull(inversions[n,5]),sep="")
  print(
    if (n_clusters == 1){
      grid.arrange( full_plot,edge_plots[[1]],edge_plots[[2]],
                    top = textGrob(title,gp=gpar(fontsize=20,font=2)),
                    layout_matrix = rbind(c(1,1),
                                          c(2,3)))
    }else if(n_clusters == 2){
      grid.arrange( full_plot,edge_plots[[1]],edge_plots[[2]],edge_plots[[3]],edge_plots[[4]],
                    top = textGrob(title,gp=gpar(fontsize=20,font=2)),
                    layout_matrix = rbind(c(1,1),
                                          c(2,3),
                                          c(4,5)))
    }else if(n_clusters == 3){
      grid.arrange( full_plot,edge_plots[[1]],edge_plots[[2]],edge_plots[[3]],edge_plots[[4]],
                    edge_plots[[5]],edge_plots[[6]],
                    top = textGrob(title,gp=gpar(fontsize=20,font=2)),
                    layout_matrix = rbind(c(1,1),
                                          c(2,3),
                                          c(4,5),
                                          c(6,7)))
    }else if (n_clusters == 4){
      grid.arrange( full_plot,edge_plots[[1]],edge_plots[[2]],edge_plots[[3]],edge_plots[[4]],
                    edge_plots[[5]],edge_plots[[6]],edge_plots[[7]],edge_plots[[8]],
                    top = textGrob(title,gp=gpar(fontsize=20,font=2)),
                    layout_matrix = rbind(c(1,1),
                                          c(2,3),
                                          c(4,5),
                                          c(6,7),
                                          c(8,9)))
    }
  )
}
dev.off()
write_tsv(fst_regions,"Ha412HO_inv.dec11.inversions.fstregions.v1.txt")

