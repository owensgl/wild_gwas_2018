library(LDna)
library(tidyverse)
library(SNPRelate)
library(forcats)
library(gridExtra)
library(grid)
genome_prefix <- "Ha412HOChr"

data_directory <- "/media/owens/Copper/wild_gwas_2018/annuus/"
data_name <- "Annuus.tranche90.snp.env.90.bi.remappedHa412HO"
ld_prefix <- ".maf5.thin5k"
samples <- read_tsv(paste(data_directory,data_name,".samplelist.txt",sep=""),col_names = c("sample"))
bcf.file <- paste(data_directory,data_name,".bcf",sep="")

clusters_tibble <- tibble(chr= character() ,pos=numeric(),cluster=character(),
                          type=character(),start_pos=numeric(),end_pos=numeric(),size=numeric())
summary_tibble <- tibble(Type=character(),Merge.at=numeric(),
                         nLoci=numeric(),nE=numeric(),lambda=numeric(),
                         Median.LD=numeric(),MAD.LD=numeric(),cluster=character())
pdf(paste(data_name,ld_prefix,".ldna.pdf",sep=""))
for (chr_n in sprintf("%02d", seq(17))){
  ld <- read_tsv(paste(data_directory,"/",data_name,ld_prefix,".Chr",chr_n,".geno.ld.gz",sep=""))
  
  ld %>% select(POS1,POS2) %>% gather(all_pos) %>% select(value) %>%distinct() %>%
    mutate(POS1=value,POS2=value,`R^2`="NA") %>%select(-value) ->self_values
  
  
  
  ld %>%select(POS1,POS2,`R^2`) %>% rbind(self_values) %>% 
    arrange(POS1,POS2) %>% spread(.,POS1,`R^2`) %>% select(-POS2) %>% data.matrix()-> ld.matrix
  
  LDnaRaw(ld.matrix) -> ldna
  clusters1 <- extractClusters(ldna, min.edges = 20,phi=10, rm.COCs = T)
  print(clusters1)
  
  summary <- summaryLDna(ldna, clusters1, ld.matrix)
  summary$cluster <- row.names(summary)
  summary_tibble <- rbind(summary_tibble,summary)
  chr_tibble <- tibble(chr= character() ,pos=character(),cluster=character(),
                            type=character(),start_pos=character(),end_pos=character(),size=character())
  #Now turn cluster list into tibble
  for (n in 1:length(clusters1)){
    name<-names(clusters1)[n]
    if (is.null(name)){
      name <- "NA"
    }
    #type <- summary %>% filter(cluster == name) %>% pull(Type) %>% as.character()
    type <- "SOC"
    tmp_tibble <- tibble(chr= paste(genome_prefix,chr_n,sep="") ,pos=clusters1[[n]],cluster=name)
    tmp_tibble <- tmp_tibble %>% group_by(cluster) %>%
      mutate(start_pos = min(as.numeric(pos)),end_pos = max(as.numeric(pos)),
             size= as.numeric(end_pos) -  as.numeric(start_pos),
             type=type) %>% ungroup()
    chr_tibble <- rbind(chr_tibble,tmp_tibble)
  }
  sites <- unique(c(ld$POS1,ld$POS2))
  print(
  ggplot(chr_tibble,aes(x=as.numeric(pos)/1000000,y=cluster,color=type)) + geom_point() + theme_bw() +
    theme(legend.position="none") + xlab("MB") + ggtitle(paste(genome_prefix,chr_n,sep=""))
  )
  clusters_tibble <- rbind(clusters_tibble,chr_tibble)
  
  # ld %>% filter(POS1 %in% chr_tibble$pos[which(chr_tibble$cluster == "635_0.72")] &
  #               POS2 %in% chr_tibble$pos[which(chr_tibble$cluster == "635_0.72")]) %>%
  #   ggplot(.,aes(x=as.factor(POS1),y=as.factor(POS2),fill=`R^2`)) + 
  #   geom_tile() + scale_fill_viridis_c(limits=c(0,1))
}
dev.off()

clusters_tibble %>% group_by(cluster) %>%
  mutate(start_pos = min(as.numeric(pos)),end_pos = max(as.numeric(pos)),size= as.numeric(end_pos) -  as.numeric(start_pos),type="main") %>% 
  ungroup() -> clusters_tibble

write_tsv(clusters_tibble, paste(data_directory,"/",data_name,ld_prefix,".ldnasites.txt",sep=""))
write_tsv(summary_tibble, paste(data_directory,"/",data_name,ld_prefix,".ldnasummary.txt",sep=""))




clusters_tibble <- read_tsv( paste(data_directory,"/",data_name,ld_prefix,".ldnasites.txt",sep=""))
full_clusters_tibble <- clusters_tibble


#Trim the lists of markers by pulling off isolated SNPs too far from others within the same cluster
# max_distance_between_outliers <- 5000000
# full_clusters_tibble %>%
#   group_by(chr, cluster) %>%
#   arrange(pos) %>%
#   mutate(ahead_n = pos - lag(pos),behind_n = abs(pos - lead(pos))) %>%
#   mutate(min_dist = pmin(ahead_n, behind_n,na.rm=T)) %>%
#   filter(min_dist < max_distance_between_outliers ) %>%
#   select(-ahead_n, -behind_n, -min_dist) %>%
#   group_by(chr,cluster) %>% 
#   mutate(start_pos = min(as.numeric(pos)),
#          end_pos = max(as.numeric(pos)),
#          size= as.numeric(end_pos) -  as.numeric(start_pos)) %>% 
#   ungroup() -> full_clusters_tibble



#Run PCA on 
full_clusters_tibble %>% 
  select(cluster,chr, start_pos, end_pos,size) %>% unique() %>% 
  arrange(desc(size)) %>% mutate(betweenSS = "NA")-> clusters

het_script <- "/home/owens/bin/pop_gen/vcf2het.pl"
#Calculate the kmeans clustering 
for (n in 1:nrow(clusters)){
  selected_cluster <- clusters$cluster[n]
  selected_chr <- clusters$chr[n]
  full_clusters_tibble %>%
    filter(cluster == selected_cluster,
           chr == selected_chr) %>%
    select(chr, pos) -> tmp.sites
  
  write_tsv(tmp.sites,"tmp.sites.txt",col_names=F)
  
  system(paste("bcftools view -R tmp.sites.txt ",bcf.file," > tmp.vcf"))
  
  system(paste("cat tmp.vcf | perl ",het_script," > tmp.het.txt",sep=""))
  het <- read_tsv("tmp.het.txt")
  snpgdsVCF2GDS("tmp.vcf", "tmp2.gds", method="biallelic.only",ignore.chr.prefix="Ha412HOChr")
  genofile <- snpgdsOpen("tmp2.gds")
  pca <- snpgdsPCA(genofile, num.thread=2)
  tab <- data.frame(sample = pca$sample.id,
                    EV1 = pca$eigenvect[,1],    # the first eigenvector
                    EV2 = pca$eigenvect[,2],    # the second eigenvector
                    stringsAsFactors = FALSE)
  snpgdsClose(genofile)
  
  try_3_clusters <-try(kmeans(tab[,2], 3, centers=c(min(tab[,2]),(min(tab[,2])+max(tab[,2]))/2,max(tab[,2]))))
  
  
  if("try-error" %in% class(try_3_clusters)){
    kmeans_cluster <-kmeans(tab[,2], 2, centers=c(min(tab[,2]),max(tab[,2])))
    out.normal <- out
    out.rotation <- out
  }else{
    
    kmeans_cluster <-kmeans(tab[,2], 3, centers=c(min(tab[,2]),(min(tab[,2])+max(tab[,2]))/2,max(tab[,2])))
    
  }
  
  tab$cluster <- kmeans_cluster$cluster - 1 
  
  betweenss <- kmeans_cluster$betweenss
  clusters$betweenSS[n] <- betweenss
  
  
}


#Choose filtering for examining clusters
filtered_clusters <- clusters %>% filter(size > 1000000,betweenSS > 0.95) 

filtered_clusters %>%
  mutate(chr = gsub("Ha412TMP","",chr)) %>%
  ggplot(.,aes(x=start_pos/1000000,xend=end_pos/1000000,y=cluster,yend=cluster)) +geom_segment() +
  facet_grid(chr~.,scales="free_y")


pdf(paste(data_name,ld_prefix,".LDna.clusteroutliers.pdf",sep=""),height=8,width=12)
for (n in 1:nrow(filtered_clusters)){
  selected_cluster <- filtered_clusters$cluster[n]
  selected_chr <- filtered_clusters$chr[n]
  selected_start <- round(filtered_clusters$start_pos[n]/1000000,3)
  selected_end <- round(filtered_clusters$end_pos[n]/1000000,3)
  selected_size <- round(filtered_clusters$size[n]/1000000,3)
  full_clusters_tibble %>%
    filter(cluster == selected_cluster,
           chr == selected_chr) %>%
    filter(type == "main") %>%
    select(chr, pos) -> tmp.sites
  
  write_tsv(tmp.sites,"tmp.sites.txt",col_names=F)
  
  system(paste("bcftools view -R tmp.sites.txt ",bcf.file," > tmp.vcf"))
  
  system(paste("cat tmp.vcf | perl ",het_script," > tmp.het.txt",sep=""))
  het <- read_tsv("tmp.het.txt")
  snpgdsVCF2GDS("tmp.vcf", "tmp3.gds", method="biallelic.only",ignore.chr.prefix="Ha412HOChr")
  genofile <- snpgdsOpen("tmp3.gds")
  pca <- snpgdsPCA(genofile, num.thread=2)
  tab <- data.frame(sample = pca$sample.id,
                    EV1 = pca$eigenvect[,1],    # the first eigenvector
                    EV2 = pca$eigenvect[,2],    # the second eigenvector
                    stringsAsFactors = FALSE)
  snpgdsClose(genofile)
  
  try_3_clusters <-try(kmeans(tab[,2], 3, centers=c(min(tab[,2]),(min(tab[,2])+max(tab[,2]))/2,max(tab[,2]))))
  
  
  if("try-error" %in% class(try_3_clusters)){
    kmeans_cluster <-kmeans(tab[,2], 2, centers=c(min(tab[,2]),max(tab[,2])))
    out.normal <- out
    out.rotation <- out
  }else{
    
    kmeans_cluster <-kmeans(tab[,2], 3, centers=c(min(tab[,2]),(min(tab[,2])+max(tab[,2]))/2,max(tab[,2])))
    
  }
  
  tab$cluster <- kmeans_cluster$cluster - 1 
  
  betweenss <- kmeans_cluster$betweenss
  
  
  plot.pca <- tab %>% 
    inner_join(.,het) %>%
    ggplot(.,aes(x=EV1,y=EV2)) + geom_point(aes(color=percent_het,shape=as.factor(cluster))) +
    theme_bw() + ggtitle(paste(selected_cluster, " BetweenSS = ",betweenss, sep="")) + 
    scale_color_viridis_c(name="Heterozygosity") +
    scale_shape(name="Kmeans\ncluster") +
    theme(legend.position="bottom")
  
  system(paste("bcftools query -H -f '%END [ %GT]\n' -R tmp.sites.txt ",bcf.file,
               '| sed s/\\#\\ //g |  sed s/\\:GT//g | sed s/END/pos/g > tmp.geno.txt',sep="")) 
  
  read_delim("tmp.geno.txt",delim=" ",col_names = c("pos","blank", as.vector(samples$sample)),skip=1) %>%
    select(-blank) %>% mutate_if(., 
                                 is.character, 
                                 str_replace_all, pattern = '0/0', replacement = "0") %>%
    mutate_if(., 
              is.character, 
              str_replace_all, pattern = '1/1', replacement = "2") %>%
    mutate_if(., 
              is.character, 
              str_replace_all, pattern = '0/1', replacement = "1") %>%
    mutate_if(., 
              is.character, 
              str_replace_all, pattern = './.', replacement = "NA") %>%
    group_by(pos) %>%gather("name","genotype",2:(ncol(.))) %>%
    rename(sample=name) -> snps
  
  #Make a plot of the SNPs themselves.
  plot.snps <- snps %>%  
    inner_join(.,tab) %>%
    filter(genotype != "NA") %>%
    arrange(pos) %>%
    ggplot(.,aes(x=as.factor(pos),y=fct_reorder(sample,EV1), fill=as.factor(genotype))) + geom_tile() +
    scale_fill_brewer(palette = "Set1",name="Genotype") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ylab("Sample") + xlab("Site")
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #scale_x_discrete(breaks = levels(as.factor(as.character(snps$pos)))[seq(1, length(levels(as.factor(as.character(snps$pos)))),2)]) 
  
  #Make a plot of SNP positions relative to total length. 
  
  plot.loci <- filtered_clusters %>% 
    filter(cluster == selected_cluster,
           chr == selected_chr) %>%
    ggplot(.) +
    geom_segment(aes(x=start_pos/1000000,xend=end_pos/1000000,y=1,yend=1)) +
    geom_point(data=tmp.sites,aes(x=as.numeric(pos)/1000000,y=1)) +
    theme_bw() + theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank()) +
    xlab("MB")  +
    ggtitle(paste(selected_chr))
  
  print(
    grid.arrange( plot.pca,plot.loci,plot.snps,
                  #nrow=3,
                  top = textGrob(paste(selected_chr,
                                              ":",
                                              selected_start,
                                              "-",
                                              selected_end,
                                              "MB Size=",
                                              selected_size,
                                              "MB",sep=""),gp=gpar(fontsize=20,font=2)),
                  layout_matrix = rbind(c(1),
                                        c(1),
                                        c(2),
                                        c(3),
                                        c(3)))
  )
  
  #Look at recombination rate in largest cluster
  
  most_common_cluster <- tab %>% group_by(cluster) %>%
    summarize(sum=n()) %>%
    filter(cluster != 1) %>%
    arrange(desc(sum)) %>% head(1) %>% pull(cluster) 
  
  write_tsv(tab %>% filter(cluster == most_common_cluster) %>% select(sample),
            "tmp.samplelist.txt",col_names = FALSE)
  
  ld <- read_tsv(paste(data_directory,"tmp.geno.ld",sep=""))
  ld %>%
    ggplot(.,aes(x=as.factor(POS1),y=as.factor(POS2),fill=`R^2`)) + geom_tile() +
    scale_fill_viridis_c()
  window_size <- 100000
  ld %>% mutate(win1 = floor(POS1/window_size)*window_size, 
                win2 = floor(POS2/window_size)*window_size) %>%
    group_by(win1,win2) %>%
    summarize(mean_ld = max(`R^2`)) %>%
    ungroup() %>%
    mutate(min_win = min(win1),max_win = max(win2)) %>%
    filter(win1 == min_win,win2==max_win)  
    
  
}
dev.off()


