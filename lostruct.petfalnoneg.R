library(tidyverse)
library(lostruct)
library(Matrix)
library(colorspace)
library(RColorBrewer)
library(ggmap)
library(scatterpie)
library(SNPRelate)
library(gridExtra)
library(ggExtra)
library(grid)

data_directory <- "/media/owens/Copper/wild_gwas_2018/petfal/"
data_name <- "Petiolaris.tranche90.snp.petfalnoneg.90.bi"

labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T)
pop_loc <- read_tsv("pop_loc_allnum.txt")
pop_loc %>% rename(population = pop) %>% inner_join(.,labels) -> labels

#Prepare map data
usa <- map_data('state')
states <- map_data("state")
target_state <- map_data('state')
lat_range <- c(30, 40)
long_range <- c(-115,-100)

window_size <- 500
bcf.file <- paste(data_directory, "/", data_name,".bcf",sep="")
samples <- read_tsv(paste(data_directory, "/", data_name,".samplelist.txt",sep=""),col_names = F)
colnames(samples) <- c("sample")
sites <- vcf_positions(bcf.file)
win.fn.snp <- vcf_windower(bcf.file, size=window_size, type="snp", sites=sites) 
system.time( snp.pca <- eigen_windows(win.fn.snp,k=2, mc.cores=5) )
#Took 19677 time to run
system.time( pcdist <- pc_dist( snp.pca ) )

pcdist_na <- which(is.na(pcdist), TRUE)

na.inds <- is.na(pcdist[,1]) 
k_kept <- 20
mds <- cmdscale( pcdist[!na.inds,!na.inds], eig=TRUE, k=k_kept )
mds.coords <- mds$points
colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))
win.regions <- region(win.fn.snp)()
win.regions <- win.regions[!na.inds,]
win.regions %>% mutate(mid = (start + end) / 2) ->  win.regions
#Add the columns for all the MDS coordinates
for (k in 1:k_kept){
  str_pad(k, 2, pad = "0")
  
  name = paste("mds",str_pad(k, 2, pad = "0"),sep="")
  win.regions$tmp <- "NA"
  win.regions <- win.regions %>% rename(!!name := tmp) 
}

#Add the MDS coordinates to each window.
for (i in 1:k_kept){
  j = i + 4
  win.regions[,j] <- mds.coords[,i]
}
win.regions$n <- 1:nrow(win.regions)

#Make plots of all the MDS to visualize patterns
pdf(paste(data_name,".",window_size,".MDSplots.pdf",sep=""), height=16,width=25)
for (k in seq(1,k_kept,2)){
  
  name1 = paste("mds",str_pad(k, 2, pad = "0"),sep="")
  name2 = paste("mds",str_pad((k+1), 2, pad = "0"),sep="")
  print(
    win.regions %>% ggplot(.,aes_string(x=name1,y=name2)) + geom_point(aes(color=chrom),size=5) + theme_bw()
  )
}
win.regions %>%
  gather(., mds, value, colnames(win.regions)[5:(ncol(win.regions)-1)]) %>%
  ggplot(.,aes(value)) + geom_histogram() + facet_wrap(~mds,scales = "free")

win.regions %>%
  gather(., mds, value, colnames(win.regions)[5:(ncol(win.regions)-1)]) %>%
  ggplot(aes(x=mid,y=value)) + geom_point(alpha=0.5) +
  facet_grid(mds~chrom,scales="free") + theme_bw() 
dev.off()


saveRDS(win.regions, file = paste(data_name,".",window_size,".windows.rds",sep=""))


#Trying to select outliers for each MDS PC
mds_pcs <- colnames(win.regions)[5:(ncol(win.regions)-1)]
min_windows <- 4


mds_clustering <- tibble(mds_coord = character(),direction = character(),clust_pvalue = numeric(), outliers = numeric(),n1_outliers=numeric(),
                         high_cutoff = numeric(),lower_cutoff=numeric(),chr=character())
for (mds_chosen in mds_pcs){
  print(paste("Processing",mds_chosen))
  win.regions %>%
    mutate_(the_mds = mds_chosen ) %>%
    mutate(sd_mds = sd(the_mds)) %>%
    filter(the_mds > (sd_mds *4)) -> pos_windows
  
  win.regions %>%
    mutate_(the_mds = mds_chosen ) %>%
    summarize(sd_mds = sd(the_mds)) %>% pull() -> sd_mds
  mds_high_cutoff <- sd_mds * 4
  mds_low_cutoff <- sd_mds * 3

  n_permutations <- 1000
  if (nrow(pos_windows) >= min_windows){
    permutations <- matrix( nrow = n_permutations, ncol = 1)
    for (i in 1:n_permutations){
      max_1 <- win.regions %>% sample_n(nrow(pos_windows)) %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
        ungroup() %>% summarize(sum=sum(count)) %>% pull(sum)
      permutations[i,1] <- max_1
    }
    pos_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
      ungroup() %>% summarize(sum=sum(count)) %>% pull(sum) -> sampled_max_1
    pos_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>%
      pull(chrom) -> clustered_chr
    
    x <- as.tibble(permutations) %>%  filter(V1 >= sampled_max_1) %>% nrow() 
    x <- x+1
    pvalue <- x/(n_permutations+1)
    tmp <- tibble(mds_coord = as.character(mds_chosen),direction = as.character("pos"),clust_pvalue = as.numeric(pvalue), outliers=as.numeric(nrow(pos_windows)),
                  n1_outliers=as.numeric(sampled_max_1), high_cutoff = as.numeric(mds_high_cutoff),lower_cutoff=as.numeric(mds_low_cutoff),
                  chr=as.character(clustered_chr))
    mds_clustering <- rbind(mds_clustering, tmp)
  }else{
    tmp <- tibble(mds_coord = as.character(mds_chosen),direction = as.character("pos"),clust_pvalue = as.numeric(NA), outliers=as.numeric(nrow(pos_windows)),
                  n1_outliers=as.numeric(NA), high_cutoff = as.numeric(NA),lower_cutoff=as.numeric(NA),
                  chr=as.numeric(NA))
    mds_clustering <- rbind(mds_clustering, tmp)
    
  }
  
  
  win.regions %>%
    mutate_(the_mds = mds_chosen ) %>%
    mutate(sd_mds = sd(the_mds)) %>%
    filter(the_mds < -(sd_mds *4)) -> neg_windows
  
  if (nrow(neg_windows) >= min_windows){
    permutations <- matrix( nrow = n_permutations, ncol = 1)
    for (i in 1:n_permutations){
      max_1 <- win.regions %>% sample_n(nrow(neg_windows)) %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
        ungroup() %>% summarize(sum=sum(count)) %>% pull(sum)
      permutations[i,1] <- max_1
    }
    neg_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
      ungroup() %>% summarize(sum=sum(count)) %>% pull(sum) -> sampled_max_1
    neg_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>%
      pull(chrom) -> clustered_chr
    
    x <- as.tibble(permutations) %>%  filter(V1 >= sampled_max_1) %>% nrow() 
    x <- x+1
    pvalue <- x/(n_permutations+1)
    tmp <- tibble(mds_coord = as.character(mds_chosen),direction = as.character("neg"),clust_pvalue = as.numeric(pvalue), outliers=as.numeric(nrow(neg_windows)),
                  n1_outliers=as.numeric(sampled_max_1), high_cutoff = as.numeric(-mds_high_cutoff),lower_cutoff=as.numeric(-mds_low_cutoff),
                  chr=as.character(clustered_chr))
    mds_clustering <- rbind(mds_clustering, tmp)
  }else{
    tmp <- tibble(mds_coord = as.character(mds_chosen),direction = as.character("neg"),clust_pvalue = as.numeric(NA), outliers=as.numeric(nrow(neg_windows)),
                  n1_outliers=as.numeric(NA), high_cutoff = as.numeric(NA),lower_cutoff=as.numeric(NA),
                  chr=as.numeric(NA))
    mds_clustering <- rbind(mds_clustering, tmp)
    
  }
}
mds_clustering %>% filter(clust_pvalue < 0.01) -> sig_mds_clusters


outlier_windows <- tibble(chrom=character(),start=numeric(),end=numeric(),mid=numeric(),the_mds=numeric(),mds_coord=character(),outlier=character(),n=numeric())  
cluster_genotypes <- tibble(mds_coord=character(),name=character(),PC1=numeric(),genotype=character())

#For each mds outlier set that is chromosomally clustered, print a bunch of stuff about it.
pdf(paste(data_name,".",window_size,".mdsoutliers.pdf",sep=""),height=8,width=18)
for (i in 1:nrow(sig_mds_clusters)){
  coord <- pull(sig_mds_clusters[i,1])
  direction <- pull(sig_mds_clusters[i,2])
  high_cutoff <- pull(sig_mds_clusters[i,6])
  low_cutoff <- pull(sig_mds_clusters[i,7])
  cluster_chr <- pull(sig_mds_clusters[i,8])
  coord_direction <- paste(coord, "-",direction,sep="")
  
  
  if (direction == "pos"){
    windows <- win.regions %>%
      mutate_(the_mds = coord ) %>% 
      filter(the_mds > high_cutoff & chrom == cluster_chr) %>% pull(n)
    
    #Check if there are low outlier.
    if (win.regions %>%
        mutate_(the_mds = coord ) %>% 
        filter(the_mds < high_cutoff & the_mds > low_cutoff & chrom == cluster_chr) %>%
        nrow() > 0){
      genome_plot <- win.regions %>%
        mutate_(the_mds = coord ) %>% 
        mutate(outlier = case_when((the_mds > high_cutoff) & (chrom == cluster_chr) ~ "Higher",
                                   (the_mds > low_cutoff) & (chrom == cluster_chr) ~ "Lower",
                                   TRUE ~ "Non-outlier")) %>%
        ggplot(.,aes(x=mid/1000000,y=the_mds,color=outlier)) + geom_point() + facet_wrap(~chrom,scales= "free_x",nrow=1) +
        theme_bw() + scale_color_manual(breaks=c("Higher", "Lower", "Non-outlier"),values=c("#E41A1C","#377EB8","black")) + xlab("MB") + ylab(paste(coord))
    }else{
      genome_plot <- win.regions %>%
        mutate_(the_mds = coord ) %>% 
        mutate(outlier = case_when((the_mds > high_cutoff) & (chrom == cluster_chr) ~ "Higher",
                                   (the_mds > low_cutoff) & (chrom == cluster_chr) ~ "Lower",
                                   TRUE ~ "Non-outlier")) %>%
        ggplot(.,aes(x=mid/1000000,y=the_mds,color=outlier)) + geom_point() + facet_wrap(~chrom,scales= "free_x",nrow=1) +
        theme_bw() + scale_color_manual(values=c("#E41A1C","black")) + xlab("MB") + ylab(paste(coord))
    }
    
    
    
    outlier_windows <- rbind( win.regions %>%
                                mutate_(the_mds = coord ) %>% 
                                mutate(outlier = case_when((the_mds > high_cutoff) & (chrom == cluster_chr) ~ "Higher",
                                                           (the_mds > low_cutoff) & (chrom == cluster_chr) ~ "Lower",
                                                           TRUE ~ "Non-outlier"))  %>%
                                filter(outlier != "Non-outlier") %>%
                                select(chrom,start,end,mid,the_mds,outlier,n) %>%
                                mutate(mds_coord = coord_direction), outlier_windows)
    
  }else{
    windows <- win.regions %>%
      mutate_(the_mds = coord ) %>% 
      filter(the_mds < high_cutoff & chrom == cluster_chr) %>% pull(n)
    #Check if there are low outlier.
    if (win.regions %>%
        mutate_(the_mds = coord ) %>% 
        filter(the_mds > high_cutoff & the_mds < low_cutoff & chrom == cluster_chr) %>%
        nrow() > 0){
      
      genome_plot <- win.regions %>%
        mutate_(the_mds = coord ) %>% 
        mutate(outlier = case_when((the_mds < high_cutoff) & (chrom == cluster_chr) ~ "Higher",
                                   (the_mds < low_cutoff) & (chrom == cluster_chr) ~ "Lower",
                                   TRUE ~ "Non-outlier")) %>%
        ggplot(.,aes(x=mid/1000000,y=the_mds,color=outlier)) + geom_point() + facet_wrap(~chrom,scales= "free_x",nrow=1) +
        theme_bw() + scale_color_manual(values=c("#E41A1C","#377EB8","black")) + xlab("MB") + ylab(paste(coord))
    }else{
      genome_plot <- win.regions %>%
        mutate_(the_mds = coord ) %>% 
        mutate(outlier = case_when((the_mds < high_cutoff) & (chrom == cluster_chr) ~ "Higher",
                                   (the_mds < low_cutoff) & (chrom == cluster_chr) ~ "Lower",
                                   TRUE ~ "Non-outlier")) %>%
        ggplot(.,aes(x=mid/1000000,y=the_mds,color=outlier)) + geom_point() + facet_wrap(~chrom,scales= "free_x",nrow=1) +
        theme_bw() + scale_color_manual(values=c("#E41A1C","black")) + xlab("MB") + ylab(paste(coord))
    }
    outlier_windows <- rbind(win.regions %>%
                               mutate_(the_mds = coord ) %>% 
                               mutate(outlier = case_when((the_mds < high_cutoff) & (chrom == cluster_chr) ~ "Higher",
                                                          (the_mds < low_cutoff) & (chrom == cluster_chr) ~ "Lower",
                                                          TRUE ~ "Non-outlier")) %>%
                               filter(outlier != "Non-outlier") %>%
                               select(chrom,start,end,mid,the_mds,outlier,n) %>%
                               mutate(mds_coord = coord_direction), outlier_windows)
    
    
  }
  out <- cov_pca(win.fn.snp(windows),k=2)
  out <- matrix(out[4:length(out)],ncol=nrow(samples),byrow=T) %>% as.tibble() 
  colnames(out) <- pull(samples)
  
  out <- as_tibble(cbind(nms = names(out), t(out))) %>% rename(name=nms,PC1=V1,PC2=V2) 
  
  
  
  pca_plot <- out %>%
    mutate(PC1 = as.double(PC1), PC2 = as.double(PC2)) %>%
    ggplot(.,aes(x=PC1,y=PC2)) + geom_point() + theme_bw()
  
  out %>% pull(PC1) %>%  as.numeric(.) -> tmp_PC1
  
  kmeans_cluster <-kmeans(tmp_PC1, 3, centers=c(min(tmp_PC1),(min(tmp_PC1)+max(tmp_PC1))/2,max(tmp_PC1)))
  
  out$cluster <- kmeans_cluster$cluster - 1 
  out$cluster <- as.character(out$cluster)
  out$mds_coord <- paste(coord, direction,sep="_")
  
  genotype.out <- out %>% select(mds_coord,name,PC1,cluster) %>%
    rename(genotype = cluster) 
  cluster_genotypes <- rbind(cluster_genotypes, genotype.out)
  
  hist_plot <- out %>%
    mutate(PC1 = as.double(PC1), PC2 = as.double(PC2)) %>%
    ggplot(.,aes(PC1)) + geom_histogram(aes(fill=as.character(cluster))) +
    theme_bw() + scale_fill_brewer(palette = "Set1",name="Cluster")
  
  
  win.fn.snp(windows) %>% as.tibble() -> snps
  colnames(snps) <- pull(samples)
  
  snps %>% gather("name","genotype",1:ncol(snps)) %>%group_by(name, genotype) %>%
    summarize(count=n()) %>%
    spread(genotype, count) %>%
    summarize(het=`1`/(`0` + `1` + `2`)) -> heterozygosity
  
  het_plot <- inner_join(out, heterozygosity) %>% 
    ggplot(.,aes(x=as.character(cluster),y=het,fill=as.character(cluster))) + 
    geom_boxplot() + scale_fill_brewer(palette = "Set1",name="Cluster") + theme_bw() + xlab("Cluster") + ylab("Heterozygosity")
  
  map_plot <- ggplot(target_state, aes(long, lat)) +
    geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
    coord_quickmap() +
    geom_scatterpie(data=inner_join(out, labels) %>% 
                      group_by(population, lat, long, cluster) %>% 
                      tally() %>%
                      spread(., cluster, n,fill=0),
                    aes(x=long, y=lat, r=0.4), 
                    cols=c("0","1","2"), color=NA, alpha=.8) +
    scale_fill_brewer(name="Cluster",palette = "Set1") +theme_bw() +
    xlab("Longitude") + ylab("Latitude") +
    scale_x_continuous(limits = long_range, expand = c(0, 0)) +
    scale_y_continuous(limits = lat_range, expand = c(0, 0))
  
  print(
    grid.arrange(
      pca_plot, hist_plot, het_plot ,map_plot, genome_plot,
      widths = c(1, 1, 1, 1),
      layout_matrix = rbind(c(1, 2, 3, 4),
                            c(5, 5, 5, 5)),
      top = textGrob(paste(coord, direction),gp=gpar(fontsize=20,font=1))
    )
  )
  
  
}
dev.off()


#Save cluster text output
write_tsv(outlier_windows, paste(data_name,".",window_size,"mds_cluster_windows.txt",sep=""))
write_tsv(cluster_genotypes, paste(data_name,".",window_size,"mds_cluster_genotyped.txt",sep=""))
