

i <- 1
clustering_pcs <- tibble(chrom=character(),start=numeric(),end=numeric(),mid=numeric(),
                  betweenss=numeric())
for (i in 1:nrow(win.regions)){
  tmp <- as.data.frame(snp.pca[i,4:length(snp.pca[i,])])
  tibble::rownames_to_column(tmp) -> tmp
  colnames(tmp) <- c("x","value")
  tmp %>% separate(x, into = c('PC', 'name'), sep = 5) -> tmp
  tmp$PC <- str_sub(tmp$PC, 1, 4)
  tmp <- tmp %>% spread(PC, value) 

  tmp %>% ggplot(.,aes(x=PC_1,y=PC_2)) + geom_point()
  
  
  tmp %>% pull(PC_1) %>%  as.numeric(.) -> tmp_PC1
  output <- tibble(chrom=win.regions$chrom[i],
                   start=win.regions$start[i],
                   end=win.regions$end[i],
                   mid=win.regions$mid[i],
                   betweenss=kmeans(tmp_PC1, 3, centers=c(min(tmp_PC1),(min(tmp_PC1)+max(tmp_PC1))/2,max(tmp_PC1)))$betweenss)
  clustering_pcs <- rbind(clustering_pcs, output)
  print(paste("Loaded",i))
  
}

clustering_pcs %>% 
  ggplot(.,aes(x=mid,y=betweenss)) + geom_point() + facet_wrap(~chrom,scales="free_x")


