#FST between 
library(tidyverse)
library(ggthemes)



directory <- "/media/owens/Copper/wild_gwas_2018/petpet/"
fst_files <- list.files(directory, pattern="fst.txt.gz")


fst.data <- data_frame(filename = fst_files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = purrr::map(filename, ~ read_tsv(file.path(directory, .))) # a new data column
  ) %>% unnest(.) %>% mutate(filename = gsub("Petiolaris.tranche90.snp.petpet.90.bi.","",filename)) %>%
  mutate(filename = gsub(".fst.txt.gz","",filename))


fst.data %>% head()


w_size <- 1000000
pdf("Petiolaris.tranche90.snp.petpet.90.bi.fst.pdf")
for (mds in unique(fst.data$filename)){
  print(
    fst.data %>%
      filter(filename == mds) %>%
      mutate(window = floor(pos/w_size)*w_size) %>%
      group_by(chr,window,filename) %>%
      summarize(w_fst = sum(FstNum,na.rm=T)/sum(FstDenom,na.rm=T)) %>%
      ggplot(.,aes(x=window/1000000,y=w_fst)) +
      geom_line()+
      facet_wrap(~chr,scales = "free_x") +
      theme_few() +
      xlab("MB") +
      ylab(paste("Fst Window = ",w_size/1000000, " MB",sep="")) +
      ggtitle(paste(mds, " PetPet",sep=""))
  )
}
dev.off()




