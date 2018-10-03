#FST between 
library(tidyverse)
library(ggthemes)

####PETPET

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

#####PETFAL



directory <- "/media/owens/Copper/wild_gwas_2018/petfal/"
fst_files <- list.files(directory, pattern="fst.txt.gz")


fst.data <- data_frame(filename = fst_files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = purrr::map(filename, ~ read_tsv(file.path(directory, .))) # a new data column
  ) %>% unnest(.) %>% mutate(filename = gsub("Petiolaris.tranche90.snp.petfal.90.bi.","",filename)) %>%
  mutate(filename = gsub(".fst.txt.gz","",filename))


fst.data %>% head()


w_size <- 1000000
pdf("Petiolaris.tranche90.snp.petfal.90.bi.fst.pdf")
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
      ggtitle(paste(mds, " PetFal",sep=""))
  )
}
dev.off()

####Anomalus


directory <- "/media/owens/Copper/wild_gwas_2018/anomalus/"
fst_files <- list.files(directory, pattern="fst.txt.gz")


fst.data <- data_frame(filename = fst_files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = purrr::map(filename, ~ read_tsv(file.path(directory, .))) # a new data column
  ) %>% unnest(.) %>% mutate(filename = gsub("Anomalus.tranche90.snp.gwas.90.bi.","",filename)) %>%
  mutate(filename = gsub(".fst.txt.gz","",filename))


fst.data %>% head()


w_size <- 1000000
pdf("Anomalus.tranche90.snp.gwas.90.bi.fst.pdf")
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
      ggtitle(paste(mds, " Anomalus",sep=""))
  )
}
dev.off()


####Niveus


directory <- "/media/owens/Copper/wild_gwas_2018/niveus/"
fst_files <- list.files(directory, pattern="fst.txt.gz")


fst.data <- data_frame(filename = fst_files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = purrr::map(filename, ~ read_tsv(file.path(directory, .))) # a new data column
  ) %>% unnest(.) %>% mutate(filename = gsub("Niveus.tranche90.snp.gwas.90.bi.","",filename)) %>%
  mutate(filename = gsub(".fst.txt.gz","",filename))


fst.data %>% head()


w_size <- 1000000
pdf("Niveus.tranche90.snp.gwas.90.bi.fst.pdf")
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
      ggtitle(paste(mds, " Niveus",sep=""))
  )
}
dev.off()

###Argophyllus


directory <- "/media/owens/Copper/wild_gwas_2018/argophyllus/"
fst_files <- list.files(directory, pattern="fst.txt.gz")


fst.data <- data_frame(filename = fst_files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = purrr::map(filename, ~ read_tsv(file.path(directory, .))) # a new data column
  ) %>% unnest(.) %>% mutate(filename = gsub("Argophyllus.tranche90.snp.gwas.90.bi.","",filename)) %>%
  mutate(filename = gsub(".fst.txt.gz","",filename))


fst.data %>% head()


w_size <- 1000000
pdf("Argophyllus.tranche90.snp.gwas.90.bi.fst.pdf")
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
      ggtitle(paste(mds, " Argophyllus",sep=""))
  )
}
dev.off()

###Annuus


directory <- "/media/owens/Copper/wild_gwas_2018/annuus/"
fst_files <- list.files(directory, pattern="fst.txt.gz")


fst.data <- data_frame(filename = fst_files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = purrr::map(filename, ~ read_tsv(file.path(directory, .))) # a new data column
  ) %>% unnest(.) %>% mutate(filename = gsub("Annuus.tranche90.snp.env.90.bi.","",filename)) %>%
  mutate(filename = gsub(".fst.txt.gz","",filename))


fst.data %>% head()


w_size <- 1000000
pdf("Annnuus.tranche90.snp.env.90.bi.fst.pdf")
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
      ggtitle(paste(mds, " Annuus",sep=""))
  )
}
dev.off()

