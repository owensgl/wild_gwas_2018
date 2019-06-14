#This is for visualizing some 10X data
library(tidyverse)
library(forcats)
library(gganimate)
library(gifski)
library(transformr)
directory <- "/media/owens/Fuchs/10x_ann1238/"

file <- "ann1238.Ha412HOChr17.neg2.linkedreads.txt"

data <- read_tsv(paste(directory,file,sep=""),col_names = c("chr","pos","barcode"))

window_size <- 1000000
data %>%
  mutate(window = floor(pos/window_size)*window_size) %>%
  group_by(barcode) %>%
  summarize(start=min(window),end=max(window)) %>%
  group_by(start, end) %>%
  tally() %>% 
  mutate(distance = end - start) %>%
  filter(distance > 2*window_size) %>% 
  filter(n > 10) %>% View()
  ggplot(.,aes(x=start,y=end,color=n),alpha=0.7) + geom_point(aes(size=n)) +
  scale_color_viridis_c() + theme_bw()

  
#Chr01.neg1 Candidates
  
  window_size <- 100000
  full_data <- tibble(chr=character(),pos=numeric(),targeted_zone=numeric())
  for (loc in seq(4900000, 10000000,window_size)){
    
    candidates <- c(loc)
    data %>%
      mutate(window = floor(pos/window_size)*window_size) %>%
      filter(window %in% candidates) %>% pull(barcode) %>% unique() -> candidate_barcodes
    if (length(candidate_barcodes) == 0){next}
    tmp_data <- data %>% 
      filter(barcode %in% candidate_barcodes) %>%
      mutate(targeted_zone = loc) %>%
      select(-barcode)
    full_data <- rbind(full_data,tmp_data)
  }
  pdf("ann1238.Ha412HOChr01.neg1.linkedreads.density.pdf",height=25,width=25)
  full_data %>%
    ggplot(.,aes(pos/1000000)) + geom_density() +
    #geom_vline(xintercept = (((targeted_zone+(window_size/2))/1000000)), linetype="dotted") +
    coord_cartesian(xlim=c(4, 10)) +
    facet_wrap(.~targeted_zone,scales="free_y") +
    xlab("MB") + ylab("Read density") +
    ggtitle("ANN1238 10X, Ha412HOChr01")
  dev.off()
  
  
  
  
#Candidates for Chr05.neg1
candidates <- c(148200000,148210000)
candidates <- c(148200000)
data %>%
  mutate(window = floor(pos/window_size)*window_size) %>%
  filter(window %in% candidates) %>% pull(barcode) %>% unique() -> candidate_barcodes

window_size <- 100000
full_data <- tibble(chr=character(),pos=numeric(),targeted_zone=numeric())
for (loc in seq(148000000, 177000000,window_size)){
  
  candidates <- c(loc)
  data %>%
    mutate(window = floor(pos/window_size)*window_size) %>%
    filter(window %in% candidates) %>% pull(barcode) %>% unique() -> candidate_barcodes
  tmp_data <- data %>% 
    filter(barcode %in% candidate_barcodes) %>%
    mutate(targeted_zone = loc) %>%
    select(-barcode)
  full_data <- rbind(full_data,tmp_data)
}
full_data %>%
  ggplot(.,aes(pos/1000000)) + geom_density() +
  #geom_vline(xintercept = (((targeted_zone+(window_size/2))/1000000)), linetype="dotted") +
  coord_cartesian(xlim=c(140, 190)) +
  facet_wrap(.~targeted_zone,scales="free_y")



####Chr17.neg2
window_size <- 100000
full_data <- tibble(chr=character(),pos=numeric(),targeted_zone=numeric())
for (loc in seq(180000000, 200000000,window_size)){
  
  candidates <- c(loc)
  data %>%
    mutate(window = floor(pos/window_size)*window_size) %>%
    filter(window %in% candidates) %>% pull(barcode) %>% unique() -> candidate_barcodes
  if (length(candidate_barcodes) == 0){next}
  tmp_data <- data %>% 
    filter(barcode %in% candidate_barcodes) %>%
    mutate(targeted_zone = loc) %>%
    select(-barcode)
  full_data <- rbind(full_data,tmp_data)
}
pdf("ann1238.Ha412HOChr17.neg2.linkedreads.density.pdf",height=25,width=25)
full_data %>%
  ggplot(.,aes(pos/1000000)) + geom_density() +
  #geom_vline(xintercept = (((targeted_zone+(window_size/2))/1000000)), linetype="dotted") +
  coord_cartesian(xlim=c(180, 200)) +
  facet_wrap(.~targeted_zone,scales="free_y") +
  xlab("MB") + ylab("Read density") +
  ggtitle("ANN1238 10X, Ha412HOChr07")
dev.off()