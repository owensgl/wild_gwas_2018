#Inversion fst break point analysis 
#Tried stuff but nothing actually worked. Oct 16, 2018
library(tidyverse)
library(ggthemes)
library(depmixS4)

directory <- "/media/owens/Copper/wild_gwas_2018/petpet/"
fst_files <- list.files(directory, pattern="fst.txt.gz")  %>% .[24]


fst.data <- data_frame(filename = fst_files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = purrr::map(filename, ~ read_tsv(file.path(directory, .))) # a new data column
  ) %>% unnest(.) %>% mutate(filename = gsub("Petiolaris.tranche90.snp.petpet.90.bi.","",filename)) %>%
  mutate(filename = gsub(".fst.txt.gz","",filename))


fst.data %>% 
  mutate(window = floor(pos/window_size)*window_size,
         outlier = case_when(freq_dif > 0.75 ~ 1,
                             TRUE ~ 0)) %>% 
  group_by(chr, window, outlier) %>%
  tally() %>%
  spread(outlier,n, fill=0) %>%
  mutate(perc_outlier = `1` / (`0` + `1`)) %>%
  filter(perc_outlier > 0.1) -> tmp

for (row in 1:nrow(tmp)){
  current_chr <- tmp$chr[row]
  current_start <- tmp$window[row]
  tmp %>% filter(window != current_start) %>% filter(chr == current_chr)
}


fst.data %>% 
  #filter(chr == "HanXRQChr02") %>%
  #filter(pos > 100000000, pos < 102000000) %>%
  mutate(window = floor(pos/window_size)*window_size,
         outlier = case_when(freq_dif > 0.75 ~ 1,
                             TRUE ~ 0)) %>% 
  group_by(chr, window, outlier) %>%
  tally() %>%
  spread(outlier,n, fill=0) %>%
  mutate(perc_outlier = `1` / (`0` + `1`)) %>%
  ggplot(.,aes(x=window,y=perc_outlier)) + geom_point() +
  facet_wrap(~chr)
  
 window_size = 10000
tmp_window <- fst.data %>% 
  filter(chr == "HanXRQChr17") 
  #filter(pos > 5000000, pos < 12000000) %>%
  mutate(window = floor(pos/window_size)*window_size,
         outlier = case_when(freq_dif > 0.75 ~ 1,
                             TRUE ~ 0)) %>% 
  group_by(chr, window, outlier) %>%
  tally() %>%
  spread(outlier,n, fill=0) %>%
  mutate(perc_outlier = `1` / (`0` + `1`))
  
mod <- depmix(freq_dif ~ 1, family = poisson(), nstates = 2, data = tmp_window )
set.seed(1)
fm2 <- fit(mod, verbose = FALSE)


tmp_window <- cbind(tmp_window, posterior(fm2)) 

tmp_window %>%
  ggplot(.,aes(x=pos,y=freq_dif)) + geom_point() +
  geom_line(aes(x=pos,y=state))

data(speed)
# 2-state model on rt and corr from speed data set
# with Pacc as covariate on the transition matrix
# ntimes is used to specify the lengths of 3 separate series
mod1 <- depmix(list(rt~1,corr~1),data=speed,transition=~Pacc,nstates=2,
               family=list(gaussian(),multinomial("identity")),ntimes=c(168,134,137))
fb <- forwardbackward(mod1)
all.equal(-sum(log(fb$sca)),fb$logLike)

trstart=c(0.99,0.01,0.01,0.99)
transInit(~1, 2, pstart=trstart)


