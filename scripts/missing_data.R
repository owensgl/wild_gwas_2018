#Visualizing missing data example
library(tidyverse)

missing <- read_tsv("texanus_chr2_head_missingdata.txt")

labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T) %>%
  rename(sample = name)
read_depth <- read_delim("read_depths_summary.txt",delim=" ")

pdf("texanus_chr2_head_missingdata.pdf")
inner_join(missing,labels) %>%
  ggplot(.,aes(y=percent_missing,x=species)) + geom_boxplot() +
  theme_bw() 
inner_join(missing, read_depth) %>%
  ggplot(.,aes(x=kernel_max_depth,y=percent_missing,color=species)) + geom_point() +
  theme_bw()
dev.off()

         