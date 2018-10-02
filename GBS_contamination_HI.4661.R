#Plotting contamination of HI.4661
library(tidyverse)
contamination <- read.delim("HI.4661_contamination.txt",header=T)


pdf("HI.4661_contamination.v1.pdf")
contamination %>% filter(GBS_lane != 1) %>%
  group_by(file, same_lane, Lane) %>%
  summarize(mean_contam = mean(percent_barcoded)) %>%
  ggplot(.) +
  geom_jitter(aes(x=as.factor(same_lane),y=mean_contam,color=as.factor(Lane)),
              width=0.2,size=3) +
  scale_x_discrete(labels = c('Different_lane','Same_lane')) + 
  scale_color_brewer(palette = "Set1",name="Lane") +
  ylab("Percent contaminated reads") +
  xlab("Potential source of contamination") +
  ggtitle("Percent GBS contamination in HI.4661") +
  theme_bw()
dev.off()
  
  
