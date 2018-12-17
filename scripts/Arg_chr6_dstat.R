library(tidyverse)
library(gridExtra)
dstat <- read_tsv("/home/owens/working/texanus/WGS/arg_chr6_test/Texanus.tranche90.snp.gwas.90.bi.remappedHa412.chr6.dstats.txt",
                  skip=1,col_names = c("chr","pos","AAAA", "AABA", "ABAA", "ABBA", "BAAA", "BABA", "BBAA", "BBBA"))
pdf("Texanus.tranche90.snp.gwas.90.bi.remappedHa412.chr6.dstats.pdf")

plot1 <- dstat %>%
  mutate(inversion= case_when(pos > 30000000 ~ "normal_region",
                              pos <= 30000000 ~ "chr6_inversion")) %>%
  group_by(inversion) %>%
  summarize(Dstat = (sum(ABBA)-sum(BABA))/(sum(ABBA)+sum(BABA)),
            baba = sum(BABA), abba = sum(ABBA)) %>%
  gather(.,pattern,count, abba:baba) %>%
  ggplot(.,aes(y=count,fill=pattern,x=inversion)) + geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set1") + theme_bw()

plot2 <- dstat %>%
  ggplot(.) + geom_bar(aes(x=pos/1000000,y=ABBA),color="red",stat="identity") +
  geom_bar(aes(x=pos/1000000,y=-BABA),color="blue",stat="identity") + theme_bw() +
  ylab("red=ABBA blue=BABA") + xlab("MB")

grid.arrange(plot1, plot2, nrow = 2)

dev.off()
