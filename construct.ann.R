library(conStruct)
library(tidyverse)
library(fields)
##HELP SECTIONS
# vignette(topic="format-data",package="conStruct")
# # how to run a conStruct analysis
# vignette(topic="run-conStruct",package="conStruct")
# 
# # how to visualize the output of a conStruct model
# vignette(topic="visualize-results",package="conStruct")
# 
# # how to compare and select between different conStruct models
# vignette(topic="model-comparison",package="conStruct")

ann.genotypes <- t(as.matrix(read.delim("/media/owens/Copper/wild_gwas_2018/annuus/Annuus.tranche90.snp.env.90.bi.0p01subset.flipconstruct.txt",
                     row.names = 1,header=T,sep=" ")))


pop_loc <- read_tsv("pop_loc_allnum.txt")
pop_loc %>% filter(taxon == "annuus") %>% select(long,lat) %>% as.matrix() -> ann.geography


rdist.earth(ann.geography, ann.geography) -> ann.distance




xvalid <- x.validation(train.prop = 0.9,
                       n.reps = 5, 
                       K = 1:5, 
                       freqs = ann.genotypes, 
                       geoDist = ann.distance, 
                       coords = ann.geography,
                       prefix = "Annuus.tranche90.snp.env.90.bi.0p01subset",
                       n.iter = 1000, 
                       make.figs = TRUE, 
                       save.files = FALSE,
                       parallel = TRUE,
                       n.nodes = 5)
