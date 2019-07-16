#Plotting DTB GWAS peak example
library(tidyverse)
library(viridis)
library(ggbeeswarm)



inv_associations <- read_tsv(paste0("/media/owens/Copper/wild_gwas/gwas/","annuus","/results/","Ha412HO_inv.v3.pcasites.","annuus",
                                    ".shuffled.","Days_to_budding",".ps.gz"),
                             col_names=c("id","beta","sd","pvalue")) %>%
  mutate(variable = "Days_to_budding",
         logp = abs(log10(pvalue)))


#Load up inversion locations
inv_locations <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt") %>%
  filter(chr == "Ha412HOChr13",spe=="annuus",mds=="pos1")

inv_associations$start <- inv_locations$start[1]
inv_associations$end <- inv_locations$end[1]
inv_associations$type <- "No SVs"
min_freq <- 0.03
freqs <- read_tsv(paste0("/media/owens/Copper/wild_gwas/gwas/","annuus","/","Annuus",".tranche90.snp.","gwas",".90.bi.remappedHa412HO.beagle.freq"),
                  col_names=c("chr","pos","freq")) %>%
  filter(freq > min_freq, freq < (1-min_freq)) %>%
  select(-freq)

gwas_noSV <- read_tsv(paste0("/media/owens/Copper/wild_gwas/gwas/","annuus","/results/","Annuus",".tranche90.snp.","gwas",
                        ".90.bi.remappedHa412HO.","Days_to_budding",".beagle.fullgenome.v3noinv.ldfilter.ps.gz"),
                 col_names=c("id","beta","sd","pvalue")) %>%
  separate(id,c("chr","pos"),"_") %>%
  mutate(pos=as.numeric(pos)) %>%
  semi_join(.,freqs) %>%
  filter(chr == "13") %>%
  mutate(type = "No SVs",logp=abs(log10(pvalue)))

gwas_SV <- read_tsv(paste0("/media/owens/Copper/wild_gwas/gwas/","annuus","/results/","Annuus",".tranche90.snp.","gwas",
                             ".90.bi.remappedHa412HO.","Days_to_budding",".beagle.ps.gz"),
                      col_names=c("id","beta","sd","pvalue")) %>%
  separate(id,c("chr","pos"),"_") %>%
  mutate(pos=as.numeric(pos)) %>%
  semi_join(.,freqs) %>%
  filter(chr == "13") %>%
  mutate(type = "Full genome",logp=abs(log10(pvalue)))

all_gwas <- rbind(gwas_noSV,gwas_SV)

pdf("gwas/Ann_chr13dtb_peak.pdf",height=6,width=10)
all_gwas %>%
  filter(logp > 1) %>%
  ggplot(.,aes(x=pos/1000000,y=logp,color=type)) +
  geom_point(alpha=0.3) + 
  facet_wrap(~type,nrow=2) +
  ylab("Log p-value") + xlab("Mbp") +
  theme_minimal() +
  geom_segment(data=inv_associations %>% 
                 filter(id == "13_1"),
               aes(x=start/1000000,xend=end/1000000,
                   y=abs(log10(pvalue)),yend=abs(log10(pvalue))),
               color="#FFB31C",size=2,alpha=0.8) +
  scale_color_manual(values=c("grey","black")) +
  theme(legend.position = "none") +
  ggtitle("Chr13 - Days to bud")
dev.off()

sample_info <- read_tsv("sample_info_apr_2018.tsv") %>% rename(pop = population,sample=name)
pop_loc <- read_tsv("pop_loc_allnum.txt")
genotypes <- read_tsv("MDS_outliers/Ha412HO/annuus/Ha412HO_inv.v3.pcasites.Ha412HOChr13.pos1.genotypes.txt")
traits <- read_tsv("gwas/annuus/Days_to_budding.txt",col_names = c("sample","spacer","DTB"))
pca <- read_delim("PCA/Annuus.tranche90.snp.gwas.90.bi.remappedHa412HO.ldr0p2.eigenvectors.txt",delim=" ",
                  col_names = paste0("PCA",seq(1:614)))
cbind(traits %>% select(sample),pca ) %>% as_tibble() -> pca
af <- anova(lm(DTB ~ PCA1 + PCA2 +PCA3+triangle_genotype,
               genotypes %>%
                 inner_join(traits) %>%
                 inner_join(pca) %>%
                 inner_join(sample_info)
                 #filter(pop != "ANN_48", pop != "ANN_49")
                 ))

afss <- af$"Sum Sq"
#print(cbind(af,PctExp=afss/sum(afss)*100))

genotypes %>%
  inner_join(traits) %>%
  inner_join(pca) %>%
  inner_join(sample_info) %>%
  group_by(triangle_genotype) %>%
  summarize(mean_DTB = mean(DTB,na.rm=T),
            sd_DTB = sd(DTB,na.rm=T),
            DTB=mean(DTB,na.rm = T)) -> genotypes_summary

pdf("gwas/Ann_chr13dtb_effectsize.pdf",height=5,width=5)

genotypes %>%
  inner_join(traits) %>%
  inner_join(pca) %>%
  inner_join(sample_info) %>%
  ggplot(.,aes(x=as.factor(triangle_genotype),y=DTB,color=as.factor(triangle_genotype))) +
  geom_quasirandom() +
  geom_errorbar(data=genotypes_summary, 
                aes(x=as.factor(triangle_genotype),ymin = mean_DTB, ymax = mean_DTB),
                color="black",size=2,alpha=0.25) +
  scale_color_manual(values=c("#FFCE6E","#FFB31C","#D39117")) +
  ylab("Days to bud") + xlab("SV Genotype") +
  theme_bw() + 
  theme(legend.position = "none") 

genotypes %>%
  inner_join(traits) %>%
  inner_join(pca) %>%
  inner_join(sample_info) %>%
  ggplot(.,aes(x=as.factor(triangle_genotype),y=DTB,color=as.factor(triangle_genotype))) +
  geom_boxplot() +
  scale_color_manual(values=c("#FFCE6E","#FFB31C","#D39117")) +
  ylab("Days to bud") + xlab("SV Genotype") +
  theme_bw() + 
  theme(legend.position = "none")
  
  
dev.off()
  

usa <- map_data('state')
states <- map_data("state")
target_state <- map_data('state')
lat_range <- c(25, 50)
long_range <- c(-125,-93)
pie_size <- 0.4

pdf("gwas/Ann_chr13dtb_genomap.pdf",height=6,width=6)

ggplot(target_state, aes(long, lat)) +
  geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
  coord_quickmap() +
  geom_scatterpie(data=genotypes %>%
                    inner_join(sample_info) %>%
                    filter(species == "Ann") %>%
                    inner_join(pop_loc) %>%
                    group_by(pop, lat, long, triangle_genotype) %>% 
                    tally() %>%
                    spread(., triangle_genotype, n,fill=0),
                  aes(x=long, y=lat, r=pie_size), 
                  cols=c("0","1","2"), color=NA, alpha=.8) +
  scale_fill_manual(name="Genotype",values=c("#FFCE6E","#FFB31C","#D39117")) +theme_bw() +
  xlab("Longitude") + ylab("Latitude") +
  scale_x_continuous(limits = long_range, expand = c(0, 0)) +
  scale_y_continuous(limits = lat_range, expand = c(0, 0))  +
  ggtitle("H. annuus Chr13.pos1 SV") +
  theme_bw() 

dev.off()





