library(tidyverse)

env_associations <- read_tsv("/media/owens/Copper/wild_gwas/env_associations/for_Greg/var_out_annuus_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.txt",
                               col_names = c("chr_pos","chr","pos", "Latitude_e","Longitude_e","Elevation_e","MAT_e",
                                             "MWMT_e","MCMT_e","TD_e","MAP_e","MSP_e","AHM_e",
                                             "SHM_e","DD_0_e","DD5_e","DD_18_e","DD18_e","NFFD_e",
                                             "bFFP_e","eFFP_e","FFP_e","PAS_e","EMT_e","EXT_e","Eref_e",
                                             "CMD_e","MAR_e","RH_e")) %>%
  select(-chr,-pos) %>%
  separate(chr_pos, c("chr","pos"),"__") %>%
  mutate(pos = as.numeric(pos))

inversion_regions <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt")
chr_lengths <- read_tsv("Ha412HO.chrlengths.txt")

ann_inversions <- chr_lengths %>% 
  select(chr,end) %>%
  rename(chr_len = end) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(inversion_regions, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate( startcum=start+tot,endcum=end+tot) %>%
  filter(species == "Annuus") 

ann_inversions_length <- ann_inversions %>%
  mutate(length = end - start) %>%
  group_by(chr, mds) %>%
  summarize(length=round(sum(length)/1000000))

ann.env.cum <- chr_lengths %>% 
  select(chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(env_associations, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( poscum=pos+tot)

ann.env.cum %>%
  gather(.,variable,bayesfactor,Latitude_e:RH_e) -> ann.env.cum.long

axisdf = ann.env.cum %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2 )

unique(ann.env.cum.long$variable) -> variables

pdf("var_out_annuus_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.pdf",height=6,width=12)
for (i in 1:length(variables)){
  chosen_variable = variables[i]

plot <- ann.env.cum.long %>%
  filter(variable == chosen_variable) %>%
  filter(bayesfactor > 9) %>%
  ggplot(.) +
  
  # Show all points
  geom_point( aes(x=poscum, y=bayesfactor,color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  geom_hline(yintercept=10,linetype="dotted") +
  geom_rect(data=ann_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0,ymax=5)) +
  
  # custom X axis:
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  ggtitle(chosen_variable) +
  xlab("Chr")

print(plot)
}
dev.off()


ann.env.cum.long %>%
  filter(chr == "Ha412HOChr15") %>%
  filter(pos > 103772702, pos < 176634307) %>%
  filter(variable == "Latitude_e") %>%
  mutate(sig = case_when(bayesfactor > 10 ~ "sig",
                         TRUE ~ "nonsig")) %>%
  group_by(sig) %>%
  summarize(n=n()) %>%
  mutate(freq = n / sum(n))


pdf("var_out_annuus_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD.perchr.pdf",height=30,width=8)
for (i in 1:nrow(chr_lengths)){
  
  
  chosen_chr <- chr_lengths$chr[i]
  chosen_inversions <- ann_inversions %>%
    filter(chr == chosen_chr)
  
  
  if (nrow(chosen_inversions) > 0){
    plot <- ann.env.cum.long %>%
      filter(chr == chosen_chr) %>%
      filter(bayesfactor > 8) %>%
      ggplot(.) +
      
      # Show all points
      geom_point( aes(x=pos/1000000, y=bayesfactor),size=0.9, alpha=0.5) +
      geom_hline(yintercept=10,linetype="dotted") +
      geom_rect(data=chosen_inversions,aes(xmin=start/1000000,xmax=end/1000000,ymin=-0,ymax=5)) +
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
      xlab("Position") +
      facet_wrap(~variable,ncol=1,scales="free_y") +
      ggtitle(chosen_chr)
  }
  else{
    plot <- ann.env.cum.long %>%
      filter(chr == chosen_chr) %>%
      filter(bayesfactor > 8) %>%
      ggplot(.) +
      
      # Show all points
      geom_point( aes(x=pos/1000000, y=bayesfactor),size=0.9, alpha=0.5) +
      geom_hline(yintercept=10,linetype="dotted") +
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
      xlab("Position") +
      facet_wrap(~variable,ncol=1,scales="free_y") +
      ggtitle(chosen_chr)
    
  }
  print(plot)
}
dev.off()

