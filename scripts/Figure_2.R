library(tidyverse)
library(statebins)
library(grid)
library(gridExtra)

species <- tibble(chosen_species=c("annuus","argophyllus","petiolaris","petiolaris"),
                  chosen_species_long=c("annuus","argophyllus","petiolaris.fallax","petiolaris.petiolaris"),
                  chosen_abbreviation=c("ann","arg","petfal","petpet"))

myplots <- list()
for (x in 1:nrow(species)){
  chosen_species <- species$chosen_species[x]
  chosen_species_long <- species$chosen_species_long[x]
  chosen_abbreviation <- species$chosen_abbreviation[x]



inv_locations <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt")

chr_lengths <- read_tsv("Ha412HO.chrlengths.txt")
chr_lengths %>% 
  select(chr,end) %>%
  rename(chr_len = end) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(inv_locations, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate( startcum=start+tot,endcum=end+tot) %>%
  filter(species != "Niveus") %>%
  separate(chr, c("ref","chr_n"), "Chr") %>%
  mutate(chr_n = as.numeric(chr_n)) -> inv_locations_formatted


read_tsv(paste0("gwas/",chosen_abbreviation,"/Ha412HO_inv.v3.allgwas.popstr.txt")) %>%
  select(trait_types) %>% unique() %>%
  mutate(color=case_when(trait_types == "veg" ~"#609732",
                         trait_types == "flower" ~"#AFAD0A",
                         trait_types == "ft" ~"#FFD103",
                         trait_types == "seed"~"#BC3908")) %>%
  mutate(name=case_when(trait_types == "flower" ~"Flower",
                        trait_types == "ft" ~"Phenology",
                        trait_types == "veg"~"Vegetatitive",
                        trait_types == "seed"~"Seed"),
         y = case_when(trait_types == "flower" ~2.5,
                       trait_types == "ft" ~3.5,
                       trait_types == "veg"~2.0,
                       trait_types == "seed"~3.0)) -> trait_looks



read_tsv("gea/climate_variables.txt") %>%
  select(type) %>% unique() %>%
  rename(variable_type= type) %>%
  mutate(color=case_when(variable_type == "Temperature" ~ "#60D5EB",
                         variable_type == "Precipitation" ~ "#7486CE",
                         variable_type == "Combined" ~ "#7588AE",
                         variable_type == "Location" ~ "#C0FDFB"),
         y = case_when(variable_type == "Temperature" ~ 0,
                       variable_type == "Precipitation" ~ -1,
                       variable_type == "Combined" ~ -0.5,
                       variable_type == "Location" ~ -1.5)) %>%
  filter(variable_type != "Location") -> gea_looks

read_tsv("gea/climate_variables.txt") %>%
  rename(variable_type= type,climate_variable = variable) %>%
  mutate(color=case_when(variable_type == "Temperature" ~ "#60D5EB",
                         variable_type == "Precipitation" ~ "#7486CE",
                         variable_type == "Combined" ~ "#7588AE",
                         variable_type == "Location" ~ "#C0FDFB"),
         y = case_when(variable_type == "Temperature" ~ 0,
                       variable_type == "Precipitation" ~ -1,
                       variable_type == "Combined" ~ -0.5,
                       variable_type == "Location" ~ -1.5)) -> gea_traits

max_length <- chr_lengths %>%
  mutate(cumstart=cumsum(end)-end,
         cumend=cumsum(end)) %>%
  summarize(end=max(cumend)) %>% pull(end)
inv_windows <- inv_locations_formatted %>%
  filter(spe == chosen_species) %>%
  select(chr_n, mds) %>%
  unique() %>%
  mutate(n=0:(n()-1)) %>% 
  mutate(window_start = (max_length/n())*n,
         window_end = (max_length/n())*(n+1)) 
  
#Load in GWAS data
read_tsv(paste0("gwas/",chosen_abbreviation,"/Ha412HO_inv.v3.allgwas.popstr.txt")) %>% 
  mutate(mds2 = paste0(direction,mds)) %>%
  select(-direction,-mds) %>%
  rename(chr_n = chr) %>%
  rename(mds=mds2) %>%
  inner_join(inv_locations_formatted %>%
               filter(spe == chosen_species) %>%
               group_by(chr_n,mds) %>%
               summarize(chr_start=min(startcum),chr_end=max(endcum))) %>%
  inner_join(inv_windows) %>%
  group_by(chr_n,mds,trait_types,window_start,window_end,chr_start,chr_end) %>%
  summarize(sig_count = sum(sig == "sig"),
            max_sig = max(abs(log10(q)))) %>% 
  mutate(row = case_when(trait_types == "veg" ~ 2,
                         trait_types == "flower" ~2.5,
                         trait_types == "seed" ~ 3,
                         trait_types == "ft" ~ 3.5),
         any_sig = case_when(sig_count > 0 ~ "significant",
                             TRUE ~ "nonsignificant"))-> gwas_structure

gwas_structure %>%
  ungroup() %>%
  filter(trait_types == "veg") %>%
  select(chr_n,mds, window_start,window_end,chr_start,chr_end) %>%
  mutate(pt1=paste(chr_start,1.5,sep=" "),
         pt2=paste(window_start,2,sep=" "),
         pt3=paste(window_end,2,sep=" "),
         pt4=paste(chr_end,1.5,sep=" ")) %>%
  select(chr_n,mds, pt1,pt2,pt3,pt4) %>%
  gather(row,value,pt1:pt4) %>%
  separate(value, c("x","y")," ") %>%
  mutate(x = as.numeric(x),y=as.numeric(y),
         id=paste0(chr_n,mds))-> gwas_polygons


#Load in GEA data


gea <- read_delim(paste0("gea/H_",chosen_species_long,"_HA412_remapped/baypass/tables/varout_",chosen_species_long,"_HA412_remapped_Baypass_noinv_matrix_25_vars.tab.table"),
                  delim=" ") %>% 
  separate(inversion_ID, c("chr","mds")) %>%
  mutate(chr_n = gsub("Ha412HOChr","",chr)) %>% 
  select(-chr) %>%
  select(-latitude_e,-longitude_e,-elevation_e) %>%
  gather(climate_variable, bayesfactor, MAT_e:RH_e) %>%
  inner_join(gea_traits) %>%
  group_by(chr_n,mds,variable_type,color,y) %>%
  mutate(bayesfactor = as.numeric(bayesfactor)) %>%
  summarize(sig_count = sum(abs(bayesfactor) > 10),
            max_sig = max(abs(bayesfactor))) %>%
  mutate(any_sig = case_when(sig_count > 0 ~ "significant",
                             TRUE ~ "nonsignificant")) %>%
  ungroup() %>%
  mutate(chr_n = as.numeric(chr_n)) %>%
  inner_join(inv_locations_formatted %>%
               filter(spe == chosen_species) %>%
               group_by(chr_n,mds) %>%
               summarize(chr_start=min(startcum),chr_end=max(endcum))) %>%
  inner_join(inv_windows) %>% ungroup()
  

gea %>%
  ungroup() %>%
  select(chr_n,mds, window_start,window_end,chr_start,chr_end) %>%
  mutate(pt1=paste(chr_start,1,sep=" "),
         pt2=paste(window_start,0.5,sep=" "),
         pt3=paste(window_end,0.5,sep=" "),
         pt4=paste(chr_end,1,sep=" ")) %>%
  select(chr_n,mds, pt1,pt2,pt3,pt4) %>%
  gather(row,value,pt1:pt4) %>%
  separate(value, c("x","y")," ") %>%
  mutate(x = as.numeric(x),y=as.numeric(y),
         id=paste0(chr_n,mds))-> gea_polygons





p <- chr_lengths %>%
  mutate(cumstart=cumsum(end)-end,
         cumend=cumsum(end)) %>%
  ggplot(.) +
  geom_rect(data=inv_locations_formatted %>%
              filter(spe == chosen_species),aes(xmin=startcum/1000000,xmax=endcum/1000000,ymin=1,ymax=1.5),fill="#F6511D",color=NA) +
  statebins:::geom_rrect(aes(xmin=cumstart/1000000,xmax=cumend/1000000,ymin=1,ymax=1.5),fill="NA",color="black",radius=grid::unit(3, "pt")) +
  scale_fill_manual(values=c("white","black")) +
  geom_polygon(data=gwas_polygons, aes(x=x/1000000,y=y,group=id),fill="grey",color="black") +
  geom_polygon(data=gea_polygons, aes(x=x/1000000,y=y,group=id),fill="grey",color="black") +
  theme_bw() +
  coord_cartesian(xlim=c(-1500,3300),ylim=c(-2,7)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  geom_text(data=inv_windows,aes(x=((window_start+(window_end-window_start)/2)/1000000),y=4.1,label=paste0("Chr",chr_n,".",mds)),
            angle = 90,hjust=0) 



#Add in the chromosome labels


#Add in the GWAS labels and colors  
for (i in 1:4){
  p <- p + annotate("text", x = -100, y = trait_looks$y[i]+0.25, label = trait_looks$name[i],
                  hjust = 1,color=trait_looks$color[i])
  p <- p + annotate("rect",xmin=0,xmax=max_length/1000000,ymin=trait_looks$y[i],ymax=trait_looks$y[i]+0.5,
                  fill=trait_looks$color[i],color="black")
}
#Add in the GEA label and colors
for (i in 1:nrow(gea_looks)){
  p <- p + annotate("text", x = -100, y = gea_looks$y[i]+0.25, label = gea_looks$variable_type[i],
                    hjust = 1,color=gea_looks$color[i])
  p <- p + annotate("rect",xmin=0,xmax=max_length/1000000,ymin=gea_looks$y[i],ymax=gea_looks$y[i]+0.5,
                    fill=gea_looks$color[i],color="black")
}
p <- p +  geom_rect(data=gwas_structure,aes(xmin=window_start/1000000,xmax=window_end/1000000,
                                                 ymin=row,ymax=row+0.5,fill=any_sig),color="black",alpha=0.5)
p <- p +   geom_rect(data=gea,
                     aes(xmin=window_start/1000000,xmax=window_end/1000000,
                         ymin=y,ymax=y+0.5,fill=any_sig),color="black",alpha=0.5) 
p <- p +   geom_text(data=gwas_structure %>% filter(sig_count > 1) %>% 
                       mutate(string = gsub("(.{3})", "\\1\n", paste(rep('*',floor(max_sig/5)),collapse = ""))),
                     aes(((x=window_start+(window_end-window_start)/2)/1000000),y=row+0.25,
                         label=string),lineheight = 0.5,color="white",size=5)
p <- p +   geom_text(data=gea %>% filter(sig_count > 1) %>% 
                       group_by(chr_n, mds, variable_type, window_start, window_end, chr_start) %>%
                       mutate(string = gsub("(.{3})", "\\1\n", paste(rep('*',floor((max_sig-10)/5)),collapse = ""))),
                     aes(((x=window_start+(window_end-window_start)/2)/1000000),y=y+0.25,
                         label=string),lineheight = 0.5,color="white",size=5)
p <- p + annotate("text", x = sum(chr_lengths$end)/2000000, y = -1.5, label = paste("H.",chosen_species_long),
                  hjust = 0.5,color="black")
myplots[[x]] <-p

}
pdf("inversion_effects_test.v1.pdf",height=10,width=10)
grid.arrange(myplots[[1]],myplots[[3]],
             myplots[[2]],myplots[[4]],
             nrow = 2)
dev.off()
