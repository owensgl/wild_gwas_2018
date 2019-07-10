library(tidyverse)
library(ggrepel)

species_list <- c("annuus","argophyllus","petfal","petpet")
chr_lengths <- read_tsv("Ha412HO.chrlengths.txt")
for (x in 1:4){
  chosen_species <- species_list[x]
  column_names <- colnames(read_tsv(paste0("/media/owens/Copper/wild_gwas/env_associations/soil/column.names.",chosen_species,".varout.txt")))
  
  env_associations <- read_tsv(paste0("/media/owens/Copper/wild_gwas/env_associations/soil/",chosen_species,"_HA412_all_covariates_baypass.BF.txt"),
                               col_names = column_names) %>%
    select(-chrom) %>%
    separate(snp_id, c("chr","pos"),"__") %>%
    mutate(pos = as.numeric(pos))
  
  env.cum <- chr_lengths %>% 
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
  
  env.cum %>%
    gather(.,variable,bayesfactor,OM:SOL_SALTS) -> env.cum.long
  
  axisdf = env.cum %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2 )
  
  # ##Load in candidate genes for labelling. For the Shaghayegh based filtering.
  # genes <- read_tsv("gea/GEA_Genes_With_Description.txt")
  # location <- read_tsv("gea/Han412-HO-only-Genes_with_IDs_shaghayegh.gff",col_names = c("chr","type","type2","start","end","spacer1",
  #                                                                                        "spacer2","spacer3","info","geneID")) %>%
  #   select(chr,info,start,end) %>%
  #   separate(info,c("gene_name","name","extra"),";") %>%
  #   select(-name,-extra) %>%
  #   mutate(gene_name = gsub("ID=gene:","",gene_name)) %>%
  #   mutate(mid = ((end-start)/2)+start) 
  # 
  # genes %>% select(gene_name,P.Val,test_name) %>%
  #   inner_join(.,location) %>% unique()-> genes
  # 
  # ann.genes.cum <- chr_lengths %>% 
  #   select(chr,end) %>%
  #   rename(chr_len = end) %>%
  #   
  #   # Calculate cumulative position of each chromosome
  #   mutate(tot=cumsum(chr_len)-chr_len) %>%
  #   select(-chr_len) %>%
  #   
  #   # Add this info to the initial dataset
  #   left_join(genes, ., by=c("chr"="chr")) %>%
  #   
  #   # Add a cumulative position of each SNP
  #   arrange(chr, mid) %>%
  #   mutate( poscum=mid+tot)
  
  genes <- read_tsv(paste0("/media/owens/Copper/wild_gwas/env_associations/soil/",chosen_species,
                           "_HA412_all_covariates_baypass.BF.genehits.txt")) %>%
    mutate(mid = ((end-start)/2)+start) 
  genes.cum <- chr_lengths %>%
    select(chr,end) %>%
    rename(chr_len = end) %>%
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(genes, ., by=c("chr"="chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(chr, mid) %>%
    mutate( poscum=mid+tot)
  
  unique(env.cum.long$variable) -> variables
  #Figure out top hit for each gene
  genes.cum %>%
    arrange(desc(second_hit)) %>%
    group_by(variable) %>%
    slice(1:10) -> top_hits
  
  pdf(paste0("gea/",chosen_species,"_HA412_all_covariates_baypass.BF.pdf"),height=6,width=12)
  for (i in 1:length(variables)){
    chosen_variable = variables[i]
    
    plot <- env.cum.long %>%
      filter(variable == chosen_variable) %>%
      filter(bayesfactor > 9) %>%
      ggplot(.) +
      
      # Show all points
      geom_point( aes(x=poscum, y=bayesfactor,color=as.factor(chr)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
      geom_hline(yintercept=10,linetype="dotted") +
      #geom_rect(data=ann_inversions,aes(xmin=startcum,xmax=endcum,ymin=-0,ymax=5)) +
      
      # custom X axis:
      scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
      #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
      
      # Custom the theme:
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
      ggtitle(chosen_variable) +
      xlab("Chr")  +
      geom_text_repel(
        data = top_hits %>%filter(tolower(variable) == tolower(chosen_variable)),
        aes(x=poscum,y=as.numeric(top_hit),label=gene),
        segment.size  = 0.2,
        segment.color = "grey50",
        direction     = "y"
      ) +
      geom_point(data=top_hits %>%filter(tolower(variable) == tolower(chosen_variable)),
                 aes(x=poscum, y=as.numeric(top_hit)), color="black", alpha=1, size=2)
    
    print(plot)
  }
  dev.off()
  
}



