library(ggtree)
library(treeio)

inversion_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt")
sample_info <- read_tsv("sample_info_apr_2018.tsv") %>% rename(sample = name)

#presentation_colors <- c("#FF0000","#A6A2A1","#000000","#000000","#EF36DC","#1294AB","#133EAA")
presentation_colors <- c("#FF0000","#A6A2A1","#000000","#000000","#133EAA","#133EAA","#133EAA")
pdf("Ha412HO_inv.jan09.genetrees.petgroup.v0.pdf")
for (i in 1:nrow(inversion_list)){
  spe <- pull(inversion_list[i,1])
  chosen_chr <- paste("Ha412HOChr",sprintf("%02d",pull(inversion_list[i,4])),sep="")
  chosen_mds <- paste(pull(inversion_list[i,6]),pull(inversion_list[i,5]),sep="")
  tree <- read.newick(paste("/media/owens/Copper/wild_gwas_2018/inv_gen_dist/",
                            spe,"/",chosen_chr,".",chosen_mds,"/",chosen_chr,".",chosen_mds,".contree",sep=""))
                            
  inversion_genotypes <- read_tsv(paste("MDS_outliers/Ha412HO/",
                                        spe,"/Ha412HO_inv.jan09.pcasites.",
                                        chosen_chr,".",chosen_mds,".genotypes.txt",sep=""))
  sample_info %>% filter(sample %in% tree$tip.label) %>%
    select(sample, species)  %>% inner_join(inversion_genotypes) -> species_tree
      
  print(                
  ggtree(tree)  %<+% species_tree +
    geom_tiplab(aes(color=species)) +
    geom_text2(aes(subset = !isTip, label=label)) +
    geom_tippoint(aes(shape=as.factor(cluster_genotype), color=species), size=5,alpha=0.25) +
    theme(legend.position="right") +
    scale_color_manual(values = presentation_colors) +
    xlim(0,0.03) +
    ggtitle(paste(spe," ",chosen_chr,".",chosen_mds,sep=""))
  )
}

dev.off()



