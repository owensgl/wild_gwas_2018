library(ggtree)
library(treeio)
library(tidyverse)
library(phytools)


inversion_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt")
sample_info <- read_tsv("sample_info_apr_2018.tsv") %>% rename(sample = name)

kate_set <- read_tsv("dune_sample_info.txt") %>%
  mutate(location = case_when(spe == "H. neglectus" ~ "MON",
                              spe == "H. petiolaris" ~ "GSD")) %>%
  rename(sample = name) %>%
  select(sample, location, ecotype)


dune_MON <- c("PET_49","PET_47")
sample_info %>% filter(population %in% dune_MON) %>%
  select(sample) %>% mutate(location = "MON",ecotype="dune") -> dune_mon_set

all_dune <- rbind(kate_set,dune_mon_set) %>% filter(ecotype == "dune")



#presentation_colors <- c("#FF0000","#A6A2A1","#000000","#000000","#EF36DC","#1294AB","#133EAA")
presentation_colors <- c("#0B2273","#455CAF","#803A00","#F7801E","#133EAA","#133EAA","#133EAA")
pdf("Testtree.pdf",height=25,width=10)
for (i in 1:nrow(inversion_list)){
  spe <- pull(inversion_list[i,1])
  if (spe != "petiolaris"){next}
  chosen_chr <- paste("Ha412HOChr",sprintf("%02d",pull(inversion_list[i,4])),sep="")
  chosen_mds <- paste(pull(inversion_list[i,6]),pull(inversion_list[i,5]),sep="")
  tree <- read.newick(paste("trees/",chosen_chr,".",chosen_mds,"/Petiolaris.tranche90.snp.petPF.90.bi.remappedHa412HO.",
                            chosen_chr,".",chosen_mds,".fasta.contree",sep=""))




  inversion_genotypes <- read_tsv(paste("MDS_outliers/Ha412HO/",
                                        spe,"/Ha412HO_inv.v3.pcasites.",
                                        chosen_chr,".",chosen_mds,".genotypes.txt",sep=""))
  sample_info %>% filter(sample %in% tree$tip.label) %>%
    select(sample, species,population)  %>% inner_join(inversion_genotypes) %>%
    mutate(ecotype = case_when(sample %in% all_dune$sample ~ "dune",
                               TRUE ~ "non-dune"))     %>%
    mutate(species_type = paste0(species,"-",triangle_genotype)) -> species_tree

  
  midroot.tree <- midpoint.root(tree)
  print(
   ggtree(midroot.tree)  %<+% species_tree +
      #geom_tiplab(size=2, align=TRUE, linesize=.5) +
      geom_tiplab(aes(color=species_type,label=population)) +
      #geom_text2(aes(subset = !isTip, label=label)) +
      geom_tippoint(aes(shape=as.factor(ecotype), fill=species), size=5,alpha=0.25) +
      theme(legend.position="right") +
      scale_color_manual(values = presentation_colors) +
     # xlim(0,0.03) +
      ggtitle(paste(spe," ",chosen_chr,".",chosen_mds,sep=""))  +theme_tree2()
      # gheatmap(p, species_heatmap[,2:3], offset=0, width=0.5, font.size=3, colnames_angle=-45, hjust=0) %>%
      #   scale_x_ggtree()
  )

}

dev.off()
