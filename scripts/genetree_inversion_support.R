#Script for trying to plot support for different trees in inversions
library(tidyverse)
library(RColorBrewer)

inv_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt")

directory <- "/media/owens/Copper/wild_gwas_2018/inv_gen_dist"

treetests <- tibble(tree=character(),logL=numeric(),bprell=numeric(),bprellconf=character(),
                    pKH=numeric(),pKHconf=character(),pSH=numeric(),pSHconf=character(),
                    cELW=numeric(),cELWconf=character(),pAU=numeric(),pAUconf=character(),
                    n_gene=numeric(),species=character(),inv=character(),rel_position=numeric())
for (n in 1:18){
  species <- pull(inv_list[n,1])
  chosen_chr <- paste("Ha412HOChr",sprintf("%02d",pull(inv_list[n,4])),sep="")
  chosen_mds <- paste(pull(inv_list[n,6]),pull(inv_list[n,5]),sep="")
  inv <- paste(chosen_chr,".",chosen_mds,sep="")
  
  
  

  gene_n <- sort(as.numeric(gsub(".treetests.txt","",gsub("merged.","",list.files(paste(directory,"/",species,"/",inv,sep=""),"treetests.txt")))))
  treetest_columns <- c("tree","logL","deltaL","bprell","bprellconf",
                        "pKH","pKHconf","pSH","pSHconf",
                        "cELW","cELWconf","pAU","pAUconf")
  
  rel_position = 1
  for (i in 1:length(gene_n)){
    data <- read_tsv(paste(directory,"/",species,"/",inv,"/","merged.",i,".treetests.txt",sep=""),col_names = treetest_columns,skip=1)
    if (nrow(data) == 0){next}
    data$n_gene <- i
    data$species <- species
    data$inv <- inv
    data$rel_position <- rel_position
    rel_position <- rel_position + 1
    treetests <- rbind(treetests, data)
  }
}
arg_labels = c("((((Arg1,Arg2),Ann),Pet),Per)",
               "((((Arg2,Ann),Arg1),Pet),Per)",
               "((((Arg2,Ann),Pet),Arg1),Per)",
               "((((Arg1,Ann),Arg2),Pet),Per)",
               "((((Arg1,Ann),Pet),Arg2),Per)")

ann_labels = c("((((Ann1,Ann2),Arg),Pet),Per)",
               "((((Ann2,Arg),Ann1),Pet),Per)",
               "((((Ann2,Arg),Pet),Ann1),Per)",
               "((((Ann1,Arg),Ann2),Pet),Per)",
               "((((Ann1,Arg),Pet),Ann2),Per)")

pdf("Ha412HO_inv.jan09.genetreesupport.v0.pdf",height=10,width=20)
for (chosen_species in c("annuus","argophyllus")){
  
  
  if(chosen_species == "annuus"){
    chosen_labels = ann_labels
  }else if(chosen_species == "argophyllus"){
    chosen_labels = arg_labels
  }
  print(
    treetests %>%
      filter(species == chosen_species) %>%
      ggplot(.,aes(x=rel_position,y=pAU,fill=as.character(tree))) + geom_bar(stat="identity") +
      facet_grid(tree~inv,scales="free_x") + theme_bw() + 
      scale_fill_brewer(palette = "Set1",
                        labels=chosen_labels,
                        name="Tree") +
      xlab("Gene_position") +
      ylab("AU support") +
      ggtitle(paste(chosen_species))
  )
  set1 <-brewer.pal(5,"Set1")
  print(
    treetests %>%
      filter(species == chosen_species) %>%
      filter(pAU > 0.80) %>%
      ggplot(.,aes(x=inv,fill=as.character(tree))) + geom_bar(stat="count",position="fill",width=1) +
      theme_bw() + 
      scale_fill_manual(values=set1,
                        labels=chosen_labels,
                        limits=c("1", "2", "3","4","5"),
                        name="Tree",drop=FALSE) +
      ylab("Proportion supporting phylogeny") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      ggtitle(paste(chosen_species))
  )
}
dev.off()
  


     