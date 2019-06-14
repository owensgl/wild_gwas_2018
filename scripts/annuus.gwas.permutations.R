library(tidyverse)
library(gridExtra)
library(grid)
#library(fitdistrplus)
directory <- "/media/owens/Copper/wild_gwas_2018/gwas"
prefix <- "Ha412HO_inv.jan09.pcasites.ANN."
#traits <- list.files(directory, paste(prefix,"*ps$",sep="")) %>%
#  str_remove(prefix) %>% str_remove(".ps$")
traits <- c( "TLN", "LIR", "Days_to_budding", "Stem_diamater_at_flowering",
             "Plant_height_at_flowering", "Primary_branches",
             "SLA_mm2_per_mg", "Leaf_total_N", "Leaf_total_C",
             "Disk_diameter", "Petal_length", "Petal_width",
             "Guides_3_petals", "Stem_colour", "Leaf_Area",
             "Leaf_Maximum_height", "Phyllaries_length", 
             "Phyllaries_width", "Seed_area", "Seed_HW_ratio")

inv_values <- tibble(snp_id=character(),beta=numeric(),se=numeric(),
                     pvalue=numeric(),type=character(),trait=character())
sig_snps <-  tibble(snp_id=character(),beta=numeric(),se=numeric(),
                    pvalue=numeric(),type=character(),trait=character())

perm_pvalues <- tibble(trait=character(),pvalue=numeric())
plot_list = list()
perm_plot_list = list()
for (i in 1:length(traits)){
  trait <- traits[i]
  
  
  #trait <- "Leaf_total_N"
  snps <- read_tsv(paste(directory,"/annuus/emmax/", "Annuus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.jan09noinv.",trait,".znorm.ps.gz",sep=""),
                   col_names = c("snp_id","beta","se","pvalue"))
  
  snps$type <- "snp"
  snps$trait <- trait
  sig_snps <- rbind(sig_snps,snps %>% filter(pvalue < 0.05))
  
  inv <- read_tsv(paste(directory,"/annuus/emmax/", "Ha412HO_inv.jan09.pcasites.annuus.",trait,".znorm.ps",sep=""),
                  col_names = c("snp_id","beta","se","pvalue"))
  
  
  inv$type <- "inv"
  inv$trait <- trait
  inv_values <- rbind(inv_values,inv)
  
  
  #permutation test 
  n_perm <- 5000
  min_inv_pvalue <- min(inv$pvalue)
  min_pvalue <- vector() 
  for (x in 1:n_perm){
    min_pvalue[x] <- sample_n(snps, nrow(inv)) %>% pull(pvalue) %>% min()
  }
  perm_pvalue <- sum(min_pvalue < min_inv_pvalue)/(n_perm+1)
  
  # perm_plot_list[[i]] <- ggplot(tibble(min_pvalue)) + 
  #   geom_density(aes(abs(log10(min_pvalue)))) +
  #   geom_vline(xintercept = abs(log10(min_inv_pvalue)),linetype="dashed") +
  #   theme_bw() +
  #   ylab("Density") + xlab("Minimum permuted log10(p-value)") +
  #   ggtitle(paste(trait,"\np-value =",round(perm_pvalue,3),sep=""))
  # 
  # 
  # plot_list[[i]] <- rbind(snps,inv) %>%
  #   ggplot(.,aes(abs(log10(pvalue)),fill=type)) + geom_density(alpha=0.5) +
  #   theme_bw() + scale_fill_brewer(palette = "Set1") +
  #   ggtitle(paste(trait,"\np-value =",round(perm_pvalue,3),sep=""))
  tmp <- tibble(trait=trait,pvalue=perm_pvalue)
  perm_pvalues <- rbind(perm_pvalues, tmp)
}





plot_inv_1 <- inv_values %>%
  mutate(significant = case_when(pvalue < 0.05 ~ 1,
                                 TRUE ~ 0)) %>%
  ggplot(.,aes(x=trait,y=snp_id,fill=as.factor(significant))) + geom_tile() +
  scale_fill_brewer(palette = "Set2",name="Pvalue < 0.05") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Inversion_ID")
plot_inv_2 <- inv_values %>%
  ggplot(.,aes(x=trait,y=snp_id,fill=abs(log10(pvalue)))) + geom_tile() +
  scale_fill_viridis(option="magma") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Inversion_ID")

write_tsv(inv_values, paste("gwas/",prefix,"invgwas.txt",sep=""))

pdf(paste("gwas/",prefix,"invgwas.permpvalues.pdf",sep=""),height=4,width=6)
perm_pvalues %>%
  ggplot(.,aes(pvalue)) + geom_histogram() +
  theme_bw() + ylab("Frequency") + xlab("p")
dev.off()

###Testing for pleotropy
sig_snps %>% 
  group_by(snp_id) %>%
  summarize(n_sig = n()) %>%
  mutate(type = "snp")-> sig_snp_counts

p.adjust(perm_pvalues$pvalue,method="BH")

inv_values %>%
  filter(pvalue < 0.05) %>%
  group_by(snp_id) %>%
  summarize(n_sig = n()) %>%
  mutate(type ="inv") -> sig_inv_counts

pleotropy_pvalue <- wilcox.test(sig_inv_counts$n_sig, sig_snp_counts$n_sig)

plot_pleotropy <- rbind(sig_inv_counts, sig_snp_counts) %>%
  ggplot(.,aes(x=type, y=n_sig,fill=type)) + geom_boxplot() +
  theme_bw() + 
  scale_fill_brewer(palette = "Set1") +
  ylab("Traits affected") +
  ggtitle(paste("Pleiotropy | p-value =", round(pleotropy_pvalue$p.value,3)))

pdf(paste("gwas/",prefix,"invgwas.permpleiotropy.pdf",sep=""),height=4,width=6)
sig_snp_counts %>%
  sample_n(10000) %>%
  rbind(., sig_inv_counts) %>%
  ggplot(.,aes(x=type, y=n_sig,fill=type)) + geom_boxplot() +
  theme_bw() + 
  scale_fill_brewer(palette = "Set1") +
  ylab("Traits affected") +
  ggtitle(paste("Pleiotropy | p-value =", round(pleotropy_pvalue$p.value,3)))
dev.off()


pdf(paste(prefix,"gwascomparison.pdf",sep=""),height=20,width=20)
lay <- rbind(c(1,2,3,4,5),
             c(6,7,8,9,10),
             c(11,12,13,14,15),
             c(16,17,18,19,20))
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],
             perm_plot_list[[1]],perm_plot_list[[2]],perm_plot_list[[3]],perm_plot_list[[4]],perm_plot_list[[5]],
             plot_list[[6]],plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],
             perm_plot_list[[6]],perm_plot_list[[7]],perm_plot_list[[8]],perm_plot_list[[9]],perm_plot_list[[10]],
             layout_matrix = lay)

grid.arrange(plot_list[[11]],plot_list[[12]],plot_list[[13]],plot_list[[14]],plot_list[[15]],
             perm_plot_list[[11]],perm_plot_list[[12]],perm_plot_list[[13]],perm_plot_list[[14]],perm_plot_list[[15]],
             plot_list[[16]],plot_list[[17]],plot_list[[18]],plot_list[[19]],plot_list[[20]],
             perm_plot_list[[16]],perm_plot_list[[17]],perm_plot_list[[18]],perm_plot_list[[19]],perm_plot_list[[20]],
             layout_matrix = lay)

grid.arrange(
  plot_inv_1, plot_inv_2,plot_pleotropy,
  layout_matrix = rbind(c(1,3),
                        c(2,3)))
dev.off()


