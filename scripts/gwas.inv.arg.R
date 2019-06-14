library(tidyverse)
library(lme4)
library(visreg)
library(car)
library(emmeans)
library(jtools)
library(qvalue)
prefix <- "Ha412HO_inv.v3.pcasites.ARG."
chosen_species <- "argophyllus"
inv_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt")
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T) %>% rename(sample=name)
traits <- c( "TLN", "LIR", "DTF","DTF_corrected","DTF_deduced",
             "late_early","Stem_diamater_at_flowering",
             "Plant_height_at_flowering", "Primary_branches",
             "SLA_mm2_per_mg", "Leaf_total_N", "Leaf_total_C",
             "Disk_diameter", "Petal_length", "Petal_width",
             "Stem_colour", "Leaf_Area",
             "Leaf_Maximum_height", "Phyllaries_length", 
             "Phyllaries_width", "Seed_area", "Seed_HW_ratio")

trait_types <- c("ft","ft","ft","ft","ft","ft","veg","veg","veg","veg","veg","veg","flower","flower","flower",
                 "veg","veg","veg","flower","flower","seed","seed")
trait_guide <- tibble(trait=as.character(traits),trait_types=as.character(trait_types))

pcavalues <- read_delim("PCA/Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO.v3noinv.ldfilter.ldr0p2.eigenvectors.txt",
                        col_names = paste("PC",1:614,sep=""),
                        delim=" ")

phenotypes <- read_tsv("gwas/arg/ARG_all_phenotypes_jan_2019.txt") %>%
  rename(sample = FID)
pcavalues$sample <- phenotypes$sample

sample_info <- read_tsv("sample_info_apr_2018.tsv") %>%
  rename(sample=name) %>%
  select(sample,population)

inversion_gwas <- tibble(estimate=character(),topCI=numeric(),botCI=numeric(),
                         tval=numeric(),p=numeric(),chr=character(),
                         direction=character(),mds=character(),trait=character())


for (i in 1:length(traits)){
  chosen_trait <- traits[i]
  if(!any(colnames(phenotypes) %in% chosen_trait)){
    next
  }
  mean_phenotype <- phenotypes %>% pull(chosen_trait) %>% mean(.,na.rm=T)
  sd_phenotype <- phenotypes %>% pull(chosen_trait) %>% sd(.,na.rm=T)
  scaled_phenotype <- phenotypes %>% pull(chosen_trait) %>% scale(.)
  phenotypes$chosen_phenotype <- scaled_phenotype
  
  
  inversion_genotypes <- phenotypes %>% select(sample)
  n_inv <- nrow(inv_list %>% filter(spe == !!chosen_species))
  counter <- 1
  inversion_names <- tibble(inv_n = character(),inv_list_n = numeric())
  for (n in 1:nrow(inv_list)){
    if (pull(inv_list[n,1]) != chosen_species){next}
    chosen_noncapital_species <- pull(inv_list[n,1])
    chosen_chr <- paste("Ha412HOChr",sprintf("%02d",pull(inv_list[n,4])),sep="")
    chosen_mds <- paste(pull(inv_list[n,6]),pull(inv_list[n,5]),sep="")
    inversion_n <- paste("inv_",counter,sep="")
    tmp <- tibble(inv_n =  paste("inv_",counter,sep=""), inv_list_n = n)
    inversion_names <- rbind(inversion_names, tmp)
    counter <- counter+ 1
    
    info <- read_tsv(paste("MDS_outliers/Ha412HO/",chosen_noncapital_species,"/Ha412HO_inv.v3.pcasites.",chosen_chr,".",chosen_mds,".genotypes.txt",sep="")) %>%
      select(sample, triangle_genotype) %>%
      mutate( triangle_genotype = triangle_genotype) %>%
      rename( !!inversion_n := triangle_genotype)
    
    inversion_genotypes <- inner_join(inversion_genotypes,info)
  }
  
  for (n in 1:n_inv){
    inversion_genotypes %>% 
      inner_join(.,pcavalues) %>%
      inner_join(.,phenotypes) %>%
      inner_join(.,sample_info) %>% 
      lm(as.formula(paste("chosen_phenotype ~ ",
                          paste(paste("PC",1:5,sep=""),collapse=" + "), 
                          " + ",
                          paste(paste("inv_",n,sep=""),collapse=" + "), 
                          sep="")),.)  -> lm_result
    popstr_result <- summ(lm_result,robust=F,confint = TRUE)
    
    popstr_inv <- tail(popstr_result$coeftable, 1)
    
    
    
    inversion_genotypes %>% 
      inner_join(.,pcavalues) %>%
      inner_join(.,phenotypes) %>%
      inner_join(.,sample_info) %>%
      lm(as.formula(paste("chosen_phenotype ~ ",
                          paste(paste("inv_",n,sep=""),collapse=" + "),
                          " -1", 
                          sep="")),.)  -> lm_result
    nostr_result <- summ(lm_result,robust=F,confint = TRUE)
    
    nostr_inv <- tail(nostr_result$coeftable, 1)
    
    result <- as.tibble(rbind(popstr_inv,nostr_inv))
    colnames(result) <- c("estimate","topCI","botCI","tval","p")
    result$chr <- inv_list %>% filter(spe == !!chosen_species) %>% pull(chr) %>% .[n]
    result$mds <- inv_list %>% filter(spe == !!chosen_species) %>% pull(mds) %>% .[n]
    result$direction <- inv_list %>% filter(spe == !!chosen_species) %>% pull(direction) %>% .[n]
    result$type <- c("popstr","nostr")
    result$trait <- chosen_trait
    inversion_gwas <- rbind(inversion_gwas, result)
  }
}
inversion_gwas_popstr %>%
  inner_join(.,trait_guide) %>%
  mutate(sig = case_when(q < 0.05 ~ "sig",
                         TRUE ~ "nonsig")) %>%
  inner_join(.,trait_guide) %>%
  write_tsv("gwas/arg/Ha412HO_inv.v3.allgwas.popstr.txt")

inversion_gwas %>%
  filter(type == "popstr") -> inversion_gwas_popstr
inversion_gwas_popstr$q <-  qvalue(inversion_gwas_popstr$p,lambda=0)$q
inversion_gwas_popstr %>%
  inner_join(.,trait_guide) %>%
  mutate(sig = case_when(q < 0.05 ~ "sig",
                         TRUE ~ "nonsig")) %>%
  ggplot(.) + geom_tile(aes(x=paste0(chr,"-",mds),y=trait,fill=sig))

inversion_gwas %>%
  filter(type == "nostr") -> inversion_gwas_nostr
inversion_gwas_nostr$q <-  qvalue(inversion_gwas_nostr$p, lambda=0)$q

inversion_gwas_nostr %>%
  inner_join(.,trait_guide) %>%
  mutate(sig = case_when(q < 0.05 ~ "sig",
                         TRUE ~ "nonsig")) %>%
  inner_join(.,trait_guide) %>%
  write_tsv("gwas/arg/Ha412HO_inv.v3.allgwas.nostr.txt")

inversion_gwas_nostr %>%
  inner_join(.,trait_guide) %>%
  mutate(trait = fct_reorder(trait, desc(trait_types))) %>%
  group_by(trait_types,trait) %>%
  count() %>%
  group_by(trait_types) %>%
  mutate(trait = fct_rev(trait)) %>%
  count() %>%
  ungroup() %>%
  mutate(trait_types = fct_rev(trait_types)) %>% map_df(rev)-> trait_groups

trait_groups %>%
  mutate(color=case_when(trait_types == "arc" ~"#609732",
                         trait_types == "color" ~"#A8383B",
                         trait_types == "flower" ~"#AFAD0A",
                         trait_types == "ft" ~"#FFD103",
                         trait_types == "leaf"~"#03C03C",
                         trait_types == "seed"~"#736C57")) %>%
  mutate(name=case_when(trait_types == "arc" ~"Architecture",
                        trait_types == "color" ~"Color",
                        trait_types == "flower" ~"Flower\nmorphology",
                        trait_types == "ft" ~"Flowering time",
                        trait_types == "leaf"~"Leaf traits",
                        trait_types == "seed"~"Seed traits")) -> trait_groups

trait_groups %>%
  ungroup() %>%
  mutate(cumsum=cumsum(n)) -> trait_groups
n_traits <- nrow(trait_groups)

n_inversions <- inversion_gwas_nostr %>%
  mutate(inversion = paste0("Chr",sprintf("%02d", chr),".",direction,mds)) %>%
  group_by(inversion) %>%
  count() %>% nrow()

trait_guide <- trait_guide[trait_guide$trait %in% inversion_gwas$trait,]
trait_guide$trait_n <- as.numeric(as.factor(trait_guide$trait_types))
inversion_gwas_nostr %>%
  inner_join(.,trait_guide) %>%
  mutate(sig = case_when(q < 0.05 ~ "sig",
                         TRUE ~ "nonsig")) %>%
  mutate(trait = fct_reorder(trait, trait_n)) %>%
  mutate(trait = fct_rev(trait)) %>%
  mutate(inversion = paste0("Chr",sprintf("%02d", chr),".",direction,mds)) %>%
  ggplot(.) + geom_tile(aes(x=inversion,y=trait),fill="white",color="black") +
  theme_bw() +
  geom_point(data = inversion_gwas_nostr %>%
               inner_join(.,trait_guide) %>%
               mutate(sig = case_when(q < 0.05 ~ "sig",
                                      TRUE ~ "nonsig")) %>%
               mutate(trait = fct_reorder(trait, trait_n)) %>%
               mutate(inversion = paste0("Chr",sprintf("%02d", chr),".",direction,mds)) %>%
               filter(sig == "sig"),
             aes(x=inversion,y=trait),color="blue",shape=1,size=4) +
  scale_fill_manual(values=c("black","white")) +
  geom_point(data=inversion_gwas_popstr %>%
               inner_join(.,trait_guide) %>%
               mutate(sig = case_when(q < 0.05 ~ "sig",
                                      TRUE ~ "nonsig")) %>%
               filter(sig == "sig"),
             aes(x=paste0("Chr",sprintf("%02d", chr),".",direction,mds),y=trait),shape='*',color="red",size=5) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        panel.grid.major.y = element_blank()) +
  coord_cartesian(xlim=c(0.5,n_inversions+4))+
  ylab("Trait") + xlab("Inversion") +
  ggtitle("Arg inversion GWAS") -> p

for (i in 1:n_traits){
  p <- p + annotate("rect", xmin=0.5,xmax=(n_inversions + 0.5),ymin=trait_groups$cumsum[i] - trait_groups$n[i] + 0.5,ymax=trait_groups$cumsum[i] + 0.5, 
                    alpha=0.2, fill=trait_groups$color[i],color="black")
  p <- p + annotate("text", x = (n_inversions + 0.75), y = trait_groups$cumsum[i] - (trait_groups$n[i]/2) +0.5, label = trait_groups$name[i],
                    hjust = 0,color=trait_groups$color[i])
}
p
pdf("gwas/Ha412HO_inv.v3.pcasites.ARG.gwas.v1.pdf",height=4,width=6)
p
dev.off()


