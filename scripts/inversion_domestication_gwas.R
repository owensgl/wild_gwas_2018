#This is for checking direction of inversion effect.
library(tidyverse)
library(lme4)
library(visreg)
library(car)
library(emmeans)
library(jtools)
prefix <- "Ha412HO_inv.jan09.pcasites.ANN."
chosen_species <- "annuus"
inv_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt")
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T) %>% rename(sample=name)
traits <- c( "TLN", "LIR", "Days_to_budding", "Stem_diamater_at_flowering",
             "Plant_height_at_flowering", "Primary_branches",
             "SLA_mm2_per_mg", "Leaf_total_N", "Leaf_total_C",
             "Disk_diameter", "Petal_length", "Petal_width",
             "Guides_3_petals", "Stem_colour", "Leaf_Area",
             "Leaf_Maximum_height", "Phyllaries_length", 
             "Phyllaries_width", "Seed_area", "Seed_HW_ratio")

pcavalues <- read_delim("Annuus.tranche99.snp.gwas.90.bi.remappedHa412HO.jan09noinv.ldr0p2.eigenvectors.txt",
                        col_names = paste("PC",1:614,sep=""),
                        delim=" ")


phenotypes <- read_tsv("/media/owens/Copper/wild_gwas_2018/gwas/annuus/ANN_all_phenotypes_jan_2019.txt") %>%
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
  mean_phenotype <- phenotypes %>% pull(chosen_trait) %>% mean(.,na.rm=T)
  sd_phenotype <- phenotypes %>% pull(chosen_trait) %>% sd(.,na.rm=T)
  scaled_phenotype <- phenotypes %>% pull(chosen_trait) %>% scale(.)
  phenotypes$chosen_phenotype <- scaled_phenotype
  
  
  inversion_genotypes <- phenotypes %>% select(sample)
  n_inv <- nrow(inv_list %>% filter(spe == !!chosen_species))
  for (n in 1:n_inv){
    chosen_noncapital_species <- pull(inv_list[n,1])
    chosen_chr <- paste("Ha412HOChr",sprintf("%02d",pull(inv_list[n,4])),sep="")
    chosen_mds <- paste(pull(inv_list[n,6]),pull(inv_list[n,5]),sep="")
    inversion_n <- paste("inv_",n,sep="")
    
    info <- read_tsv(paste("MDS_outliers/Ha412HO/",chosen_noncapital_species,"/Ha412HO_inv.jan09.pcasites.",chosen_chr,".",chosen_mds,".genotypes.txt",sep="")) %>%
      select(sample, triangle_genotype) %>%
      mutate( triangle_genotype = triangle_genotype) %>%
      rename( !!inversion_n := triangle_genotype)
    
    inversion_genotypes <- inner_join(inversion_genotypes,info)
  }

  for (n in 1:n_inv){
    inversion_genotypes %>% inner_join(.,pcavalues) %>%
      inner_join(.,phenotypes) %>%
      inner_join(.,sample_info) %>%
      lm(as.formula(paste("chosen_phenotype ~ ",
                          paste(paste("PC",1:20,sep=""),collapse=" + "), 
                          " + ",
                          paste(paste("inv_",n,sep=""),collapse=" + "), 
                          sep="")),.)  -> lm_result
    popstr_result <- summ(lm_result,robust=T,confint = TRUE)
    
    popstr_inv <- tail(popstr_result$coeftable, 1)
    
    inversion_genotypes %>% inner_join(.,pcavalues) %>%
      inner_join(.,phenotypes) %>%
      inner_join(.,sample_info) %>%
      lm(as.formula(paste("chosen_phenotype ~ ",
                          paste(paste("inv_",n,sep=""),collapse=" + "), 
                          sep="")),.)  -> lm_result
    nostr_result <- summ(lm_result,robust=T,confint = TRUE)
    
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


domestic <- read_tsv("MDS_outliers/Ha412HO/annuus/Ha412HO_inv.v3.domesticalleles.txt")
pdf(paste0(prefix,"domestic_effect_size.pdf"),height=6,width=8)
for (i in 1:length(traits)){
  chosen_trait <- traits[i]
  print(
    inversion_gwas %>%
      filter(trait == chosen_trait) %>%
      inner_join(., domestic) %>%
      filter(chr != "1", chr != "5") %>%
      mutate(estimate = case_when(domestic_allele == 0 ~ -estimate,
                                  TRUE ~ estimate),
             topCI = case_when(domestic_allele == 0 ~ -topCI,
                               TRUE ~ topCI),
             botCI = case_when(domestic_allele == 0 ~ -botCI,
                               TRUE ~ botCI)) %>%
      mutate(type = case_when(type == "nostr" ~ "No controls",
                              type == "popstr" ~ "Population structure controls")) %>%
      mutate(inversion_name = paste0("Ha412HOChr",chr,"-",direction,mds)) %>%
      ggplot(.,aes(x=inversion_name,y=estimate)) + geom_point() +
      geom_errorbar(aes(ymin=botCI, ymax=topCI), colour="black", width=.1) +
      facet_wrap(~type,nrow=2) +
      geom_hline(yintercept=0,linetype="dotted") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      ggtitle(chosen_trait) +
      ylab("LM Effect size") + xlab("Inversion")
  )
}
dev.off()





