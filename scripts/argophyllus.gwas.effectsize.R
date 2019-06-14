library(tidyverse)
library(lme4)
library(visreg)
library(car)
library(emmeans)
prefix <- "Ha412HO_inv.jan09.pcasites.argophyllus."
chosen_species <- "argophyllus"
chosen_abbreviation <- "ARG"
chosen_noncapital_species <- "argophyllus"
chosen_capital_species <- "Argophyllus"
inv_list <- read_tsv("MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt") %>%
  filter(spe == chosen_noncapital_species)
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv",col_names = T) %>% rename(sample=name)
traits <- c( "TLN", "LIR", "DTF", "Stem_diamater_at_flowering",
             "Plant_height_at_flowering", "Primary_branches",
             "SLA_mm2_per_mg", "Leaf_total_N", "Leaf_total_C",
             "Disk_diameter", "Petal_length", "Petal_width",
             "Guides_3_petals", "Stem_colour", "Leaf_Area",
             "Leaf_Maximum_height", "Phyllaries_length", 
             "Phyllaries_width", "Seed_area", "Seed_HW_ratio")


pcavalues <- read_delim("Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO.jan09noinv.ldfilter.ldr0p2.eigenvectors.txt",
                        col_names = paste("PC",1:614,sep=""),
                        delim=" ")
rerun_snp_grab <- TRUE

phenotypes <- read_tsv(paste("/media/owens/Copper/wild_gwas_2018/gwas/",chosen_species,"/",chosen_abbreviation,"_all_phenotypes_jan_2019.txt",sep="")) %>%
  rename(sample = FID)
pcavalues$sample <- phenotypes$sample

sample_info <- read_tsv("sample_info_apr_2018.tsv") %>%
  rename(sample=name) %>%
  select(sample,population)

all_top_snps <- tibble(snp_id=character(),chr=numeric(),
                       pos=numeric(),beta=numeric(),se=numeric(),
                       pvalue=numeric(),trait=character())
PVE <- tibble(type = character(),
              value = numeric(),
              trait = character(),
              controlled=character())

for (i in 1:length(traits)){
  chosen_trait <- traits[i]
  if (!any(names(phenotypes) == chosen_trait)){
    next
  }
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
      rename( !!inversion_n := triangle_genotype)
    inversion_genotypes <- inner_join(inversion_genotypes,info)
  }
  inversion_genotypes %>% inner_join(.,pcavalues) %>%
    inner_join(.,phenotypes) %>%
    inner_join(.,sample_info) %>%
    lm(as.formula(paste("chosen_phenotype ~ ",
                        paste(paste("PC",1:20,sep=""),collapse=" + "), 
                        sep="")),.)  -> lm_result_popstr
  anova(lm_result_popstr) -> anova_result_popstr
  #inv_pvalue  <- formatC(anova_result$`Pr(>F)`[51],format="e",digits=3)
  total_popstr_PVE <- round(sum(anova_result_popstr$`Sum Sq`[1:(length(anova_result_popstr$`Sum Sq`) - 1)])/
                              sum(anova_result_popstr$`Sum Sq`) * 100,2)
  
  
  inversion_genotypes %>% inner_join(.,pcavalues) %>%
    inner_join(.,phenotypes) %>%
    inner_join(.,sample_info) %>%
    lm(as.formula(paste("chosen_phenotype ~ ",
                        paste(paste("PC",1:20,sep=""),collapse=" + "), 
                        " + ",
                        paste(paste("inv_",1:n_inv,sep=""),collapse=" + "), 
                        sep="")),.)  -> lm_result
  
  anova(lm_result) -> anova_result
  #inv_pvalue  <- formatC(anova_result$`Pr(>F)`[51],format="e",digits=3)
  total_inv_PVE <- round(sum(anova_result$`Sum Sq`[1:(length(anova_result$`Sum Sq`) - 1)])/
                           sum(anova_result$`Sum Sq`) * 100,2) - total_popstr_PVE
  
  #Calculate a non-controled LMM
  inversion_genotypes %>% inner_join(.,pcavalues) %>%
    inner_join(.,phenotypes) %>%
    inner_join(.,sample_info) %>%
    lm(as.formula(paste("chosen_phenotype ~ ",
                        paste(paste("inv_",1:n_inv,sep=""),collapse=" + "), 
                        sep="")),.)  -> lm_result
  
  anova(lm_result) -> anova_result
  #inv_pvalue  <- formatC(anova_result$`Pr(>F)`[51],format="e",digits=3)
  total_inv_PVE_uncontrolled <- round(sum(anova_result$`Sum Sq`[1:(length(anova_result$`Sum Sq`) - 1)])/
                                        sum(anova_result$`Sum Sq`) * 100,2)
  
  if (rerun_snp_grab){
    
    ###Read in gwas and pick topic hits in different areas
    directory <- "/media/owens/Copper/wild_gwas_2018/gwas"
    if (!file_test("-f", paste(directory,"/",chosen_species,"/", chosen_capital_species,".tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.jan09noinv.ldfilter.",chosen_trait,".znorm.ps.gz",sep=""))){
      next
    }
    snps <- read_tsv(paste(directory,"/",chosen_species,"/", chosen_capital_species,".tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.jan09noinv.ldfilter.",chosen_trait,".znorm.ps.gz",sep=""),
                     col_names = c("snp_id","beta","se","pvalue"))
    
    snps %>%
      arrange(pvalue) %>% head(n=1000) %>%
      separate(snp_id,c("chr","pos"),  "_",remove=F)-> top_snps
    
    
    min_distance <- 100000 #Minimum distance between snps looked at.
    removed_positions <- vector()
    for (x in 2:1000){
      #Check if any of the next snps are within the minimum distance, and remove them
      target_pos <- top_snps$pos[x]
      target_chr <- top_snps$chr[x]
      matches <- top_snps %>% head(n = (x-1)) %>%
        filter(chr == target_chr) %>%
        filter(abs(as.numeric(pos) - as.numeric(target_pos)) < min_distance) %>% 
        nrow() 
      if(matches > 0){
        removed_positions <- c(removed_positions,x)
      }
    }
    top_snps[-removed_positions,] %>% head(n=100) -> top_snps
    top_snps$trait <- chosen_trait
    #Save values for this loop
    all_top_snps <- rbind(top_snps, all_top_snps)
  }
  PVE_tmp <- tibble(type = c("population_structure", "inv","inv"),
                    value = c(total_popstr_PVE, total_inv_PVE,total_inv_PVE_uncontrolled),
                    trait = c(chosen_trait, chosen_trait,chosen_trait),
                    controlled = c("controlled","controlled","uncontrolled"))
  
  PVE <- rbind(PVE, PVE_tmp)
}
#Pull all the top snps (for all traits), out of the tped.

top_snp_ids <- all_top_snps %>% pull(snp_id)

read_tsv(paste(directory,"/",chosen_species,"/",chosen_capital_species, ".tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.jan09noinv.ldfilter.tped",sep=""),
         col_names=c("chr","snp_id","fill","pos",rep(phenotypes$sample, each=2))) %>%
  filter(snp_id %in% top_snp_ids) -> snp_data


snp_data %>%
  select(-chr,-fill,-pos) %>%
  gather(.,sample,genotype, -snp_id) %>%
  mutate(sample = str_remove(sample, "_1")) %>%
  group_by(snp_id,sample) %>%
  summarize(genotype = sum(genotype)) %>%
  ungroup() %>%
  mutate(genotype = case_when(genotype == "0" ~ "NA",
                              TRUE ~ as.character(genotype - 2))) -> snp_data_reformatted
#Now take all the snps that are saved and use them to calculate trait specific PVE
for (i in 1:length(traits)){
  chosen_trait <- traits[i]
  if (!any(names(phenotypes) == chosen_trait)){
    next
  }
  mean_phenotype <- phenotypes %>% pull(chosen_trait) %>% mean(.,na.rm=T)
  sd_phenotype <- phenotypes %>% pull(chosen_trait) %>% sd(.,na.rm=T)
  scaled_phenotype <- phenotypes %>% pull(chosen_trait) %>% scale(.)
  phenotypes$chosen_phenotype <- scaled_phenotype
  
  chosen_snps <- all_top_snps %>% filter(trait == chosen_trait) %>% pull(snp_id)
  snp_data_reformatted %>% 
    filter(snp_id %in% chosen_snps) %>%
    mutate(snp_id = paste("snp_",as.numeric(as.factor(snp_id)),sep="")) %>%
    spread(snp_id, genotype) %>%
    inner_join(.,pcavalues) %>%
    inner_join(.,phenotypes %>% select(sample, chosen_phenotype)) %>%
    inner_join(.,sample_info) %>%
    #select(sample,population, genotype,PC1:PC614, chosen_phenotype) %>% 
    lm(as.formula(paste("chosen_phenotype ~ ",
                        paste(paste("PC",1:20,sep=""),collapse=" + "), 
                        " + ",
                        paste(paste("snp_",1:100,sep=""),collapse=" + "), 
                        sep="")),.)  -> snp_lm_result
  
  anova(snp_lm_result) -> snp_anova_result
  #snp_pvalue  <- formatC(snp_anova_result$`Pr(>F)`[51],format="e",digits=3)
  total_popstr_PVE <- PVE %>% filter(trait == chosen_trait, type == "population_structure") %>% pull(value) 
  total_snp_PVE <- round(sum(snp_anova_result$`Sum Sq`[1:(length(snp_anova_result$`Sum Sq`) - 1)])/
                           sum(snp_anova_result$`Sum Sq`) * 100,2) - total_popstr_PVE
  
  #Check for uncontrolled PVE
  snp_data_reformatted %>% 
    filter(snp_id %in% chosen_snps) %>%
    mutate(snp_id = paste("snp_",as.numeric(as.factor(snp_id)),sep="")) %>%
    spread(snp_id, genotype) %>%
    inner_join(.,pcavalues) %>%
    inner_join(.,phenotypes %>% select(sample, chosen_phenotype)) %>%
    inner_join(.,sample_info) %>%
    #select(sample,population, genotype,PC1:PC614, chosen_phenotype) %>% 
    lm(as.formula(paste("chosen_phenotype ~ ",
                        paste(paste("snp_",1:100,sep=""),collapse=" + "), 
                        sep="")),.)  -> snp_lm_result
  
  anova(snp_lm_result) -> snp_anova_result
  #snp_pvalue  <- formatC(snp_anova_result$`Pr(>F)`[51],format="e",digits=3)
  total_snp_PVE_uncontrolled <- round(sum(snp_anova_result$`Sum Sq`[1:(length(snp_anova_result$`Sum Sq`) - 1)])/
                                        sum(snp_anova_result$`Sum Sq`) * 100,2) 
  
  PVE_tmp <- tibble(type = c("snps","snps"),
                    value = c(total_snp_PVE,total_snp_PVE_uncontrolled),
                    trait = c(chosen_trait,chosen_trait),
                    controlled = c("controlled","uncontrolled"))
  
  PVE <- rbind(PVE, PVE_tmp)
}

pdf(paste("gwas/",prefix,"pve.v1.pdf",sep=""),width=12,height=8)
PVE %>% 
  ggplot(.,aes(x=trait,y=value,fill=type)) + 
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  ylab("PVE") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  scale_y_continuous(breaks=seq(0,100,10)) +
  facet_wrap(~controlled)
dev.off()
############
#Scraps
############
ref_grid(lm_result)
lm_means <- emmeans(lm_result, "triangle_genotype") 
lm_confint <- confint(lm_means)
tibble(mean = (lm_confint$emmean*sd_phenotype)+mean_phenotype,
       lower_ci = (lm_confint$lower.CL*sd_phenotype)+mean_phenotype,
       upper_ci = (lm_confint$upper.CL*sd_phenotype)+mean_phenotype,
       genotype = (lm_confint$triangle_genotype)) %>%
  ggplot(.,aes(x=genotype,y=mean)) + geom_point() +
  geom_errorbar(aes(x=genotype,ymin=lower_ci,ymax=upper_ci)) +
  theme_bw() +
  ylab(paste(chosen_trait)) +
  ggtitle(paste(chosen_chr,"-",chosen_mds,"\n",
                "p-value = ",inv_pvalue,"| PVE = ",PVE,"%",sep=""))



ref_grid(snp_lm_result)

snp_lm_means <- emmeans(snp_lm_result, "genotype") 
snp_lm_confint <- confint(snp_lm_means)
tibble(mean = (snp_lm_confint$emmean*sd_phenotype)+mean_phenotype,
       lower_ci = (snp_lm_confint$lower.CL*sd_phenotype)+mean_phenotype,
       upper_ci = (snp_lm_confint$upper.CL*sd_phenotype)+mean_phenotype,
       genotype = (snp_lm_confint$genotype)) %>%
  ggplot(.,aes(x=genotype,y=mean)) + geom_point() +
  geom_errorbar(aes(x=genotype,ymin=lower_ci,ymax=upper_ci)) +
  theme_bw() +
  ylab(paste(chosen_trait)) +
  ggtitle(paste(chosen_chr,"-",chosen_mds,"\n",
                "p-value = ",inv_pvalue,"| PVE = ",PVE,"%",sep=""))



