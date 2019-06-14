#For plotting the inversion Dstat measures
library(tidyverse)
library(gridExtra)
info <- read_tsv("sample_info_apr_2018.tsv")
files <- list.files("Dstat", "inversion*.*txt") 

pdf("Dstat/inversion_dstat_tmp.pdf")
for (i in 1:length(files)){
  file_chosen <- files[i]
  
  dstat <- read_tsv(paste0("Dstat/",file_chosen),
                    col_names = c("chr","mds","inv_species","sample1","sample2","name","genome_d","genome_d_var","inv_d",
                                  "inv_abba","inv_baba"))
  

  
  if (dstat$inv_species[1] == "petiolaris"){
    print(
    dstat %>%
      mutate(all_three = paste0(sample1,sample2,name)) %>%
      distinct(all_three, .keep_all = TRUE) %>%
      gather(.,type,d, c(genome_d,inv_d)) %>% 
      inner_join(info) %>% 
      filter(species != "PetPet",species != "PetFal") %>% 
      ggplot(.,aes(x=species,y=d,fill=type)) + geom_boxplot(position="dodge") +
      theme_bw() + geom_hline(yintercept = 0,linetype="dotted") +
      ggtitle(paste(dstat$chr[1],dstat$mds[1],dstat$inv_species[1])) +
      scale_fill_brewer(palette = "Set1",
                        name = "Region",
                        labels = c("Genome wide","Inversion"))
    )

    
  }else{
    print(
    dstat %>%
      mutate(all_three = paste0(sample1,sample2,name)) %>%
      distinct(all_three, .keep_all = TRUE) %>%
      gather(.,type,d, c(genome_d,inv_d)) %>% 
      inner_join(info) %>% 
      #filter(species != "PetPet",species != "PetFal") %>% 
      ggplot(.,aes(x=species,y=d,fill=type)) + geom_boxplot(position="dodge") +
      theme_bw() + geom_hline(yintercept = 0,linetype="dotted") +
      ggtitle(paste(dstat$chr[1],dstat$mds[1],dstat$inv_species[1])) +
      scale_fill_brewer(palette = "Set1",
                        name = "Region",
                        labels = c("Genome wide","Inversion"))
    )
  }
  
}
dev.off()
