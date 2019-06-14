library(tidyverse)

inv.div <- read_tsv("divergence_time/betweenHap.tsv")
spe.div <- read_tsv("divergence_time/betweenSpecies.tsv")

scale <- 6.9e-9 * 1000000
pdf("Inversion_dates.v1.pdf")
inv.div %>%
  mutate(inv_id = paste0(species,".",gsub("Ha412HO","",inversion))) %>%
  ggplot(.) + 
  geom_point(aes(x=fct_reorder(inv_id,height),y=height/scale,
                 color=species)) +
  geom_rect(aes(xmin=0,xmax=100,ymin=spe.div$height_HPD95_lower[1]/scale,
                              ymax=spe.div$height_HPD95_upper[1]/scale),fill="grey") +
  geom_rect(aes(xmin=0,xmax=100,ymin=spe.div$height_HPD95_lower[2]/scale,
                ymax=spe.div$height_HPD95_upper[2]/scale),fill="grey") +
  geom_rect(aes(xmin=0,xmax=100,ymin=spe.div$height_HPD95_lower[3]/scale,
                ymax=spe.div$height_HPD95_upper[3]/scale),fill="grey") +
  geom_point(aes(x=fct_reorder(inv_id,height),y=height/scale,
                             color=species)) +
  geom_linerange(aes(x=inv_id,ymin=height_HPD95_lower/scale,ymax=height_HPD95_upper/scale,
                    color=species)) +
  coord_cartesian(xlim=c(0,nrow(inv.div))) +
  scale_color_brewer(palette = "Set1",name="Species") +
  annotate("text", x = 30, y = (spe.div$height[1]/scale), label = spe.div$species[1],
           hjust = 0.5) +
  annotate("text", x = 10, y = (spe.div$height[2]/scale), label = spe.div$species[2],
           hjust = 0.5) +
  annotate("text", x = 20, y = (spe.div$height[3]/scale), label = spe.div$species[3],
           hjust = 0.5) +
  theme_bw() +
  ylab("MYA") +
  xlab("LNH") +
  theme(axis.text.x=element_text(angle=60, hjust=1))
dev.off()
