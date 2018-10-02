pdf("All_species_sampling_map.pdf",height=8,width=8)
ggplot(target_state, aes(long, lat)) +
  geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
  coord_quickmap() +
  geom_point(data=pop_loc,
                  aes(x=long, y=lat,color=taxon),size=4) +
  theme_bw() +
  scale_color_brewer(palette = "Set1",name="Species") +
  xlab("Longitude") + ylab("Latitude") +
  scale_x_continuous(limits = long_range, expand = c(0, 0)) +
  scale_y_continuous(limits = lat_range, expand = c(0, 0)) 
dev.off()
