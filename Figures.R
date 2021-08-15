# extrapolation of richness at city level
png(filename = "Out_sac_city.png", 
    width = 1500, height = 900, 
    type = c("cairo"), res = 300)
fun_accum(all_plant_data, 348, 350, method = "city") + 
  theme(legend.position = "none")
dev.off()

# rank abundance curves at city level 
png(filename = "Out_rac_city.png", 
    width = 1500, height = 900, 
    type = c("cairo"), res = 300)
ggarrange(plotlist = list(
  fun_rank_plot(city_tree_rank, title = "(a)", method = "city"), 
  fun_rank_plot(city_shrub_rank, title = "(b)", method = "city")
), nrow = 1, common.legend = TRUE)
dev.off()

# extrapolation of richness at land use level 
png(filename = "Out_sac_land_use.png", 
    width = 1500, height = 900, 
    type = c("cairo"), res = 300)
fun_accum(all_plant_data, 348, 50, method = "land_use") + 
  labs(color = "Land use") + 
  scale_color_manual(values = c("#FF0000", "#FF7800", "#DF73FF", 
                                "#BFBF30", "#6BE400", "#00733E"))
dev.off()
