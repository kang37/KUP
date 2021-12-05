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

# rank abundance curves at land use level
png(filename = "Out_rac_land_use.png", 
    width = 2200, height = 1500, 
    type = c("cairo"), res = 300)
ggarrange(fun_rank_plot(lu_tree_rank, "(a)", method = "land_use"),
          fun_rank_plot(lu_shrub_rank, "(b)", method = "land_use"), 
          nrow = 2, common.legend = TRUE, legend = "bottom")
dev.off()

# richness ~ land use at quadrat level
png(filename = "Out_richness_quadrat.png", 
    width = 1500, height = 900, 
    type = c("cairo"), res = 300)
ggplot(qua_plant_div) + 
  geom_boxplot(aes(Land_use_type, Richness)) + 
  labs(x = "Land use type", y = "Quadrat richness") + 
  geom_text(data = data.frame(
    Land_use_type = Land_use_type_faclev, 
    Label = c(rep("a", 4), "ab", "b")
  ), aes(x = Land_use_type, y = Inf, label = Label), vjust = 2) + 
  scale_y_continuous(limits = c(0, max(qua_plant_div$Richness)), 
                     expand = expansion(mult = c(0, 0.1))) + 
  theme_bw()
dev.off()

# indexes ~ land use at quadrat level 
png(filename = "Out_indexes_quadrat.png", 
    width = 2400, height = 1500, 
    type = c("cairo"), res = 300)
ggarrange(fun_box_plot(qua_tree_div_long, tree_box_pvalue, "(a)"), 
          fun_box_plot(qua_shrub_div_long, shrub_box_pvalue, "(b)"))
dev.off()

# nMDS at quadrat level 
png(filename = "Out_nmds_quadrat.png", 
    width = 2400, height = 1500, 
    type = c("cairo"), res = 300)
ggarrange(
  fun_nmds_plot(tree_mds_selected, tree_hulls, "(a)", 
                tree_mds_meta, tree_anosim), 
  fun_nmds_plot(shrub_mds_selected, shrub_hulls, "(b)", 
                shrub_mds_meta, shrub_anosim),
  common.legend = T, legend = "right"
)
dev.off()

