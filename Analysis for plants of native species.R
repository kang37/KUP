## general description
# the number of species
cat("total species:", length(unique(subset(all_plant_data, Nt_ex == "nt")$Species_CN)), "\n", 
    "total genera:", length(unique(subset(all_plant_data, Nt_ex == "nt")$Genus)), "\n", 
    "total families:", length(unique(subset(all_plant_data, Nt_ex == "nt")$Family)), "\n", "\n", 
    "total species of trees:", length(unique(tree_native_data$Species_CN)), "\n", 
    "total species of shrubs:", length(unique(shrub_native_data$Species_CN)), "\n", 
    "common species of trees and shrubs:", length(intersect(unique(tree_native_data$Species_CN), 
                                                            unique(shrub_native_data$Species_CN))), "\n", 
    "species solely for trees:", length(setdiff(unique(tree_native_data$Species_CN), 
                                                unique(shrub_native_data$Species_CN))), "\n", 
    "Species solely for shrubs:", length(setdiff(unique(shrub_native_data$Species_CN), 
                                                 unique(tree_native_data$Species_CN))))

# the number of trees or area of shrubs
cat("number of trees:", nrow(tree_native_data), "\n", 
    "number of tree-plot:", nrow(tree_native_diversity), "\n", 
    "area of shrubs:", sum(shrub_native_data$Area), "\n", 
    "number of shrub-plot:", nrow(shrub_native_diversity))

# top species of trees by number
tree_native_data %>% 
  group_by(Family) %>% dplyr::summarise(Number = sum(Stem), Prop = Number/sum(tree_native_data$Stem)) %>% 
  arrange(desc(Prop))

# top species of shrubs by area
shrub_native_data %>% 
  group_by(Family) %>% dplyr::summarise(SArea = sum(Area), Prop = SArea/sum(shrub_native_data$Area)) %>% 
  arrange(desc(Prop))

