setwd("C:/Rdata/KUP")
library(tidyverse)
library(Rmisc)
library(vegan)
library(PerformanceAnalytics)
library(digest)
library(car)
library(dunn.test)
library(leaps)
library(gvlma)
library(BiodiversityR)
opar <- par(no.readonly = TRUE)

# define the factor levels
Ward_faclev <- c("Ukyo-ku", "Sakyo-ku", "Kita-ku", "Kamigyo-ku", "Nakagyo-ku", "Shimogyo-ku", "Higashiyama-ku", "Yamashina-ku", "Fushimi-ku", "Minami-ku", "Nishikyo-ku")
Landuse_class_faclev <- c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")
Land_ownership_faclev <- c("Private area", "Propriater area", "Public area")

## get data
# data of all_plot_info
all_plot_info <- read.csv("In_plot_info.csv", stringsAsFactors = FALSE) %>%
  mutate(Dist_level = NA)
all_plot_info$Dist_level[all_plot_info$Dist < 2000] <- "1"
all_plot_info$Dist_level[all_plot_info$Dist >= 2000 & all_plot_info$Dist < 4000] <- "2"
all_plot_info$Dist_level[all_plot_info$Dist >= 4000 & all_plot_info$Dist < 6000] <- "3"
all_plot_info$Dist_level[all_plot_info$Dist >= 6000 & all_plot_info$Dist < 8000] <- "4"
all_plot_info$Dist_level[all_plot_info$Dist >= 8000] <- "5"
all_plot_info <- all_plot_info %>% 
  mutate(Dist_level = factor(Dist_level, levels = c("1", "2", "3", "4", "5")), 
         Ward = factor(Ward, levels = Ward_faclev), 
         Landuse_class = factor(Landuse_class, levels = Landuse_class_faclev))
all_plot_info$Land_ownership <- NA
all_plot_info$Land_ownership[all_plot_info$Landuse_detail %in% c("道路", "高层商业", "公园", "河流", "机构", "政府团地")] <- "Public area"
all_plot_info$Land_ownership[all_plot_info$Landuse_detail %in% c("低层私宅", "私宅")] <- "Private area"
all_plot_info$Land_ownership[all_plot_info$Landuse_detail %in% c("低层公寓", "高层公寓", "低层商业", "工厂", "农田", "其他", "球场", "寺庙神社", "学校")] <- "Propriater area"

# data of all_plant_info
all_plant_info <- read.csv("In_plant_info.csv", stringsAsFactors = FALSE)
all_plant_data <- read.csv("In_plant_data.csv", stringsAsFactors = FALSE) %>%
  left_join(all_plant_info, by = "Species_CN") %>%
  left_join(all_plot_info,by = "Plot_ID")
all_plant_data$Landuse_subclass <- NULL
all_plant_data$Landuse_agg <- NULL

# data of trees
tree_data <- all_plant_data %>% 
  subset(Tree_shrub == "Tree")

# data of shrubs
shrub_data <- all_plant_data %>% 
  subset(Tree_shrub == "Shrub")

# data of native trees
tree_native_data <- subset(tree_data, Nt_ex == "nt")

# shrub_native_data
shrub_native_data <- subset(shrub_data, Nt_ex == "nt")

# data of tree_plant_info
tree_plant_info <- subset(all_plant_data, Tree_shrub == "Tree", select = "Species_CN") %>% 
  unique() %>% 
  left_join(all_plant_info, by = "Species_CN")

# data of shrub info
shrub_plant_info <- subset(all_plant_data, Tree_shrub == "Shrub", select = "Species_CN") %>% 
  unique() %>% 
  left_join(all_plant_info, by = "Species_CN")

# data of tree_diversity
tree_diversity <- subset(tree_data, select = c("Plot_ID", "Species_CN", "Stem")) %>%
  pivot_wider(names_from = Species_CN, values_from = Stem, 
              values_fn = list(Stem = sum), values_fill = list(Stem = 0)) %>% 
  mutate(Sum_stem = rowSums(.[2:ncol(.)]), 
         Richness = apply(.[2:ncol(.)]>0, 1, sum),
         Shannon = diversity(.[2:ncol(.)], index = "shannon"), 
         Simpson = diversity(.[2:ncol(.)], index = "simpson"),
         Evenness = Shannon / log(Richness)) %>%
  left_join(all_plot_info,by = "Plot_ID") %>% 
  as.data.frame()

tree_diversity_perc_planted <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pla_spo == "Planted", Stem, 0)/sum(Stem))) 
tree_diversity_perc_nonpot <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pot == "Non_pot", Stem, 0)/sum(Stem)))
tree_diversity_perc_private <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pub_pri == "Private", Stem, 0)/sum(Stem)))
tree_diversity_perc_nonstreet <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Street == "Non_street", Stem, 0)/sum(Stem)))
tree_diversity_perc_native <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Nt_ex == "nt", Stem, 0)/sum(Stem)))

tree_diversity <- tree_diversity %>% mutate(
  perc_planted = tree_diversity_perc_planted$perc, 
  perc_nonpot = tree_diversity_perc_nonpot$perc, 
  perc_private = tree_diversity_perc_private$perc, 
  perc_nonstreet = tree_diversity_perc_nonstreet$perc, 
  perc_native = tree_diversity_perc_native$perc
)
rm(tree_diversity_perc_planted, tree_diversity_perc_nonpot, tree_diversity_perc_private, tree_diversity_perc_nonstreet, tree_diversity_perc_native)

# tree_native_diversity
tree_native_diversity <- subset(tree_native_data, select = c("Plot_ID", "Species_CN", "Stem")) %>%
  pivot_wider(names_from = Species_CN, values_from = Stem, 
              values_fn = list(Stem = sum), values_fill = list(Stem = 0)) %>% 
  mutate(Sum_stem = rowSums(.[2:ncol(.)]), 
         Richness = apply(.[2:ncol(.)]>0, 1, sum),
         Shannon = diversity(.[2:ncol(.)], index = "shannon"), 
         Simpson = diversity(.[2:ncol(.)], index = "simpson"),
         Evenness = Shannon / log(Richness)) %>%
  left_join(all_plot_info,by = "Plot_ID") %>% 
  as.data.frame()

tree_native_diversity_perc_planted <- tree_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pla_spo == "Planted", Stem, 0)/sum(Stem))) 
tree_native_diversity_perc_nonpot <- tree_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pot == "Non_pot", Stem, 0)/sum(Stem)))
tree_native_diversity_perc_private <- tree_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pub_pri == "Private", Stem, 0)/sum(Stem)))
tree_native_diversity_perc_nonstreet <- tree_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Street == "Non_street", Stem, 0)/sum(Stem)))
tree_native_diversity_perc_native <- tree_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Nt_ex == "nt", Stem, 0)/sum(Stem)))

tree_native_diversity <- tree_native_diversity %>% mutate(
  perc_planted = tree_native_diversity_perc_planted$perc, 
  perc_nonpot = tree_native_diversity_perc_nonpot$perc, 
  perc_private = tree_native_diversity_perc_private$perc, 
  perc_nonstreet = tree_native_diversity_perc_nonstreet$perc, 
  perc_native = tree_native_diversity_perc_native$perc
)
rm(tree_native_diversity_perc_planted, tree_native_diversity_perc_nonpot, tree_native_diversity_perc_private, tree_native_diversity_perc_nonstreet, tree_native_diversity_perc_native)

# data of shrub_diversity
shrub_diversity <- subset(shrub_data, select = c("Plot_ID", "Species_CN", "Area")) %>%
  pivot_wider(names_from = Species_CN, values_from = Area, 
              values_fn = list(Area = sum), values_fill = list(Area = 0)) %>% 
  mutate(Sum_area = rowSums(.[2:ncol(.)]), 
         Richness = apply(.[2:ncol(.)]>0, 1, sum),
         Shannon = diversity(.[2:ncol(.)], index = "shannon"), 
         Simpson = diversity(.[2:ncol(.)], index = "simpson"),
         Evenness = Shannon / log(Richness)) %>%
  left_join(all_plot_info,by = "Plot_ID") %>%
  as.data.frame()

shrub_diversity_perc_planted <- shrub_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pla_spo == "Planted", Area, 0)/sum(Area))) 
shrub_diversity_perc_nonpot <- shrub_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pot == "Non_pot", Area, 0)/sum(Area)))
shrub_diversity_perc_private <- shrub_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pub_pri == "Private", Area, 0)/sum(Area)))
shrub_diversity_perc_nonstreet <- shrub_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Street == "Non_street", Area, 0)/sum(Area)))
shrub_diversity_perc_native <- shrub_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Nt_ex == "nt", Area, 0)/sum(Area)))

shrub_diversity <- shrub_diversity %>% mutate(
  perc_planted = shrub_diversity_perc_planted$perc, 
  perc_nonpot = shrub_diversity_perc_nonpot$perc, 
  perc_private = shrub_diversity_perc_private$perc, 
  perc_nonstreet = shrub_diversity_perc_nonstreet$perc, 
  perc_native = shrub_diversity_perc_native$perc
)
rm(shrub_diversity_perc_planted, shrub_diversity_perc_nonpot, shrub_diversity_perc_private, shrub_diversity_perc_nonstreet, shrub_diversity_perc_native)

# shrub_native_diversity
shrub_native_diversity <- subset(shrub_native_data, select = c("Plot_ID", "Species_CN", "Area")) %>%
  pivot_wider(names_from = Species_CN, values_from = Area, 
              values_fn = list(Area = sum), values_fill = list(Area = 0)) %>% 
  mutate(Sum_area = rowSums(.[2:ncol(.)]), 
         Richness = apply(.[2:ncol(.)]>0, 1, sum),
         Shannon = diversity(.[2:ncol(.)], index = "shannon"), 
         Simpson = diversity(.[2:ncol(.)], index = "simpson"),
         Evenness = Shannon / log(Richness)) %>%
  left_join(all_plot_info,by = "Plot_ID") %>% 
  as.data.frame()

shrub_native_diversity_perc_planted <- shrub_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pla_spo == "Planted", Area, 0)/sum(Area))) 
shrub_native_diversity_perc_nonpot <- shrub_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pot == "Non_pot", Area, 0)/sum(Area)))
shrub_native_diversity_perc_private <- shrub_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pub_pri == "Private", Area, 0)/sum(Area)))
shrub_native_diversity_perc_nonstreet <- shrub_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Street == "Non_street", Area, 0)/sum(Area)))
shrub_native_diversity_perc_native <- shrub_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Nt_ex == "nt", Area, 0)/sum(Area)))

shrub_native_diversity <- shrub_native_diversity %>% mutate(
  perc_planted = shrub_native_diversity_perc_planted$perc, 
  perc_nonpot = shrub_native_diversity_perc_nonpot$perc, 
  perc_private = shrub_native_diversity_perc_private$perc, 
  perc_nonstreet = shrub_native_diversity_perc_nonstreet$perc, 
  perc_native = shrub_native_diversity_perc_native$perc
)
rm(shrub_native_diversity_perc_planted, shrub_native_diversity_perc_nonpot, shrub_native_diversity_perc_private, shrub_native_diversity_perc_nonstreet, shrub_native_diversity_perc_native)

# tree diversity longer and shrub diversity longer dataset
tree_diversity_long <- 
  subset(tree_diversity, select = c("Sum_stem", "Richness", "Shannon", "Evenness", "Landuse_class", "Land_ownership", "Dist")) %>% 
  pivot_longer(cols = c("Sum_stem", "Richness", "Shannon", "Evenness"), 
               names_to = "index", values_to = "index_value") %>% 
  pivot_longer(cols = c("Landuse_class", "Land_ownership"), 
               names_to = "attr", values_to = "attr_value") %>% 
  mutate(index = factor(index, levels = c("Sum_stem", "Richness", "Shannon", "Evenness")), 
         attr = factor(attr, levels = c("Landuse_class", "Land_ownership")), 
         attr_value = factor(attr_value, levels = c(Landuse_class_faclev, Land_ownership_faclev)))

shrub_diversity_long <- 
  subset(shrub_diversity, select = c("Sum_area", "Richness", "Shannon", "Evenness", "Landuse_class", "Land_ownership", "Dist")) %>% 
  pivot_longer(cols = c("Sum_area", "Richness", "Shannon", "Evenness"), 
               names_to = "index", values_to = "index_value") %>% 
  pivot_longer(cols = c("Landuse_class", "Land_ownership"), 
               names_to = "attr", values_to = "attr_value") %>% 
  mutate(index = factor(index, levels = c("Sum_area", "Richness", "Shannon","Evenness")), 
         attr = factor(attr, levels = c("Landuse_class", "Land_ownership")), 
         attr_value = factor(attr_value, levels = c(Landuse_class_faclev, Land_ownership_faclev)))



### analysis begins



## general description
# the number of species
cat("total species:", length(unique(all_plant_data$Species_CN)), "\n", 
    "total genera:", length(unique(all_plant_data$Genus)), "\n", 
    "total families:", length(unique(all_plant_data$Family)), "\n", "\n", 
    "total species of trees:", length(unique(tree_data$Species_CN)), "\n", 
    "total species of shrubs:", length(unique(shrub_data$Species_CN)), "\n", 
    "common species of trees and shrubs:", length(intersect(unique(tree_data$Species_CN), 
                                                            unique(shrub_data$Species_CN))), "\n", 
    "species solely for trees:", length(setdiff(unique(tree_data$Species_CN), 
                                               unique(shrub_data$Species_CN))), "\n", 
    "Species solely for shrubs:", length(setdiff(unique(shrub_data$Species_CN), 
                                                 unique(tree_data$Species_CN))))

# the number of trees or area of shrubs
cat("number of trees:", nrow(tree_data), "\n", 
    "number of tree-plot:", nrow(tree_diversity), "\n", 
    "area of shrubs:", sum(shrub_data$Area), "\n", 
    "number of shrub-plot:", nrow(shrub_diversity))

# top species of trees by number
tree_data %>% 
  group_by(Family) %>% dplyr::summarise(Number = sum(Stem), Prop = Number/sum(tree_data$Stem)) %>% 
  arrange(desc(Prop))

# top species of shrubs by area
shrub_data %>% 
  group_by(Family) %>% dplyr::summarise(SArea = sum(Area), Prop = SArea/sum(shrub_data$Area)) %>% 
  arrange(desc(Prop))



## attributes of the tree and shrub
# the exotic vs. native regarding species
all_plant_info %>% group_by(Nt_ex) %>% 
  dplyr::summarise(n()/nrow(all_plant_info))

# while regarding the number or area
tree_data %>% group_by(Nt_ex) %>% 
  dplyr::summarise(n()/sum(tree_data$Stem))
shrub_data %>% group_by(Nt_ex) %>% 
  dplyr::summarise(sum(Area)/sum(shrub_data$Area))

# the graph of attributes of the plants
par(mfrow= c(2,3), cex.axis = 1.5)
j <- 0
for (i in c("Pla_spo", "Pub_pri", "Nt_ex")) {
  j <- j + 1
  barplot(tapply(subset(all_plant_data, Tree_shrub == "Tree")[, "Stem"], 
                 subset(all_plant_data, Tree_shrub == "Tree")[, i], 
                 sum), ylim = c(0, 1400))
  title(main = paste("(", letters[j], ")"), adj = 0)
}
for (i in c("Pla_spo", "Pub_pri", "Nt_ex")) {
  j <- j + 1
  barplot(tapply(subset(all_plant_data, Tree_shrub == "Shrub")[, "Area"], 
                 subset(all_plant_data, Tree_shrub == "Shrub")[, i], 
                 sum), ylim = c(0, 1400))
  title(main = paste("(", letters[j], ")"), adj = 0)
}
par(opar)



## Species accumulation curve 
par(mfrow = c(1,2))
# Species accumulation curve for trees
plot(specaccum(subset(tree_diversity, Landuse_class == "Com")[,2:143]), col = "red", lty = 1, 
     ci.lty = 0, xlim = c(0, 40), ylim = c(0, 110))
plot(specaccum(subset(tree_diversity, Landuse_class == "Com neigh")[,2:143]),col = "red", lty = 2, 
     ci.lty = 0, add = T)
plot(specaccum(subset(tree_diversity, Landuse_class == "R low")[,2:143]), col = "blue", lty = 1, 
     ci.lty = 0, add = T)
plot(specaccum(subset(tree_diversity, Landuse_class == "R high")[,2:143]), col = "blue", lty = 2, 
     ci.lty = 0, add = T)
plot(specaccum(subset(tree_diversity, Landuse_class == "R resi")[,2:143]), col = "blue", lty = 3, 
     ci.lty = 0, add = T)
plot(specaccum(subset(tree_diversity, Landuse_class == "Ind")[,2:143]), col = "black", lty = 1, 
     ci.lty = 0, add = T)
# Species accumulation curve for shrubs
plot(specaccum(subset(shrub_diversity, Landuse_class == "Com")[,2:143]), col = "red", lty = 1, 
     ci.lty = 0, xlim = c(0, 40), ylim = c(0, 110))
plot(specaccum(subset(shrub_diversity, Landuse_class == "Com neigh")[,2:143]),col = "red", lty = 2, 
     ci.lty = 0, add = T)
plot(specaccum(subset(shrub_diversity, Landuse_class == "R low")[,2:143]), col = "blue", lty = 1, 
     ci.lty = 0, add = T)
plot(specaccum(subset(shrub_diversity, Landuse_class == "R high")[,2:143]), col = "blue", lty = 2, 
     ci.lty = 0, add = T)
plot(specaccum(subset(shrub_diversity, Landuse_class == "R resi")[,2:143]), col = "blue", lty = 3, 
     ci.lty = 0, add = T)
plot(specaccum(subset(shrub_diversity, Landuse_class == "Ind")[,2:143]), col = "black", lty = 1, 
     ci.lty = 0, add = T)
par(opar)



## rank ahundance plot
# the rank abundance plot by land use types of trees: doesn't omit site 279
tree_rankabun_list <- vector("list", 5)
for (i in c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")) {
  tree_rankabun_ori <- as.data.frame(rankabundance(subset(tree_diversity, Landuse_class == i)[2:143]))
  tree_rankabun_list[[1]] <- c(tree_rankabun_list[[1]], rownames(tree_rankabun_ori))
  tree_rankabun_list[[2]] <- c(tree_rankabun_list[[2]], tree_rankabun_ori$rank)
  tree_rankabun_list[[3]] <- c(tree_rankabun_list[[3]], tree_rankabun_ori$abundance)
  tree_rankabun_list[[4]] <- c(tree_rankabun_list[[4]], tree_rankabun_ori$proportion)
  tree_rankabun_list[[5]] <- c(tree_rankabun_list[[5]], rep(i, nrow(tree_rankabun_ori)))
}
tree_rankabun_df <- data.frame(
  Species_CN = tree_rankabun_list[[1]], 
  rank = tree_rankabun_list[[2]], 
  abundance = tree_rankabun_list[[3]], 
  proportion = tree_rankabun_list[[4]], 
  Landuse_class = tree_rankabun_list[[5]]
) %>% 
  left_join(all_plant_info[, c("Species_CN", "Nt_ex")], by = "Species_CN")
tree_rankabun_df[tree_rankabun_df == 0] <- NA 

# the rank abundance plot by land use types of shrubs
shrub_rankabun_list <- vector("list", 5)
for (i in c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")) {
  shrub_rankabun_ori <- as.data.frame(rankabundance(subset(shrub_diversity, Landuse_class == i)[2:143]))
  shrub_rankabun_list[[1]] <- c(shrub_rankabun_list[[1]], rownames(shrub_rankabun_ori))
  shrub_rankabun_list[[2]] <- c(shrub_rankabun_list[[2]], shrub_rankabun_ori$rank)
  shrub_rankabun_list[[3]] <- c(shrub_rankabun_list[[3]], shrub_rankabun_ori$abundance)
  shrub_rankabun_list[[4]] <- c(shrub_rankabun_list[[4]], shrub_rankabun_ori$proportion)
  shrub_rankabun_list[[5]] <- c(shrub_rankabun_list[[5]], rep(i, nrow(shrub_rankabun_ori)))
}
shrub_rankabun_df <- data.frame(
  Species_CN = shrub_rankabun_list[[1]], 
  rank = shrub_rankabun_list[[2]], 
  abundance = shrub_rankabun_list[[3]], 
  proportion = shrub_rankabun_list[[4]], 
  Landuse_class = shrub_rankabun_list[[5]]
) %>% 
  left_join(all_plant_info[, c("Species_CN", "Nt_ex")], by = "Species_CN")
shrub_rankabun_df[shrub_rankabun_df == 0] <- NA

# rearrange and plot
tree_rankabun_df <- tree_rankabun_df %>% 
  mutate(Landuse_class = factor(Landuse_class, levels = c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")))
shrub_rankabun_df <- shrub_rankabun_df %>% 
  mutate(Landuse_class = factor(Landuse_class, levels = c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")))
Rmisc::multiplot(plotlist = c(
  list(ggplot(tree_rankabun_df, aes(rank, proportion, label = Species_CN)) + 
         geom_line() + 
         geom_point(aes(color = Nt_ex), alpha = 0.3, size = 2) + 
         geom_text(aes(label = ifelse(rank<4, as.character(Species_CN), "")), hjust = -0.5, vjust = 0) + 
         facet_wrap(~Landuse_class, nrow = 1) + labs(title = "(a)")), 
  list(ggplot(shrub_rankabun_df, aes(rank, proportion, label = Species_CN)) + 
         geom_line() + 
         geom_point(aes(color = Nt_ex), alpha = 0.3, size = 2) + 
         geom_text(aes(label = ifelse(rank<4, as.character(Species_CN), "")), hjust = -0.5, vjust = 0) + 
         facet_wrap(~Landuse_class, nrow = 1) + labs(title = "(b)"))
))
rm(tree_rankabun_list, tree_rankabun_ori, tree_rankabun_df, 
   shrub_rankabun_list, shrub_rankabun_ori, shrub_rankabun_df)

## mds analysis
set.seed(1234)
#
# nmds calculation for tree
tree_mds_selected_ID <- tree_diversity$Plot_ID[!(tree_diversity$Plot_ID %in% c(214, 261, 313, 244, 67))]
tree_mds_selected <- tree_diversity %>% filter(Plot_ID %in% tree_mds_selected_ID)
tree_mds_metaMDS <- tree_mds_selected %>% 
  select(2:143) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
tree_mds_metaMDS$stress
stressplot(tree_mds_metaMDS)
tree_mds_selected <- cbind(tree_mds_selected, tree_mds_metaMDS$points)
#
# nmds calculation for shrub
shrub_mds_selected_ID <- shrub_diversity$Plot_ID[!(shrub_diversity$Plot_ID %in% c(269, 214, 75, 244, 164))]
shrub_mds_selected <- shrub_diversity %>% filter(Plot_ID %in% shrub_mds_selected_ID)
shrub_mds_metaMDS <- shrub_mds_selected %>% 
  select(2:196) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
shrub_mds_metaMDS$stress
stressplot(shrub_mds_metaMDS)
shrub_mds_selected <- cbind(shrub_mds_selected, shrub_mds_metaMDS$points)
#
# mds plot for trees and shrubs
# general Anosim of trees and shrubs
tree_landuseclass_anosim_result <- anosim(tree_mds_selected[2:143], tree_mds_selected$Landuse_class)
tree_landownership_anosim_result <- anosim(tree_mds_selected[2:143], tree_mds_selected$Land_ownership)
shrub_landuseclass_anosim_result <- anosim(shrub_mds_selected[2:196], shrub_mds_selected$Landuse_class)
shrub_landownership_anosim_result <- anosim(shrub_mds_selected[2:196], shrub_mds_selected$Land_ownership)
# get statistic results as labels for the mds plots
tree_landuseclass_mds_lab <- paste("stress=", round(tree_mds_metaMDS$stress, digits = 3), 
                                   ", R=", round(tree_landuseclass_anosim_result$statistic, digits = 3), 
                                   ", p=", round(tree_landuseclass_anosim_result$signif, digits = 3), 
                                   sep = "")
tree_landownership_mds_lab <- paste("stress=", round(tree_mds_metaMDS$stress, digits = 3), 
                                    ", R=", round(tree_landownership_anosim_result$statistic, digits = 3), 
                                    ", p=", round(tree_landownership_anosim_result$signif, digits = 3), 
                                    sep = "")
shrub_landuseclass_mds_lab <- paste("stress=", round(shrub_mds_metaMDS$stress, digits = 3), 
                                    ", R=", round(shrub_landuseclass_anosim_result$statistic, digits = 3), 
                                    ", p=", round(shrub_landuseclass_anosim_result$signif, digits = 3), 
                                    sep = "")
shrub_landownership_mds_lab <- paste("stress=", round(shrub_mds_metaMDS$stress, digits = 3), 
                                     ", R=", round(shrub_landownership_anosim_result$statistic, digits = 3), 
                                     ", p=", round(shrub_landownership_anosim_result$signif, digits = 3), 
                                     sep = "")
# mds plots for trees and shrubs by land use types and land ownership 
Rmisc::multiplot(plotlist = list(
  ggplot(tree_mds_selected, aes(MDS1, MDS2, color = Landuse_class)) + geom_point(alpha = 0.7) + 
    labs(title = "Tree - Land use type", subtitle = tree_landuseclass_mds_lab), 
  ggplot(tree_mds_selected, aes(MDS1, MDS2, color = Land_ownership)) + geom_point(alpha = 0.7) + 
    labs(title = "Tree - Land ownership", subtitle = tree_landownership_mds_lab), 
  ggplot(shrub_mds_selected, aes(MDS1, MDS2, color = Landuse_class)) + geom_point(alpha = 0.7) + 
    labs(title = "Shrub - Land use type", subtitle = shrub_landuseclass_mds_lab), 
  ggplot(shrub_mds_selected, aes(MDS1, MDS2, color = Land_ownership)) + geom_point(alpha = 0.7) + 
    labs(title = "Shrub - Land ownership", subtitle = shrub_landownership_mds_lab)
), layout = matrix(1:4, nrow = 2, byrow = T))
#
# pairwise result of ANOSIM of trees by landuse_class
tree_pair_anosim_list <- vector("list",3)
tree_pair_anosim_list[[1]] <- c(combn(levels(factor(tree_mds_selected$Landuse_class)),2)[1,])
tree_pair_anosim_list[[2]] <- c(combn(levels(factor(tree_mds_selected$Landuse_class)),2)[2,])
set.seed(1234)
for (i in 1:length(tree_pair_anosim_list[[1]])) {
  tree_mds_selected_sub <- subset(tree_mds_selected, Landuse_class == tree_pair_anosim_list[[1]][i] |
                                    Landuse_class == tree_pair_anosim_list[[2]][i])
  tree_pair_anosim_list[[3]] <- c(tree_pair_anosim_list[[3]], 
                                  anosim(tree_mds_selected_sub[2:143], 
                                         tree_mds_selected_sub$Landuse_class)$signif)
}
tree_pair_anosim_df <- data.frame(comp_1 = tree_pair_anosim_list[[1]], 
                                  comp_2 = tree_pair_anosim_list[[2]], 
                                  p = tree_pair_anosim_list[[3]]) %>% subset(p < 0.05) %>% print()

# pairwise result of ANOSIM of trees by landuse_ownership
tree_pair_anosim_list <- vector("list",3)
tree_pair_anosim_list[[1]] <- c(combn(levels(factor(tree_mds_selected$Land_ownership)),2)[1,])
tree_pair_anosim_list[[2]] <- c(combn(levels(factor(tree_mds_selected$Land_ownership)),2)[2,])
set.seed(1234)
for (i in 1:length(tree_pair_anosim_list[[1]])) {
  tree_mds_selected_sub <- subset(tree_mds_selected, Land_ownership == tree_pair_anosim_list[[1]][i] |
                                    Land_ownership == tree_pair_anosim_list[[2]][i])
  tree_pair_anosim_list[[3]] <- c(tree_pair_anosim_list[[3]], 
                                  anosim(tree_mds_selected_sub[2:143], 
                                         tree_mds_selected_sub$Land_ownership)$signif)
}
tree_pair_anosim_df <- data.frame(comp_1 = tree_pair_anosim_list[[1]], 
                                  comp_2 = tree_pair_anosim_list[[2]], 
                                  p = tree_pair_anosim_list[[3]]) %>% subset(p < 0.05) %>% print()

# pairwise result of ANOSIM of shrubs by landuse_class
shrub_pair_anosim_list <- vector("list",3)
shrub_pair_anosim_list[[1]] <- c(combn(levels(factor(shrub_mds_selected$Landuse_class)),2)[1,])
shrub_pair_anosim_list[[2]] <- c(combn(levels(factor(shrub_mds_selected$Landuse_class)),2)[2,])
set.seed(1234)
for (i in 1:length(shrub_pair_anosim_list[[1]])) {
  shrub_mds_selected_sub <- subset(shrub_mds_selected, Landuse_class == shrub_pair_anosim_list[[1]][i] |
                                     Landuse_class == shrub_pair_anosim_list[[2]][i])
  shrub_pair_anosim_list[[3]] <- c(shrub_pair_anosim_list[[3]], 
                                   anosim(shrub_mds_selected_sub[2:196], 
                                          shrub_mds_selected_sub$Landuse_class)$signif)
}
shrub_pair_anosim_df <- data.frame(comp_1 = shrub_pair_anosim_list[[1]], 
                                   comp_2 = shrub_pair_anosim_list[[2]], 
                                   p = shrub_pair_anosim_list[[3]]) %>% subset(p < 0.05) %>% print()

# pairwise result of ANOSIM of shrubs by landuse_ownership
shrub_pair_anosim_list <- vector("list",3)
shrub_pair_anosim_list[[1]] <- c(combn(levels(factor(shrub_mds_selected$Land_ownership)),2)[1,])
shrub_pair_anosim_list[[2]] <- c(combn(levels(factor(shrub_mds_selected$Land_ownership)),2)[2,])
set.seed(1234)
for (i in 1:length(shrub_pair_anosim_list[[1]])) {
  shrub_mds_selected_sub <- subset(shrub_mds_selected, Land_ownership == shrub_pair_anosim_list[[1]][i] |
                                     Land_ownership == shrub_pair_anosim_list[[2]][i])
  shrub_pair_anosim_list[[3]] <- c(shrub_pair_anosim_list[[3]], 
                                   anosim(shrub_mds_selected_sub[2:196], 
                                          shrub_mds_selected_sub$Land_ownership)$signif)
}
shrub_pair_anosim_df <- data.frame(comp_1 = shrub_pair_anosim_list[[1]], 
                                   comp_2 = shrub_pair_anosim_list[[2]], 
                                   p = shrub_pair_anosim_list[[3]]) %>% subset(p < 0.05) %>% print()



## cor among the indexes
chart.Correlation(subset(tree_diversity, select = c("Sum_stem", "Richness", "Shannon", "Simpson", "Evenness")))
chart.Correlation(subset(shrub_diversity, select = c("Sum_area", "Richness", "Shannon", "Simpson", "Evenness")))



## kruskal test & boxplot for trees
# p-value matrix for tree
set.seed(1234)
boxplot_list_index_attr <- vector("list", 2)
#
{
  pvalue_list <- vector("list", 3)
  for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
    for (j in c("Landuse_class", "Land_ownership")) {
      pvalue_list[[1]] <- c(pvalue_list[[1]], i)
      pvalue_list[[2]] <- c(pvalue_list[[2]], j)
      pvalue_list[[3]] <- c(pvalue_list[[3]], 
                            round(kruskal.test(tree_diversity[, i] ~ tree_diversity[, j])$p.value,digits = 3))
    }
  }
}

{
  pvalue <- data.frame(index = pvalue_list[[1]],
                       attr = pvalue_list[[2]],
                       pvalue = pvalue_list[[3]])
  pvalue$label <- NA
  pvalue$label[pvalue$pvalue >= 0.05] <- paste("p=", pvalue$pvalue[pvalue$pvalue>0.05], sep = "")
  pvalue$label[pvalue$pvalue < 0.05 & pvalue$pvalue >= 0.01] <- 
    paste("p=", pvalue$pvalue[pvalue$pvalue < 0.05 & pvalue$pvalue >= 0.01], "*", sep = "")
  pvalue$label[pvalue$pvalue < 0.01 & pvalue$pvalue >= 0.001] <- 
    paste("p=", pvalue$pvalue[pvalue$pvalue < 0.01 & pvalue$pvalue >= 0.001], "**", sep = "")
  pvalue$label[pvalue$pvalue < 0.001] <- 
    paste("p=", pvalue$pvalue[pvalue$pvalue < 0.001], "***", sep = "")
}

boxplot_list_index_attr[[1]] <- tree_diversity_long %>% 
  na.omit() %>%
  ggplot(aes(attr_value, index_value)) + geom_boxplot() + 
  facet_grid(index ~ attr, scales = "free", space = "free_x", switch = "both") + 
  scale_y_continuous(expand = expand_scale(mult = c(0.05,0.3))) +
  geom_text(data = pvalue, aes(x =Inf, y = Inf, label = label), size=3.5, hjust = 1.05, vjust = 1.5) +
  theme(axis.text = element_text(angle = 90)) + 
  labs(title = "(a)", x = NULL, y = NULL)
# for shrub
{
  pvalue_list <- vector("list", 3)
  for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
    for (j in c("Landuse_class", "Land_ownership")) {
      pvalue_list[[1]] <- c(pvalue_list[[1]], i)
      pvalue_list[[2]] <- c(pvalue_list[[2]], j)
      pvalue_list[[3]] <- c(pvalue_list[[3]], round(kruskal.test(shrub_diversity[, i] ~ shrub_diversity[, j])$p.value,digits = 3))
    }
  }
  }

{
  pvalue <- data.frame(index = pvalue_list[[1]],
                       attr = pvalue_list[[2]],
                       pvalue = pvalue_list[[3]])
  pvalue$label <- NA
  pvalue$label[pvalue$pvalue >= 0.05] <- paste("p=", pvalue$pvalue[pvalue$pvalue>0.05], sep = "")
  pvalue$label[pvalue$pvalue < 0.05 & pvalue$pvalue >= 0.01] <- 
    paste("p=", pvalue$pvalue[pvalue$pvalue < 0.05 & pvalue$pvalue >= 0.01], "*", sep = "")
  pvalue$label[pvalue$pvalue < 0.01 & pvalue$pvalue >= 0.001] <- 
    paste("p=", pvalue$pvalue[pvalue$pvalue < 0.01 & pvalue$pvalue >= 0.001], "**", sep = "")
  pvalue$label[pvalue$pvalue < 0.001] <- 
    paste("p=", pvalue$pvalue[pvalue$pvalue < 0.001], "***", sep = "")
}
#
boxplot_list_index_attr[[2]] <- shrub_diversity_long %>% 
  na.omit() %>%
  ggplot(aes(attr_value, index_value)) + geom_boxplot() + 
  facet_grid(index ~ attr, scales = "free", space = "free_x", switch = "both") + 
  scale_y_continuous(expand = expand_scale(mult = c(0.05,0.3))) +
  geom_text(data = pvalue, aes(x =Inf, y = Inf, label = label), size=3.5, hjust = 1.05, vjust = 1.5) +
  theme(axis.text = element_text(angle = 90)) + 
  labs(title = "(b)", x = NULL, y = NULL)
Rmisc::multiplot(plotlist = boxplot_list_index_attr, cols = 2)
# delete the vars
rm(pvalue, pvalue_list)



## pairwise dunn test of diversity ~ I(landuse class + ownership)
pairwise_list <- vector("list", 5)
# list of pairwise test of diversity and attrs of tree
for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
  for (j in c("Landuse_class", "Land_ownership")) {
    pairwise_list[[1]] <- c(pairwise_list[[1]], 
                            rep("tree", 
                                length(dunn.test(tree_diversity[, i], tree_diversity[, j])$P.adjusted)))
    pairwise_list[[2]] <- c(pairwise_list[[2]], 
                            rep(i, length(dunn.test(tree_diversity[, i], tree_diversity[, j])$P.adjusted)))
    pairwise_list[[3]] <- c(pairwise_list[[3]], 
                            rep(j, length(dunn.test(tree_diversity[, i], tree_diversity[, j])$P.adjusted)))
    pairwise_list[[4]] <- c(pairwise_list[[4]], 
                            dunn.test(tree_diversity[, i], tree_diversity[, j])$comparisons)
    pairwise_list[[5]] <- c(pairwise_list[[5]], 
                            dunn.test(tree_diversity[, i], tree_diversity[, j])$P.adjusted)
  }
}
# list of pairwise test of diversity and attrs of shrub
for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
  for (j in c("Landuse_class", "Land_ownership")) {
    pairwise_list[[1]] <- c(pairwise_list[[1]], 
                            rep("shrub", 
                                length(dunn.test(shrub_diversity[, i], shrub_diversity[, j])$P.adjusted)))
    pairwise_list[[2]] <- c(pairwise_list[[2]], 
                            rep(i, length(dunn.test(shrub_diversity[, i], shrub_diversity[, j])$P.adjusted)))
    pairwise_list[[3]] <- c(pairwise_list[[3]], 
                            rep(j, length(dunn.test(shrub_diversity[, i], shrub_diversity[, j])$P.adjusted)))
    pairwise_list[[4]] <- c(pairwise_list[[4]], 
                            dunn.test(shrub_diversity[, i], shrub_diversity[, j])$comparisons)
    pairwise_list[[5]] <- c(pairwise_list[[5]], 
                            dunn.test(shrub_diversity[, i], shrub_diversity[, j])$P.adjusted)
    
  }
}
# data frame of pairwise test of tree and shrub
pairwise_df_up <- data.frame(taxa = pairwise_list[[1]], 
                             index = pairwise_list[[2]], 
                             attr = pairwise_list[[3]], 
                             comparison = pairwise_list[[4]], 
                             p = pairwise_list[[5]]) %>% 
  separate(comparison, into = c("comparison_1", "comparison_2"), sep = " - ")
pairwise_df_down <- data.frame(
  taxa = pairwise_df_up$taxa, 
  index = pairwise_df_up$index, 
  attr = pairwise_df_up$attr,
  comparison_1 = pairwise_df_up$comparison_2, 
  comparison_2 = pairwise_df_up$comparison_1, 
  p = pairwise_df_up$p)
pairwise_df <- rbind(pairwise_df_up, pairwise_df_down) %>% 
  mutate(index = factor(index, levels = c("Sum_stem", "Sum_area", "Richness", "Shannon", "Evenness")), 
         comparison_1 = factor(comparison_1, levels = c(Landuse_class_faclev, Land_ownership_faclev)), 
         comparison_2 = factor(comparison_2, levels = c(Landuse_class_faclev, Land_ownership_faclev)))
# plot the pairwise test results
pairwise_plot_list <- list()
for (i in c("Landuse_class", "Land_ownership")) {
  pairwise_plot_list <- c(pairwise_plot_list,
                          list(ggplot(subset(pairwise_df, taxa == "tree" & attr == i), 
                                      aes(comparison_1, comparison_2, fill = p))+
                                 geom_tile() + geom_text(aes(label = round(p*100)), size = 2.5) +
                                 scale_fill_gradient2(high = "blue", low = "red", 
                                                      midpoint = 0.05, limits = c(0, 0.05)) + 
                                 theme(axis.text.x = element_text(angle = 90)) + 
                                 xlab(NULL) + ylab(NULL) + guides(fill = FALSE) + 
                                 facet_grid(index ~ attr, scales = "free", switch = "both")))
}
for (i in c("Landuse_class", "Land_ownership")) {
  pairwise_plot_list <- c(pairwise_plot_list,
                          list(ggplot(subset(pairwise_df, taxa == "shrub" & attr == i), 
                                      aes(comparison_1, comparison_2, fill = p))+
                                 geom_tile() + geom_text(aes(label = round(p*100)), size = 2.5) +
                                 scale_fill_gradient2(high = "blue", low = "red", 
                                                      midpoint = 0.05, limits = c(0, 0.05)) + 
                                 theme(axis.text.x = element_text(angle = 90)) + 
                                 xlab(NULL) + ylab(NULL) + guides(fill = FALSE) + 
                                 facet_grid(index ~ attr, scales = "free", switch = "both")))
}
Rmisc::multiplot(plotlist = pairwise_plot_list, cols = 4)

## models for diversity indexes ~ distance
# point plot: indexes ~ distance colored by land use class
index_dist_plot_list <- vector("list", 2)
index_dist_plot_list[[1]] <- na.omit(tree_diversity_long) %>% 
  subset(attr_value %in% Landuse_class_faclev) %>%
  ggplot(aes(Dist, index_value)) + geom_point(aes(color = attr_value)) + facet_wrap( ~ index, nrow = 1, scales = "free")
index_dist_plot_list[[2]] <- na.omit(shrub_diversity_long) %>% 
  subset(attr_value %in% Landuse_class_faclev) %>%
  ggplot(aes(Dist, index_value)) + geom_point(aes(color = attr_value)) + facet_wrap( ~ index, nrow = 1, scales = "free")
Rmisc::multiplot(plotlist = index_dist_plot_list)
rm(index_dist_plot_list)

# box plot: indexes ~ distance level
index_distlev_plot_list <- vector("list", 2)
index_distlev_plot_list[[1]] <- tree_diversity %>% 
  select("Sum_stem", "Richness", "Shannon", "Evenness", "Dist_level") %>% 
  pivot_longer(cols = c("Sum_stem", "Richness", "Shannon", "Evenness"), 
               names_to = "index", values_to = "index_value") %>% 
  mutate(index = factor(index, levels = c("Sum_stem", "Richness", "Shannon", "Evenness"))) %>%
  na.omit() %>% 
  ggplot(aes(Dist_level, index_value)) + geom_boxplot() + facet_wrap(~ index, nrow = 1, scales = "free")
index_distlev_plot_list[[2]] <- shrub_diversity %>% 
  select("Sum_area", "Richness", "Shannon", "Evenness", "Dist_level") %>% 
  pivot_longer(cols = c("Sum_area", "Richness", "Shannon", "Evenness"), 
               names_to = "index", values_to = "index_value") %>% 
  mutate(index = factor(index, levels = c("Sum_area", "Richness", "Shannon", "Evenness"))) %>%
  na.omit() %>%
  ggplot(aes(Dist_level, index_value)) + geom_boxplot() + facet_wrap(~ index, nrow = 1, scales = "free")
Rmisc::multiplot(plotlist = index_distlev_plot_list)
rm(index_distlev_plot_list)

# regsubsets() for trees
par(mfrow = c(2,2))
for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
  plot(regsubsets(tree_diversity[[i]] ~ Dist + I(Dist^2) + I(Dist^3) + I(Dist^4) + I(Dist^5), 
                  data = tree_diversity), scale = "adjr2", main = i)
}
# statistic analysis
summary(lm(Sum_stem ~ I(Dist^2) + I(Dist^3) + I(Dist^4) + I(Dist^5), data = tree_diversity))
summary(lm(Richness ~ Dist + I(Dist^5), data = tree_diversity))
summary(lm(Shannon ~ I(Dist^2) +  I(Dist^5), data = tree_diversity))
summary(lm(Evenness ~ Dist + I(Dist^3) +  I(Dist^5), data = tree_diversity))
# only Sum_stem ~ Dist etc. shows weakly significant p value

# point & line plot: Sum_stem ~ Dist
lmcurve <- data.frame(x = seq(from = 0, to = 13000, by = 1000))
lmcurve$y <- 3.648e+00 + 2.522e-06*lmcurve$x^2 - 8.641e-10*lmcurve$x^3 + 1.000e-13*lmcurve$x^4 - 3.761e-18*lmcurve$x^5
ggplot(tree_diversity, aes(Dist, Sum_stem)) + geom_point() + 
  geom_smooth(aes(x, y), method = "loess", data = lmcurve, se = FALSE)
rm(lmcurve)

# regsubsets() for shrub
par(mfrow = c(2,2))
for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
  plot(regsubsets(shrub_diversity[[i]] ~ Dist + I(Dist^2) + I(Dist^3) + I(Dist^4) + I(Dist^5), 
                  data = shrub_diversity), scale = "adjr2", main = i)
}
#
# statistic analysis
summary(lm(Sum_area ~ I(Dist^2), data = shrub_diversity))
summary(lm(Evenness ~ I(Dist^2) + I(Dist^3) + I(Dist^4) + I(Dist^5), data = shrub_diversity))

# the distribute of the land use types along the distance
ggplot(tree_diversity, aes(Dist)) + geom_histogram(binwidth = 1000) + facet_wrap(~ Landuse_class)
ggplot(tree_diversity, aes(Dist)) + geom_histogram(binwidth = 1000) + facet_wrap(~ Landuse_agg)
ggplot(tree_diversity, aes(Dist_level)) + geom_bar() + facet_wrap(~ Landuse_agg)



## models for plot plant attr ~ distance & distance level
# plot plant attr ~ distance colored by land use class
index_dist_plot_list <- vector("list", 2)
index_dist_plot_list[[1]] <- tree_diversity %>% 
  select("perc_planted", "perc_private", "perc_nonstreet", "perc_native", "Dist", "Landuse_class") %>% 
  pivot_longer(cols = c("perc_planted", "perc_private", "perc_nonstreet", "perc_native"), 
               names_to = "plant_attr", values_to = "plant_attr_value") %>% 
  ggplot(aes(Dist, plant_attr_value)) + geom_point(aes(color = Landuse_class)) + facet_wrap(~ plant_attr, nrow = 1, scales = "free") 
index_dist_plot_list[[2]] <- shrub_diversity %>% 
  select("perc_planted", "perc_private", "perc_nonstreet", "perc_native", "Dist", "Landuse_class") %>% 
  pivot_longer(cols = c("perc_planted", "perc_private", "perc_nonstreet", "perc_native"), 
               names_to = "plant_attr", values_to = "plant_attr_value") %>% 
  ggplot(aes(Dist, plant_attr_value)) + geom_point(aes(color = Landuse_class)) +
  facet_wrap(~ plant_attr, nrow = 1, scales = "free")
Rmisc::multiplot(plotlist = index_dist_plot_list)
rm(index_dist_plot_list)

# box plot: plot plant attr ~ distance level
index_distlev_plot_list <- vector("list", 2)
index_distlev_plot_list[[1]] <- tree_diversity %>% 
  select("perc_planted", "perc_private", "perc_nonstreet", "perc_native", "Dist_level") %>% 
  pivot_longer(cols = c("perc_planted", "perc_private", "perc_nonstreet", "perc_native"), 
               names_to = "plant_attr", values_to = "plant_attr_value") %>% 
  mutate(plant_attr = factor(plant_attr, levels = c("perc_planted", "perc_private", "perc_nonstreet", "perc_native"))) %>%
  na.omit() %>% 
  ggplot(aes(Dist_level, plant_attr_value)) + geom_boxplot() + facet_wrap(~ plant_attr, nrow = 1, scales = "free")
index_distlev_plot_list[[2]] <- shrub_diversity %>% 
  select("perc_planted", "perc_private", "perc_nonstreet", "perc_native", "Dist_level") %>% 
  pivot_longer(cols = c("perc_planted", "perc_private", "perc_nonstreet", "perc_native"), 
               names_to = "plant_attr", values_to = "plant_attr_value") %>% 
  mutate(plant_attr = factor(plant_attr, levels = c("perc_planted", "perc_private", "perc_nonstreet", "perc_native"))) %>%
  na.omit() %>%
  ggplot(aes(Dist_level, plant_attr_value)) + geom_boxplot() + facet_wrap(~ plant_attr, nrow = 1, scales = "free")
Rmisc::multiplot(plotlist = index_distlev_plot_list)
rm(index_distlev_plot_list)

# regsubsets() for trees
par(mfrow = c(2,2))
for (i in c("perc_planted", "perc_private", "perc_native")) {
  plot(regsubsets(tree_diversity[[i]] ~ Dist + I(Dist^2) + I(Dist^3) + I(Dist^4) + I(Dist^5), 
                  data = tree_diversity), scale = "adjr2", main = i)
}
#
# statistic analysis
summary(lm(perc_planted ~ I(Dist^2), data = tree_diversity))
summary(lm(perc_native ~ Dist + I(Dist^2) + I(Dist^5), data = tree_diversity))

# regsubsets() for shrubs
par(mfrow = c(2,2))
for (i in c("perc_planted", "perc_private", "perc_native")) {
  plot(regsubsets(shrub_diversity[[i]] ~ Dist + I(Dist^2) + I(Dist^3) + I(Dist^4) + I(Dist^5), 
                  data = shrub_diversity), scale = "adjr2", main = i)
}
#
# statistic analysis
summary(lm(perc_planted ~ I(Dist^2) + I(Dist^5), data = shrub_diversity))
summary(lm(perc_private ~ I(Dist^3) + I(Dist^4) + I(Dist^5), data = shrub_diversity))
summary(lm(perc_native ~ I(Dist^5), data = shrub_diversity))
#
# point & line plot: perc_planted ~ Dist of shrubs
lmcurve <- data.frame(x = seq(from = 0, to = 13000, by = 1000))
lmcurve$y <- 9.775e-01 - 3.497e-09*lmcurve$x^2 + 2.040e-21*lmcurve$x^5
ggplot(tree_diversity, aes(Dist, perc_planted)) + geom_point() + 
  geom_smooth(aes(x, y), method = "loess", data = lmcurve, se = FALSE)
rm(lmcurve)
#
# point & line plot: perc_private ~ Dist of shrubs
lmcurve <- data.frame(x = seq(from = 0, to = 13000, by = 1000))
lmcurve$y <- 6.388e-01 + 6.419e-12*lmcurve$x^3 - 1.219e-15*lmcurve$x^4 + 5.672e-20*lmcurve$x^5
ggplot(tree_diversity, aes(Dist, perc_private)) + geom_point() + 
  geom_smooth(aes(x, y), method = "loess", data = lmcurve, se = FALSE)
rm(lmcurve)
#
# point & line plot: perc_native ~ Dist of shrubs
lmcurve <- data.frame(x = seq(from = 0, to = 13000, by = 1000))
lmcurve$y <- 5.161e-01 - 1.521e-21*lmcurve$x^5
ggplot(tree_diversity, aes(Dist, perc_private)) + geom_point() + 
  geom_smooth(aes(x, y), method = "loess", data = lmcurve, se = FALSE)
rm(lmcurve)



## model for indexes of native species ~ dist for trees
# regsubsets() for trees
par(mfrow = c(2,2))
for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
  plot(regsubsets(tree_native_diversity[[i]] ~ Dist + I(Dist^2) + I(Dist^3) + I(Dist^4) + I(Dist^5), 
                  data = tree_native_diversity), scale = "adjr2", main = i)
}
#
# statistic analysis
summary(lm(Sum_stem ~ I(Dist^2) + I(Dist^3) + I(Dist^4), data = tree_native_diversity))
summary(lm(Richness ~ Dist + I(Dist^3), data = tree_native_diversity))
summary(lm(Shannon ~ I(Dist^2) + I(Dist^4), data = tree_native_diversity))
summary(lm(Evenness ~ I(Dist^2), data = tree_native_diversity))



## model for indexes of native species ~ dist for shrubs
# regsubsets() for shrubs
par(mfrow = c(2,2))
for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
  plot(regsubsets(shrub_native_diversity[[i]] ~ Dist + I(Dist^2) + I(Dist^3) + I(Dist^4) + I(Dist^5), 
                  data = shrub_native_diversity), scale = "adjr2", main = i)
}
#
# statistic analysis
summary(lm(Richness ~ I(Dist^5), data = shrub_native_diversity))






