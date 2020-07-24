setwd("C:/Rdata/KUP")
library(tidyverse)
library(Rmisc)
library(vegan)
library(PerformanceAnalytics)
library(digest)
library(car)
library(dunn.test)
library(BiodiversityR)
opar <- par(no.readonly = TRUE)

# define the factor levels
Ward_faclev <- c("Ukyo-ku", "Sakyo-ku", "Kita-ku", "Kamigyo-ku", "Nakagyo-ku", "Shimogyo-ku", "Higashiyama-ku", "Yamashina-ku", "Fushimi-ku", "Minami-ku", "Nishikyo-ku")
Landuse_class_faclev <- c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")

## get data
# data of all_plot_info
all_plot_info <- read.csv("In_plot_info.csv", stringsAsFactors = FALSE) %>% 
  mutate(Ward = factor(Ward, levels = Ward_faclev), 
         Landuse_class = factor(Landuse_class, levels = Landuse_class_faclev))

# data of all_plant_info
all_plant_info <- read.csv("In_plant_info.csv", stringsAsFactors = FALSE)
all_plant_data <- read.csv("In_plant_data.csv", stringsAsFactors = FALSE) %>%
  left_join(all_plant_info, by = "Species_CN") %>%
  left_join(all_plot_info,by = "Plot_ID")
all_plant_data$Landuse_subclass <- NULL
all_plant_data$Landuse_agg <- NULL

# data of trees
tree_data <- all_plant_data %>% 
  subset(Tree_shrub == "tree")

# data of shrubs
shrub_data <- all_plant_data %>% 
  subset(Tree_shrub == "shrub")

# data of native trees
tree_native_data <- subset(tree_data, Nt_ex == "nt")

# shrub_native_data
shrub_native_data <- subset(shrub_data, Nt_ex == "nt")

# data of tree_diversity
tree_diversity <- subset(tree_data, select = c("Plot_ID", "Species_CN", "Stem")) %>%
  pivot_wider(names_from = Species_CN, values_from = Stem, 
              values_fn = list(Stem = sum), values_fill = list(Stem = 0)) %>% 
  mutate(Sum_stem = rowSums(.[2:ncol(.)]), 
         Richness = specnumber(.[2:ncol(.)]),
         Shannon = diversity(.[2:ncol(.)], index = "shannon"), 
         Simpson = diversity(.[2:ncol(.)], index = "simpson"),
         Evenness = Shannon / log(Richness)) %>%
  left_join(all_plot_info,by = "Plot_ID") %>% 
  as.data.frame()
tree_diversity_perc_planted <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pla_spo == "P", Stem, 0)/sum(Stem))) 
tree_diversity_perc_nonpot <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pot == "N", Stem, 0)/sum(Stem)))
tree_diversity_perc_private <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pub_pri == "N", Stem, 0)/sum(Stem)))
tree_diversity_perc_nonstreet <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Street == "N", Stem, 0)/sum(Stem)))
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
         Richness = specnumber(.[2:ncol(.)]),
         Shannon = diversity(.[2:ncol(.)], index = "shannon"), 
         Simpson = diversity(.[2:ncol(.)], index = "simpson"),
         Evenness = Shannon / log(Richness)) %>%
  left_join(all_plot_info,by = "Plot_ID") %>% 
  as.data.frame()
tree_native_diversity_perc_planted <- tree_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pla_spo == "P", Stem, 0)/sum(Stem))) 
tree_native_diversity_perc_nonpot <- tree_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pot == "N", Stem, 0)/sum(Stem)))
tree_native_diversity_perc_private <- tree_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pub_pri == "N", Stem, 0)/sum(Stem)))
tree_native_diversity_perc_nonstreet <- tree_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Street == "N", Stem, 0)/sum(Stem)))
tree_native_diversity <- tree_native_diversity %>% mutate(
  perc_planted = tree_native_diversity_perc_planted$perc, 
  perc_nonpot = tree_native_diversity_perc_nonpot$perc, 
  perc_private = tree_native_diversity_perc_private$perc, 
  perc_nonstreet = tree_native_diversity_perc_nonstreet$perc
)
rm(tree_native_diversity_perc_planted, tree_native_diversity_perc_nonpot, tree_native_diversity_perc_private, tree_native_diversity_perc_nonstreet)

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
  dplyr::summarise(perc = sum(ifelse(Pla_spo == "P", Area, 0)/sum(Area))) 
shrub_diversity_perc_nonpot <- shrub_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pot == "N", Area, 0)/sum(Area)))
shrub_diversity_perc_private <- shrub_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pub_pri == "N", Area, 0)/sum(Area)))
shrub_diversity_perc_nonstreet <- shrub_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Street == "N", Area, 0)/sum(Area)))
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
  dplyr::summarise(perc = sum(ifelse(Pla_spo == "P", Area, 0)/sum(Area))) 
shrub_native_diversity_perc_nonpot <- shrub_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pot == "N", Area, 0)/sum(Area)))
shrub_native_diversity_perc_private <- shrub_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pub_pri == "N", Area, 0)/sum(Area)))
shrub_native_diversity_perc_nonstreet <- shrub_native_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Street == "N", Area, 0)/sum(Area)))
shrub_native_diversity <- shrub_native_diversity %>% mutate(
  perc_planted = shrub_native_diversity_perc_planted$perc, 
  perc_nonpot = shrub_native_diversity_perc_nonpot$perc, 
  perc_private = shrub_native_diversity_perc_private$perc, 
  perc_nonstreet = shrub_native_diversity_perc_nonstreet$perc
)
rm(shrub_native_diversity_perc_planted, shrub_native_diversity_perc_nonpot, shrub_native_diversity_perc_private, shrub_native_diversity_perc_nonstreet)

# tree diversity longer and shrub diversity longer dataset
tree_diversity_long <- 
  subset(tree_diversity, select = c("Sum_stem", "Richness", "Shannon", "Evenness", "Landuse_class")) %>% 
  pivot_longer(cols = c("Sum_stem", "Richness", "Shannon", "Evenness"), 
               names_to = "index", values_to = "index_value") %>% 
  pivot_longer(cols = c("Landuse_class"), 
               names_to = "attr", values_to = "attr_value") %>% 
  mutate(index = factor(index, levels = c("Sum_stem", "Richness", "Shannon", "Evenness")), 
         attr = factor(attr, levels = c("Landuse_class")), 
         attr_value = factor(attr_value, levels = c(Landuse_class_faclev)))

shrub_diversity_long <- 
  subset(shrub_diversity, select = c("Sum_area", "Richness", "Shannon", "Evenness", "Landuse_class")) %>% 
  pivot_longer(cols = c("Sum_area", "Richness", "Shannon", "Evenness"), 
               names_to = "index", values_to = "index_value") %>% 
  pivot_longer(cols = c("Landuse_class"), 
               names_to = "attr", values_to = "attr_value") %>% 
  mutate(index = factor(index, levels = c("Sum_area", "Richness", "Shannon","Evenness")), 
         attr = factor(attr, levels = c("Landuse_class")), 
         attr_value = factor(attr_value, levels = c(Landuse_class_faclev)))

# tree native diversity longer and shrub native diversity longer dataset
tree_native_diversity_long <- 
  subset(tree_native_diversity, select = c("Sum_stem", "Richness", "Shannon", "Evenness", "Landuse_class")) %>% 
  pivot_longer(cols = c("Sum_stem", "Richness", "Shannon", "Evenness"), 
               names_to = "index", values_to = "index_value") %>% 
  pivot_longer(cols = c("Landuse_class"), 
               names_to = "attr", values_to = "attr_value") %>% 
  mutate(index = factor(index, levels = c("Sum_stem", "Richness", "Shannon", "Evenness")), 
         attr = factor(attr, levels = c("Landuse_class")), 
         attr_value = factor(attr_value, levels = c(Landuse_class_faclev)))

shrub_native_diversity_long <- 
  subset(shrub_native_diversity, select = c("Sum_area", "Richness", "Shannon", "Evenness", "Landuse_class")) %>% 
  pivot_longer(cols = c("Sum_area", "Richness", "Shannon", "Evenness"), 
               names_to = "index", values_to = "index_value") %>% 
  pivot_longer(cols = c("Landuse_class"), 
               names_to = "attr", values_to = "attr_value") %>% 
  mutate(index = factor(index, levels = c("Sum_area", "Richness", "Shannon","Evenness")), 
         attr = factor(attr, levels = c("Landuse_class")), 
         attr_value = factor(attr_value, levels = c(Landuse_class_faclev)))

number_tree_species <- length(unique(tree_data$Species_CN))
number_shrub_species <- length(unique(shrub_data$Species_CN))
number_tree_native_species <- length(unique(tree_native_data$Species_CN))
number_shrub_native_species <- length(unique(shrub_native_data$Species_CN))



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

# top species families of all plants by species number
all_plant_info %>% group_by(Family) %>% 
  dplyr::summarise(number = n(), prop = number/nrow(all_plant_info)) %>% 
  arrange(desc(prop))

# the number of trees or area of shrubs
cat("number of trees:", nrow(tree_data), "in", nrow(tree_diversity), "plot", "\n", 
    "area of shrubs:", sum(shrub_data$Area), "m2 in", nrow(shrub_diversity), "plot")

# top species families of trees by individual
tree_data %>% 
  group_by(Family) %>% dplyr::summarise(Number = sum(Stem), Prop = Number/sum(tree_data$Stem)) %>% 
  arrange(desc(Prop))

# top species families of shrubs by area
shrub_data %>% 
  group_by(Family) %>% dplyr::summarise(SArea = sum(Area), Prop = SArea/sum(shrub_data$Area)) %>% 
  arrange(desc(Prop))



## attributes of the tree and shrub
# the exotic vs. native by species
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
  barplot(tapply(subset(all_plant_data, Tree_shrub == "tree")[, "Stem"], 
                 subset(all_plant_data, Tree_shrub == "tree")[, i], 
                 sum), ylim = c(0, 1400))
  title(main = paste("(", letters[j], ")"), adj = 0)
}
for (i in c("Pla_spo", "Pub_pri", "Nt_ex")) {
  j <- j + 1
  barplot(tapply(subset(all_plant_data, Tree_shrub == "shrub")[, "Area"], 
                 subset(all_plant_data, Tree_shrub == "shrub")[, i], 
                 sum), ylim = c(0, 1200))
  title(main = paste("(", letters[j], ")"), adj = 0)
}
par(opar)



## Species accumulation curve 
# Species accumulation curve for trees
tree_accum_list <- vector("list", 3)
for (i in Landuse_class_faclev) {
  tree_accum_list[[1]] <- c(tree_accum_list[[1]], specaccum(subset(tree_diversity, Landuse_class == i)[,2:(number_tree_species+1)])$sites)
  tree_accum_list[[2]] <- c(tree_accum_list[[2]], specaccum(subset(tree_diversity, Landuse_class == i)[,2:(number_tree_species+1)])$richness)
  tree_accum_list[[3]] <- c(tree_accum_list[[3]], rep(i, length(specaccum(subset(tree_diversity, Landuse_class == i)[,2:(number_tree_species+1)])$sites)))
}
tree_accum_df <- data.frame(
  number_of_plot = tree_accum_list[[1]], 
  number_of_species = tree_accum_list[[2]], 
  landuse_class = tree_accum_list[[3]]
) %>% mutate(landuse_class = factor(landuse_class, levels = Landuse_class_faclev))
tree_accum_plot <- ggplot(tree_accum_df, aes(x = number_of_plot, y = number_of_species, color = landuse_class)) + geom_line(size = 1)
# Species accumulation curve for shrubs
shrub_accum_list <- vector("list", 3)
for (i in Landuse_class_faclev) {
  shrub_accum_list[[1]] <- c(shrub_accum_list[[1]], specaccum(subset(shrub_diversity, Landuse_class == i)[,2:(number_shrub_species+1)])$sites)
  shrub_accum_list[[2]] <- c(shrub_accum_list[[2]], specaccum(subset(shrub_diversity, Landuse_class == i)[,2:(number_shrub_species+1)])$richness)
  shrub_accum_list[[3]] <- c(shrub_accum_list[[3]], rep(i, length(specaccum(subset(shrub_diversity, Landuse_class == i)[,2:(number_shrub_species+1)])$sites)))
}
shrub_accum_df <- data.frame(
  number_of_plot = shrub_accum_list[[1]], 
  number_of_species = shrub_accum_list[[2]], 
  landuse_class = shrub_accum_list[[3]]
) %>% mutate(landuse_class = factor(landuse_class, levels = Landuse_class_faclev))
shrub_accum_plot <- ggplot(shrub_accum_df, aes(x = number_of_plot, y = number_of_species, color = landuse_class)) + geom_line(size = 1)
# plots of species accumulation curves for trees and shrubs 
ggarrange(tree_accum_plot, shrub_accum_plot, common.legend = T, legend = "right")



## rank ahundance plot
# the rank abundance plot by land use types of trees
tree_rankabun_list <- vector("list", 5)
for (i in c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")) {
  tree_rankabun_ori <- as.data.frame(rankabundance(subset(tree_diversity, Landuse_class == i)[2:length(unique(tree_data$Species_CN))]))
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
tree_rankabun_df <- tree_rankabun_df[which(tree_rankabun_df$abundance != 0),]

# the rank abundance plot by land use types of shrubs
shrub_rankabun_list <- vector("list", 5)
for (i in c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")) {
  shrub_rankabun_ori <- as.data.frame(rankabundance(subset(shrub_diversity, Landuse_class == i)[2:length(unique(shrub_data$Species_CN))]))
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
shrub_rankabun_df <- shrub_rankabun_df[which(shrub_rankabun_df$abundance != 0),]

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
tree_mds_selected <- tree_diversity %>% filter(Sum_stem > 3)
tree_mds_metaMDS <- tree_mds_selected %>% 
  select(2:(grep("Sum_stem", colnames(tree_mds_selected))-1)) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
tree_mds_metaMDS$stress
stressplot(tree_mds_metaMDS)
tree_mds_selected <- cbind(tree_mds_selected, tree_mds_metaMDS$points)
#
# nmds calculation for shrub
shrub_mds_selected <- shrub_diversity %>% filter(Sum_area > 5)
shrub_mds_metaMDS <- shrub_mds_selected %>% 
  select(2:(grep("Sum_area", colnames(shrub_mds_selected))-1)) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
shrub_mds_metaMDS$stress
stressplot(shrub_mds_metaMDS)
shrub_mds_selected <- cbind(shrub_mds_selected, shrub_mds_metaMDS$points)
#
# mds plot for trees and shrubs
# general Anosim of trees and shrubs
tree_landuseclass_anosim_result <- anosim(tree_mds_selected[2:(number_tree_species+1)], tree_mds_selected$Landuse_class)
shrub_landuseclass_anosim_result <- anosim(shrub_mds_selected[2:(number_shrub_species+1)], shrub_mds_selected$Landuse_class)
# get statistic results as labels for the mds plots
tree_landuseclass_mds_lab <- paste("stress=", round(tree_mds_metaMDS$stress, digits = 3), 
                                   ", R=", round(tree_landuseclass_anosim_result$statistic, digits = 3), 
                                   ", p=", round(tree_landuseclass_anosim_result$signif, digits = 3), 
                                   sep = "")
shrub_landuseclass_mds_lab <- paste("stress=", round(shrub_mds_metaMDS$stress, digits = 3), 
                                    ", R=", round(shrub_landuseclass_anosim_result$statistic, digits = 3), 
                                    ", p=", round(shrub_landuseclass_anosim_result$signif, digits = 3), 
                                    sep = "")
# get hull for mds plots of trees and shrubs
tree_find_hull <- function(tree_mds_selected) tree_mds_selected[chull(tree_mds_selected$MDS1, tree_mds_selected$MDS2), ]
tree_hulls <- ddply(tree_mds_selected, "Landuse_class", tree_find_hull)
shrub_find_hull <- function(shrub_mds_selected) shrub_mds_selected[chull(shrub_mds_selected$MDS1, shrub_mds_selected$MDS2), ]
shrub_hulls <- ddply(shrub_mds_selected, "Landuse_class", shrub_find_hull)
# mds plots for trees and shrubs by land use types and land ownership 
Rmisc::multiplot(plotlist = list(
  ggplot(tree_mds_selected, aes(MDS1, MDS2, color = Landuse_class)) + geom_point(size=3) +
    labs(title = "Tree - Land use type", subtitle = tree_landuseclass_mds_lab) +
    geom_polygon(data=tree_hulls, alpha = 0, aes(fill=Landuse_class), size=1) +theme_bw(), 
  ggplot(shrub_mds_selected, aes(MDS1, MDS2, color = Landuse_class)) + geom_point(size=3) +
    labs(title = "Shrub - Land use type", subtitle = shrub_landuseclass_mds_lab) +
    geom_polygon(data=shrub_hulls, alpha = 0, aes(fill=Landuse_class), size=1) +theme_bw()
), layout = matrix(1:2, nrow = 1, byrow = T))
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
                                  anosim(tree_mds_selected_sub[2:(number_tree_species+1)], 
                                         tree_mds_selected_sub$Landuse_class)$signif)
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
                                   anosim(shrub_mds_selected_sub[2:(number_shrub_species+1)], 
                                          shrub_mds_selected_sub$Landuse_class)$signif)
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
  pvalue_list <- vector("list", 2)
  for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
    pvalue_list[[1]] <- c(pvalue_list[[1]], i)
    pvalue_list[[2]] <- c(pvalue_list[[2]], 
                          round(kruskal.test(tree_diversity[, i] ~ tree_diversity$Landuse_class)$p.value,digits = 3))
  }
  pvalue <- data.frame(index = pvalue_list[[1]],
                       pvalue = pvalue_list[[2]])
  pvalue$label <- NA
  pvalue$label[pvalue$pvalue >= 0.05] <- 
    paste("p=", sprintf("%.3f",pvalue$pvalue[pvalue$pvalue>0.05]), sep = "")
  pvalue$label[pvalue$pvalue < 0.05 & pvalue$pvalue >= 0.01] <- 
    paste("p=", sprintf("%.3f", pvalue$pvalue[pvalue$pvalue < 0.05 & pvalue$pvalue >= 0.01]), "*", sep = "")
  pvalue$label[pvalue$pvalue < 0.01 & pvalue$pvalue >= 0.001] <- 
    paste("p=", sprintf("%.3f", pvalue$pvalue[pvalue$pvalue < 0.01 & pvalue$pvalue >= 0.001]), "**", sep = "")
  pvalue$label[pvalue$pvalue < 0.001] <- 
    paste("p=", sprintf("%.3f", pvalue$pvalue[pvalue$pvalue < 0.001]), "***", sep = "")
}
boxplot_list_index_attr[[1]] <- tree_diversity_long %>% 
  na.omit() %>%
  ggplot(aes(attr_value, index_value)) + geom_boxplot() + 
  facet_grid(index ~ attr, scales = "free", space = "free_x", switch = "both") + 
  scale_y_continuous(expand = expand_scale(mult = c(0.05,0.3))) +
  geom_text(data = pvalue, aes(x =Inf, y = Inf, label = label), size=3.5, hjust = 1.05, vjust = 1.5) +
  theme(axis.text = element_text(angle = 90)) + 
  labs(title = "(a)", x = NULL, y = NULL)
#
# for shrub
{
  pvalue_list <- vector("list", 2)
  for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
    pvalue_list[[1]] <- c(pvalue_list[[1]], i)
    pvalue_list[[2]] <- c(pvalue_list[[2]], 
                          round(kruskal.test(shrub_diversity[, i] ~ shrub_diversity$Landuse_class)$p.value,digits = 3))
  }
  pvalue <- data.frame(index = pvalue_list[[1]],
                       pvalue = pvalue_list[[2]])
  pvalue$label <- NA
  pvalue$label[pvalue$pvalue >= 0.05] <- 
    paste("p=", sprintf("%.3f",pvalue$pvalue[pvalue$pvalue>0.05]), sep = "")
  pvalue$label[pvalue$pvalue < 0.05 & pvalue$pvalue >= 0.01] <- 
    paste("p=", sprintf("%.3f", pvalue$pvalue[pvalue$pvalue < 0.05 & pvalue$pvalue >= 0.01]), "*", sep = "")
  pvalue$label[pvalue$pvalue < 0.01 & pvalue$pvalue >= 0.001] <- 
    paste("p=", sprintf("%.3f", pvalue$pvalue[pvalue$pvalue < 0.01 & pvalue$pvalue >= 0.001]), "**", sep = "")
  pvalue$label[pvalue$pvalue < 0.001] <- 
    paste("p=", sprintf("%.3f", pvalue$pvalue[pvalue$pvalue < 0.001]), "***", sep = "")
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



## pairwise dunn test of diversity ~ landuse class
pairwise_list <- vector("list", 4)
# list of pairwise test of diversity of tree
for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
  pairwise_list[[1]] <- c(pairwise_list[[1]], 
                          rep("tree", 15))
  pairwise_list[[2]] <- c(pairwise_list[[2]], 
                          rep(i, 15))
  pairwise_list[[3]] <- c(pairwise_list[[3]], 
                          dunn.test(tree_diversity[, i], tree_diversity$Landuse_class)$comparisons)
  pairwise_list[[4]] <- c(pairwise_list[[4]], 
                          dunn.test(tree_diversity[, i], tree_diversity$Landuse_class)$P.adjusted)
}
# list of pairwise test of diversity and attrs of shrub
for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
  pairwise_list[[1]] <- c(pairwise_list[[1]], 
                          rep("shrub", 15))
  pairwise_list[[2]] <- c(pairwise_list[[2]], 
                          rep(i, 15))
  pairwise_list[[3]] <- c(pairwise_list[[3]], 
                          dunn.test(shrub_diversity[, i], shrub_diversity$Landuse_class)$comparisons)
  pairwise_list[[4]] <- c(pairwise_list[[4]], 
                          dunn.test(shrub_diversity[, i], shrub_diversity$Landuse_class)$P.adjusted)
}
# data frame of pairwise test of tree and shrub
pairwise_df_up <- data.frame(taxa = pairwise_list[[1]], 
                             index = pairwise_list[[2]], 
                             comparison = pairwise_list[[3]], 
                             p = pairwise_list[[4]]) %>% 
  separate(comparison, into = c("comparison_1", "comparison_2"), sep = " - ")
pairwise_df_down <- data.frame(
  taxa = pairwise_df_up$taxa, 
  index = pairwise_df_up$index, 
  comparison_1 = pairwise_df_up$comparison_2, 
  comparison_2 = pairwise_df_up$comparison_1, 
  p = pairwise_df_up$p)
pairwise_df <- rbind(pairwise_df_up, pairwise_df_down) %>% 
  mutate(index = factor(index, levels = c("Sum_stem", "Sum_area", "Richness", "Shannon", "Evenness")), 
         comparison_1 = factor(comparison_1, levels = Landuse_class_faclev), 
         comparison_2 = factor(comparison_2, levels = Landuse_class_faclev))
# plot the pairwise test results
ggplot(subset(pairwise_df, taxa == "tree"), 
       aes(comparison_1, comparison_2, fill = p))+
  geom_tile() + geom_text(aes(label = round(p*100)), size = 2.5) +
  scale_fill_gradient2(high = "blue", low = "red", 
                       midpoint = 0.05, limits = c(0, 0.05)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  xlab(NULL) + ylab(NULL) + guides(fill = FALSE) + 
  facet_wrap(~ index, scales = "free")
ggplot(subset(pairwise_df, taxa == "shrub"), 
       aes(comparison_1, comparison_2, fill = p))+
  geom_tile() + geom_text(aes(label = round(p*100)), size = 2.5) +
  scale_fill_gradient2(high = "blue", low = "red", 
                       midpoint = 0.05, limits = c(0, 0.05)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  xlab(NULL) + ylab(NULL) + guides(fill = FALSE) + 
  facet_wrap(~ index, scales = "free")


