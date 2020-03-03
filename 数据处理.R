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
opar <- par(no.readonly = TRUE)

# define the factor levels
Ward_faclev <- c("Ukyo-ku", "Sakyo-ku", "Kita-ku", "Kamigyo-ku", "Nakagyo-ku", "Shimogyo-ku", "Higashiyama-ku", "Yamashina-ku", "Fushimi-ku", "Minami-ku", "Nishikyo-ku")
Landuse_class_faclev <- c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")
Landscaping_faclev <- c("Private area", "Propriater area", "Public area")

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
all_plot_info$Landscaping <- NA
all_plot_info$Landscaping[all_plot_info$Landuse_detail %in% c("道路", "高层商业", "公园", "河流", "机构", "政府团地")] <- "Public area"
all_plot_info$Landscaping[all_plot_info$Landuse_detail %in% c("低层私宅", "私宅")] <- "Private area"
all_plot_info$Landscaping[all_plot_info$Landuse_detail %in% c("低层公寓", "高层公寓", "低层商业", "工厂", "农田", "其他", "球场", "寺庙神社", "学校")] <- "Propriater area"

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

# data of native trees
tree_native_data <- subset(tree_data, Nt_ex == "nt")

# data of shrubs
shrub_data <- all_plant_data %>% 
  subset(Tree_shrub == "Shrub")

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

# tree diversity longer and shrub diversity longer dataset
tree_diversity_long <- 
  subset(tree_diversity, select = c("Sum_stem", "Richness", "Shannon", "Evenness", "Landuse_class", "Landscaping", "Dist")) %>% 
  pivot_longer(cols = c("Sum_stem", "Richness", "Shannon", "Evenness"), 
               names_to = "index", values_to = "index_value") %>% 
  pivot_longer(cols = c("Landuse_class", "Landscaping"), 
               names_to = "attr", values_to = "attr_value") %>% 
  mutate(index = factor(index, levels = c("Sum_stem", "Richness", "Shannon", "Evenness")), 
         attr = factor(attr, levels = c("Landuse_class", "Landscaping")), 
         attr_value = factor(attr_value, levels = c(Landuse_class_faclev, Landscaping_faclev)))

shrub_diversity_long <- 
  subset(shrub_diversity, select = c("Sum_area", "Richness", "Shannon", "Evenness", "Landuse_class", "Landscaping", "Dist")) %>% 
  pivot_longer(cols = c("Sum_area", "Richness", "Shannon", "Evenness"), 
               names_to = "index", values_to = "index_value") %>% 
  pivot_longer(cols = c("Landuse_class", "Landscaping"), 
               names_to = "attr", values_to = "attr_value") %>% 
  mutate(index = factor(index, levels = c("Sum_area", "Richness", "Shannon","Evenness")), 
         attr = factor(attr, levels = c("Landuse_class", "Landscaping")), 
         attr_value = factor(attr_value, levels = c(Landuse_class_faclev, Landscaping_faclev)))

### analysis begin
## the total species, genera and families of all plants
length(unique(all_plant_data$Species_CN))
length(unique(all_plant_data$Genus))
length(unique(all_plant_data$Family))
# of trees
length(unique(tree_data$Species_CN))
# and of shrubs
length(unique(shrub_data$Species_CN))
# common species for trees and for shrubs
length(intersect(unique(tree_data$Species_CN), 
                 unique(shrub_data$Species_CN)))
# species solely for trees
length(setdiff(unique(tree_data$Species_CN), 
                 unique(shrub_data$Species_CN)))
# Species solely for shrubs
length(setdiff(unique(shrub_data$Species_CN), 
                 unique(tree_data$Species_CN)))

## top species by species and number or area
# by species
# of all plants
all_plant_info %>% group_by(Family) %>%
  dplyr::summarise(Species =n(), Prop = n()/nrow(all_plant_info)) %>% 
  arrange(desc(Prop)) %>% 
  print() %>% ggplot(aes(reorder(Family, -Prop), Prop)) + geom_col()
# of trees by species and number/area
# tree~species
tree_plant_info %>% 
  group_by(Family) %>% dplyr::summarise(Species = n(), Prop = n()/nrow(tree_plant_info)) %>% 
  arrange(desc(Prop)) 
# tree~number
tree_data %>% 
  group_by(Family) %>% dplyr::summarise(Number = sum(Stem), Prop = Number/sum(tree_data$Stem)) %>% 
  arrange(desc(Prop))
# shrub~species
shrub_plant_info %>% 
  group_by(Family) %>% dplyr::summarise(Species = n(), Prop = n()/nrow(shrub_plant_info)) %>% 
  arrange(desc(Prop))
# shrub~area
shrub_data %>% 
  group_by(Family) %>% dplyr::summarise(SArea = sum(Area), Prop = SArea/sum(shrub_data$Area)) %>% 
  arrange(desc(Prop))
# rank species or number plots 
multiplot(
  tree_plant_info %>% 
    group_by(Family) %>% dplyr::summarise(Species = n(), Prop = n()/nrow(tree_plant_info)) %>% 
    ggplot(aes(reorder(Family, -Prop), Prop)) + geom_col() + 
    labs(title = "tree~species") + ylim(0, 0.2), 
  shrub_plant_info %>% 
    group_by(Family) %>% dplyr::summarise(Species = n(), Prop = n()/nrow(shrub_plant_info)) %>% 
    ggplot(aes(reorder(Family, -Prop), Prop)) + geom_col()+ 
    labs(title = "shrub~species") + ylim(0, 0.2), 
  tree_data %>% 
    group_by(Family) %>% dplyr::summarise(Number = sum(Stem), Prop = Number/sum(tree_data$Stem)) %>% 
    ggplot(aes(reorder(Family, -Prop), Prop)) + geom_col()+ 
    labs(title = "tree~number")+ ylim(0, 0.2),
  shrub_data %>% 
    group_by(Family) %>% dplyr::summarise(SArea = sum(Area), Prop = SArea/sum(shrub_data$Area)) %>% 
    ggplot(aes(reorder(Family, -Prop), Prop)) + geom_col()+ 
    labs(title = "shrub~area")+ ylim(0, 0.2), # why missing value? where?
  layout = matrix(1:4, ncol = 2)
)

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
par(mfrow= c(2,5), cex.axis = 1.5)
j <- 0
for (i in c("Pla_spo", "Pot", "Pub_pri", "Street", "Nt_ex")) {
  j <- j + 1
  barplot(tapply(subset(all_plant_data, Tree_shrub == "Tree")[, "Stem"], 
                 subset(all_plant_data, Tree_shrub == "Tree")[, i], 
                 sum), ylim = c(0, 1400))
  title(main = paste("(", letters[j], ")"), adj = 0)
}
for (i in c("Pla_spo", "Pot", "Pub_pri", "Street", "Nt_ex")) {
  j <- j + 1
  barplot(tapply(subset(all_plant_data, Tree_shrub == "Shrub")[, "Area"], 
                 subset(all_plant_data, Tree_shrub == "Shrub")[, i], 
                 sum), ylim = c(0, 1400))
  title(main = paste("(", letters[j], ")"), adj = 0)
}
par(opar)

## mds analysis
set.seed(1234)
all_plot_list_mds <- list()
# nmds plot for tree
tree_mds_selected_ID <- tree_diversity$Plot_ID[!(tree_diversity$Plot_ID %in% c(214, 261, 244, 313))]
tree_mds_selected <- tree_diversity %>% filter(Plot_ID %in% tree_mds_selected_ID)
tree_mds_metaMDS <- tree_mds_selected %>% 
  select(2:143) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
tree_mds_metaMDS$stress
stressplot(tree_mds_metaMDS)
tree_mds_selected <- cbind(tree_mds_selected, tree_mds_metaMDS$points)
# nmds plot for shrub
shrub_mds_selected_ID <- shrub_diversity$Plot_ID[!(shrub_diversity$Plot_ID %in% c(269, 214, 75, 164, 244))]
shrub_mds_selected <- shrub_diversity %>% filter(Plot_ID %in% shrub_mds_selected_ID)
shrub_mds_metaMDS <- shrub_mds_selected %>% 
  select(2:196) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
shrub_mds_metaMDS$stress
stressplot(shrub_mds_metaMDS)
shrub_mds_selected <- cbind(shrub_mds_selected, shrub_mds_metaMDS$points)
# mds plot for tree and shrub
Rmisc::multiplot(plotlist = list(
  ggplot(tree_mds_selected, aes(MDS1, MDS2, color = Landuse_class)) + geom_point(alpha = 0.7) + 
    labs(title = "Tree Land use"), 
  ggplot(tree_mds_selected, aes(MDS1, MDS2, color = Landscaping)) + geom_point(alpha = 0.7) + 
    labs(title = "Tree Landscaping"), 
  ggplot(shrub_mds_selected, aes(MDS1, MDS2, color = Landuse_class)) + geom_point(alpha = 0.7) + 
    labs(title = "Shrub Land use"), 
  ggplot(shrub_mds_selected, aes(MDS1, MDS2, color = Landscaping)) + geom_point(alpha = 0.7) + 
    labs(title = "Shrub Landscaping")
), layout = matrix(1:4, ncol = 2))
# ANOSIM
# general result of ANOSIM of tree and shrub
set.seed(1234)
for (i in c("Landuse_class", "Landscaping")) {
  x <- anosim(tree_mds_selected[,c(2:143)], tree_mds_selected[,i])
  cat("tree", i, "R=", x$statistic, "p=", x$signif, "\n")
}
for (i in c("Landuse_class", "Landscaping")) {
  x <- anosim(as.data.frame(shrub_mds_selected)[,c(2:196)], shrub_mds_selected[,i])
  cat("shrub", i, "R=", x$statistic, "p=", x$signif, "\n")
  rm(x)
}
# pairwise result of ANOSIM of tree
x <- list()
for (i in c("Landuse_class", "Landscaping")) {
  set.seed(1234)
  a <- combn(levels(factor(tree_mds_selected[,i])), 2)
  z <- vector("list", 3)
  for (j in 1:ncol(a)) {
    b <- tree_mds_selected %>% 
      filter(tree_mds_selected[[i]] == a[1,j] | tree_mds_selected[[i]] == a[2,j]) 
    c <- anosim(b[2:143], b[,i])
    z[[1]] <- c(z[[1]], as.character(a[1,j]))
    z[[2]] <- c(z[[2]], as.character(a[2,j]))
    z[[3]] <- c(z[[3]], c$signif)
  }
  y <- data.frame(comp_1 = z[[1]], comp_2 = z[[2]], p = z[[3]]) %>% 
    mutate(comp_1 = factor(comp_1, levels = levels(factor(tree_mds_selected[,i]))), 
           comp_2 = factor(comp_2, levels = levels(factor(tree_mds_selected[,i]))))
  print(y)
  x <- c(x, list(ggplot(y, aes(comp_1, comp_2)) + geom_tile(aes(fill = p)) + 
                   scale_fill_gradient(high = "orange", low = "red", limit = c(0, 0.05)) + 
                   theme(axis.text.x = element_text(angle = 90))))
}
Rmisc::multiplot(plotlist = x, layout = matrix(1:4, ncol = 2))
# pairwise result of ANOSIM of shrub
for (i in c("Landuse_class", "Landscaping")) {
  set.seed(1234)
  a <- combn(levels(factor(shrub_mds_selected[,i])), 2)
  z <- vector("list", 3)
  for (j in 1:ncol(a)) {
    b <- shrub_mds_selected %>% 
      filter(shrub_mds_selected[[i]] == a[1,j] | shrub_mds_selected[[i]] == a[2,j]) 
    c <- anosim(b[2:196], b[,i])
    z[[1]] <- c(z[[1]], as.character(a[1,j]))
    z[[2]] <- c(z[[2]], as.character(a[2,j]))
    z[[3]] <- c(z[[3]], c$signif)
  }
  y <- data.frame(comp_1 = z[[1]], comp_2 = z[[2]], p = z[[3]]) %>% 
    mutate(comp_1 = factor(comp_1, levels = levels(factor(shrub_mds_selected[,i]))), 
           comp_2 = factor(comp_2, levels = levels(factor(shrub_mds_selected[,i]))))
  print(y)
  x <- c(x, list(ggplot(y, aes(comp_1, comp_2)) + geom_tile(aes(fill = p)) + 
                   scale_fill_gradient(high = "orange", low = "red", limit = c(0, 0.05)) + 
                   theme(axis.text.x = element_text(angle = 90))))
}
Rmisc::multiplot(plotlist = x, layout = matrix(1:4, ncol = 2))
rm(a, b, c, x, y, z)

## cor among the indexes
chart.Correlation(subset(tree_diversity, select = c("Sum_stem", "Richness", "Shannon", "Simpson", "Evenness")))
chart.Correlation(subset(shrub_diversity, select = c("Sum_area", "Richness", "Shannon", "Simpson", "Evenness")))


# kruskal test & boxplot for trees
# p-value matrix for tree
set.seed(1234)
boxplot_list_index_attr <- vector("list", 2)
# for tree
{
  pvalue_list <- vector("list", 3)
  for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
    for (j in c("Landuse_class", "Landscaping")) {
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
  pvalue$label[pvalue$pvalue>0.05] <- paste("p=", pvalue$pvalue[pvalue$pvalue>0.05], sep = "")
  pvalue$label[pvalue$pvalue<0.05 & pvalue$pvalue>0.01] <- 
    paste("p=", pvalue$pvalue, "*", sep = "")[pvalue$pvalue<0.05 & pvalue$pvalue>0.01]
  pvalue$label[pvalue$pvalue<0.01 & pvalue$pvalue>0.001] <- 
    paste("p=", pvalue$pvalue, "**", sep = "")[pvalue$pvalue<0.01 & pvalue$pvalue>0.001]
  pvalue$label[pvalue$pvalue<0.001] <- 
    paste("p=", pvalue$pvalue, "***", sep = "")[pvalue$pvalue<0.001]
}

boxplot_list_index_attr[[1]] <- tree_diversity_long %>% 
  na.omit() %>%
  ggplot(aes(attr_value, index_value)) + geom_boxplot() + 
  facet_grid(index ~ attr, scales = "free", space = "free_x", switch = "both") + 
  scale_y_continuous(expand = expand_scale(mult = c(0.05,0.3))) +
  geom_text(data = pvalue, aes(x =Inf, y = Inf, label = label), size=3.5, hjust = 1.05, vjust = 1.5) +
  theme(axis.text = element_text(angle = 90))
# for shrub
{
  pvalue_list <- vector("list", 3)
  for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
    for (j in c("Landuse_class", "Landscaping")) {
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
  pvalue$label[pvalue$pvalue>0.05] <- paste("p=", pvalue$pvalue[pvalue$pvalue>0.05], sep = "")
  pvalue$label[pvalue$pvalue<0.05 & pvalue$pvalue>0.01] <- 
    paste("p=", pvalue$pvalue, "*", sep = "")[pvalue$pvalue<0.05 & pvalue$pvalue>0.01]
  pvalue$label[pvalue$pvalue<0.01 & pvalue$pvalue>0.001] <- 
    paste("p=", pvalue$pvalue, "**", sep = "")[pvalue$pvalue<0.01 & pvalue$pvalue>0.001]
  pvalue$label[pvalue$pvalue<0.001] <- 
    paste("p=", pvalue$pvalue, "***", sep = "")[pvalue$pvalue<0.001]
}

boxplot_list_index_attr[[2]] <- shrub_diversity_long %>% 
  na.omit() %>%
  ggplot(aes(attr_value, index_value)) + geom_boxplot() + 
  facet_grid(index ~ attr, scales = "free", space = "free_x", switch = "both") + 
  scale_y_continuous(expand = expand_scale(mult = c(0.05,0.3))) +
  geom_text(data = pvalue, aes(x =Inf, y = Inf, label = label), size=3.5, hjust = 1.05, vjust = 1.5) +
  theme(axis.text = element_text(angle = 90))
Rmisc::multiplot(plotlist = boxplot_list_index_attr, cols = 2)
# delete the vars
rm(pvalue, pvalue_list)

## pairwise dunn test of diversity ~ I(landuse class + ownership)
pairwise_list <- vector("list", 5)
# list of pairwise test of diversity and attrs of tree
for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
  for (j in c("Landuse_class", "Landscaping")) {
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
  for (j in c("Landuse_class", "Landscaping")) {
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
         comparison_1 = factor(comparison_1, levels = c(Landuse_class_faclev, Landscaping_faclev)), 
         comparison_2 = factor(comparison_2, levels = c(Landuse_class_faclev, Landscaping_faclev)))
# plot the pairwise test results
pairwise_plot_list <- list()
for (i in c("Landuse_class", "Landscaping")) {
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
for (i in c("Landuse_class", "Landscaping")) {
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
# plot indexes ~ distance colored by land use class
index_dist_plot_list <- vector("list", 2)
index_dist_plot_list[[1]] <- na.omit(tree_diversity_long) %>% 
  subset(attr_value %in% Landuse_class_faclev) %>%
  ggplot(aes(Dist, index_value)) + geom_point(aes(color = attr_value)) + facet_wrap( ~ index, nrow = 1, scales = "free")
index_dist_plot_list[[2]] <- na.omit(shrub_diversity_long) %>% 
  subset(attr_value %in% Landuse_class_faclev) %>%
  ggplot(aes(Dist, index_value)) + geom_point(aes(color = attr_value)) + facet_wrap( ~ index, nrow = 1, scales = "free")
Rmisc::multiplot(plotlist = index_dist_plot_list)
rm(index_dist_plot_list)

# plot indexes ~ distance level
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

# regsubsets() for tree
# plot
par(mfrow = c(2,2))
for (i in c("Sum_stem", "Richness", "Shannon", "Evenness", "perc_planted", "perc_private", "perc_nonstreet", "perc_native")) {
  plot(regsubsets(tree_diversity[[i]] ~ Dist + I(Dist^2) + I(Dist^3) + I(Dist^4) + I(Dist^5), 
                  data = tree_diversity), scale = "adjr2", main = i)
}
# statistic analysis
summary(lm(Sum_stem ~ I(Dist^2) + I(Dist^3) + I(Dist^4) + I(Dist^5), data = tree_diversity))
summary(lm(Richness ~ Dist + I(Dist^5), data = tree_diversity))
summary(lm(Shannon ~ I(Dist^2) +  I(Dist^5), data = tree_diversity))
summary(lm(Evenness ~ Dist + I(Dist^3) +  I(Dist^5), data = tree_diversity))
summary(lm(perc_planted ~ I(Dist^2), data = tree_diversity))
summary(lm(perc_native ~ Dist + I(Dist^2) +  I(Dist^5), data = tree_diversity))
# only Sum_stem ~ Dist etc. shows weakly significant p value
# regsubsets() for shrub
# plot
par(mfrow = c(2,2))
for (i in c("Sum_area", "Richness", "Shannon", "Evenness", "perc_planted", "perc_private", "perc_nonstreet", "perc_native")) {
  plot(regsubsets(shrub_diversity[[i]] ~ Dist + I(Dist^2) + I(Dist^3) + I(Dist^4) + I(Dist^5), 
                  data = shrub_diversity), scale = "adjr2", main = i)
}
# statistic analysis
summary(lm(Sum_area ~ I(Dist^2), data = shrub_diversity))
summary(lm(Evenness ~ I(Dist^2) + I(Dist^3) + I(Dist^4) + I(Dist^5), data = shrub_diversity))
summary(lm(perc_planted ~ I(Dist^2) + I(Dist^5), data = shrub_diversity))
summary(lm(perc_private ~ I(Dist^3) + I(Dist^4) + I(Dist^5), data = shrub_diversity))
summary(lm(perc_nonstreet ~ Dist + I(Dist^2), data = shrub_diversity))
summary(lm(perc_native ~ I(Dist^5), data = shrub_diversity))
# significant p: perc_planted ~ dist, perc_private ~ dist, perc_native ~ dist

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

# plot indexes ~ distance level
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

# regsubsets() for tree
# chose the best model for tree & shrub: perc_native% ~ distance 
par(mfrow = c(1,2))
plot(regsubsets(perc_native ~ Dist + I(Dist^2) + I(Dist^3) + I(Dist^4) + I(Dist^5), data = tree_diversity), scale = "adjr2")
plot(regsubsets(perc_native ~ Dist + I(Dist^2) + I(Dist^3) + I(Dist^4) + I(Dist^5), data = shrub_diversity), scale = "adjr2")
par(opar)

# statistic analysis for tree & shrub: perc_native% ~ distance 
summary(lm(perc_native ~ Dist + I(Dist^2) +  I(Dist^5), data = tree_diversity))
summary(lm(perc_native ~ I(Dist^5), data = shrub_diversity))

# plot of plant attr ~ distance with signifcant p value
curve_data <- data.frame(x = seq(0,12000,100))
curve_data$y <- 5.161e-01-1.521e-21*curve_data$x^5
ggplot(shrub_diversity, aes(Dist, perc_native)) + geom_point(alpha = 0.5) + geom_smooth(aes(x, y), data = curve_data, method = "loess", se = FALSE)
rm(curve_data)


# discussion: plot plant attr vs. land use
tree_data %>% group_by(Landuse_class) %>% 
  dplyr::summarise(perc = sum(ifelse(Pla_spo == "Planted", Stem, 0)/sum(Stem))) 
tree_data %>% group_by(Landuse_class) %>% 
  dplyr::summarise(perc = sum(ifelse(Pot == "Non_pot", Stem, 0)/sum(Stem)))
tree_data %>% group_by(Landuse_class) %>% 
  dplyr::summarise(perc = sum(ifelse(Pub_pri == "Private", Stem, 0)/sum(Stem)))
tree_data %>% group_by(Landuse_class) %>% 
  dplyr::summarise(perc = sum(ifelse(Street == "Non_street", Stem, 0)/sum(Stem)))
tree_data %>% group_by(Landuse_class) %>% 
  dplyr::summarise(perc = sum(ifelse(Nt_ex == "nt", Stem, 0)/sum(Stem)))


#
tree_data %>% group_by(Landuse_class) %>% 
  dplyr::summarise(total = n(), perc_street = sum(Street == "Street")/n())



# discussion: plot plant attr vs. land ownership
tree_data %>% group_by(Landscaping) %>% 
  dplyr::summarise(perc = sum(ifelse(Pla_spo == "Planted", Stem, 0)/sum(Stem))) 
tree_data %>% group_by(Landscaping) %>% 
  dplyr::summarise(perc = sum(ifelse(Pot == "Non_pot", Stem, 0)/sum(Stem)))
tree_data %>% group_by(Landscaping) %>% 
  dplyr::summarise(perc = sum(ifelse(Pub_pri == "Private", Stem, 0)/sum(Stem)))
tree_data %>% group_by(Landscaping) %>% 
  dplyr::summarise(perc = sum(ifelse(Street == "Non_street", Stem, 0)/sum(Stem)))
tree_data %>% group_by(Landscaping) %>% 
  dplyr::summarise(perc = sum(ifelse(Nt_ex == "nt", Stem, 0)/sum(Stem)))


# discussion: plot plant attr vs. distance level
tree_data %>% group_by(Dist_level) %>% 
  dplyr::summarise(perc = sum(ifelse(Pla_spo == "Planted", Stem, 0)/sum(Stem))) 
tree_data %>% group_by(Dist_level) %>% 
  dplyr::summarise(perc = sum(ifelse(Pot == "Non_pot", Stem, 0)/sum(Stem)))
tree_data %>% group_by(Dist_level) %>% 
  dplyr::summarise(perc = sum(ifelse(Pub_pri == "Private", Stem, 0)/sum(Stem)))
tree_data %>% group_by(Dist_level) %>% 
  dplyr::summarise(perc = sum(ifelse(Street == "Non_street", Stem, 0)/sum(Stem)))
tree_data %>% group_by(Dist_level) %>% 
  dplyr::summarise(perc = sum(ifelse(Nt_ex == "nt", Stem, 0)/sum(Stem)))

# following mess ...
# native species % along urban gradient? 
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

tree_absent <- subset(tree_data, select = c("Plot_ID", "Species_CN", "Stem")) %>% 
  pivot_wider(names_from = Species_CN, values_from = Stem, 
              values_fn = list(Stem = mean), values_fill = list(Stem = 0)) %>% 
  mutate(Richness = apply(.[2:ncol(.)]>0, 1, sum)) %>% 
  left_join(all_plot_info, by = "Plot_ID")

tree_absent <- subset(tree_data, select = c("Plot_ID", "Species_CN", "Nt_ex")) %>% unique() %>% 
  group_by(Plot_ID) %>% dplyr::summarise(perc_native = sum(Nt_ex == "nt")/n()) %>% 
  left_join(all_plot_info, by = "Plot_ID")

  as.data.frame() %>% 



tree_diversity_perc_native <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Nt_ex == "nt", Stem, 0)/sum(Stem)))