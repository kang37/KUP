setwd("C:/Rdata/KUP")
library(tidyverse)
library(Rmisc)
library(vegan)
library(PerformanceAnalytics)
library(digest)
library(car)
library(dunn.test)
opar <- par(no.readonly = TRUE)

# define the factor levels
Ward_faclev <- c("Uky身-ku", "Saky身-ku", "Kita-ku", "Kamigy身-ku", "Nakagy身-ku", "Shimogy身-ku", "Higashiyama-ku", "Yamashina-ku", "Fushimi-ku", "Minami-ku", "Nishiky身-ku")
Dist_level_faclev <- c("Center", "Cen-Mid", "Mid_Sub", "Suburban")
Landuse_class_faclev <- c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")
Plot_pub_pri_level_faclev <- c("Private area", "Public area")

# get the data
plot_info <- read.csv("In_Plot_info.csv", stringsAsFactors = FALSE) %>%
  mutate(Dist_level = NA)
plot_info$Dist_level[plot_info$Dist < 2000] <- "1"
plot_info$Dist_level[plot_info$Dist >= 2000 & plot_info$Dist < 4000] <- "2"
plot_info$Dist_level[plot_info$Dist >= 4000 & plot_info$Dist < 6000] <- "3"
plot_info$Dist_level[plot_info$Dist >= 6000 & plot_info$Dist < 8000] <- "4"
plot_info$Dist_level[plot_info$Dist >= 8000 & plot_info$Dist < 10000] <- "5"
plot_info <- 
  plot_info %>% mutate(Dist_level = factor(Dist_level, levels = c("1", "2", "3", "4", "5"))) %>% 
  mutate(Ward = factor(Ward, levels = Ward_faclev), 
         Landuse_class = factor(Landuse_class, levels = Landuse_class_faclev))
plant_info <- read.csv("In_plant_info.csv", stringsAsFactors = FALSE)
all_plant_data <- read.csv("In_plant_data.csv", stringsAsFactors = FALSE) %>%
  filter(Plot_ID != 67) %>% 
  left_join(plant_info, by = "Species_CN") %>%
  left_join(plot_info,by = "Plot_ID")
tree_data <- all_plant_data %>% 
  subset(Tree_shrub == "Tree")
shrub_data <- all_plant_data %>% 
  subset(Tree_shrub == "Shrub")

tree_plant_info <- all_plant_data %>% filter(Tree_shrub == "Tree") %>% 
  select(Species_CN) %>% 
  unique() %>% 
  left_join(plant_info, by = "Species_CN")
shrub_plant_info <- all_plant_data %>% filter(Tree_shrub == "Shrub") %>% 
  select(Species_CN) %>% 
  unique() %>% 
  left_join(plant_info, by = "Species_CN")

plot_info <- all_plant_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(Plot_pub_pri = sum(ifelse(Pub_pri == "Public", 1, 0)/n())) %>%
  left_join(plot_info, by = "Plot_ID") %>%
  mutate(Plot_pub_pri_level = ifelse(Plot_pub_pri ==1, "Public area", 
                                     ifelse(Plot_pub_pri ==0, "Private area", NA))) %>% 
  mutate(Plot_pub_pri_level = factor(Plot_pub_pri_level, levels = Plot_pub_pri_level_faclev))
tree_diversity <- all_plant_data %>% 
  filter(Tree_shrub == "Tree") %>% 
  select(Plot_ID, Species_CN, Stem) %>%
  pivot_wider(names_from = Species_CN, values_from = Stem, 
              values_fn = list(Stem = sum), values_fill = list(Stem = 0)) %>% 
  mutate(Sum_stem = rowSums(.[2:ncol(.)]), 
         Richness = apply(.[2:ncol(.)]>0, 1, sum),
         Shannon = diversity(.[2:ncol(.)], index = "shannon"), 
         Simpson = diversity(.[2:ncol(.)], index = "simpson"),
         Evenness = Shannon / log(Richness)) %>%
  left_join(plot_info,by = "Plot_ID")
shrub_diversity <- all_plant_data %>% 
  filter(Tree_shrub == "Shrub") %>% 
  select(Plot_ID, Species_CN, Area) %>%
  pivot_wider(names_from = Species_CN, values_from = Area, 
              values_fn = list(Area = sum), values_fill = list(Area = 0)) %>% 
  mutate(Sum_area = rowSums(.[2:ncol(.)]), 
         Richness = apply(.[2:ncol(.)]>0, 1, sum),
         Shannon = diversity(.[2:ncol(.)], index = "shannon"), 
         Simpson = diversity(.[2:ncol(.)], index = "simpson"),
         Evenness = Shannon / log(Richness)) %>%
  left_join(plot_info,by = "Plot_ID")

# analysis begin
# the total species, genera and families of all plants
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

# top species by species and number or area
# by species
# of all plants
plant_info %>% group_by(Family) %>%
  summarise(Species =n(), Prop = n()/nrow(plant_info)) %>% 
  arrange(desc(Prop)) %>% 
  print() %>% 
  ggplot(aes(reorder(Family, -Prop), Prop)) + geom_col()
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


# attributes of the plants
# the exotic vs. native regarding species
plant_info %>% group_by(Nt_ex) %>% 
  dplyr::summarise(n()/nrow(plant_info))

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

# mds graph

# build the function to draw the plot
.myplot <- function(atr) {
  ggplot(mds_plot_point_data, aes_string(x = "MDS1", y = "MDS2", color = atr)) + geom_point()
}

# nmds plot for trees 
set.seed(1234)
mds_plot_ID <- tree_diversity$Plot_ID[!(tree_diversity$Plot_ID %in% c(214, 261, 244, 313))]
mds_plot_selected <- tree_diversity %>% filter(Plot_ID %in% mds_plot_ID)
mds_plot_data <- mds_plot_selected %>% 
  select(2:143) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
mds_plot_data$stress
stressplot(mds_plot_data)
mds_plot_point_data <- mds_plot_data$points %>%
  as.data.frame() %>%
  mutate(Plot_ID = mds_plot_ID) %>% 
  left_join(plot_info, by = "Plot_ID")
mds_plot <- mapply(.myplot, atr = c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level"), SIMPLIFY = FALSE)
Rmisc::multiplot(plotlist = mds_plot, layout = matrix(1:4, nrow = 2))

# general result of ANOSIM of trees
par(mfrow = c(2,2))
for (i in c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level")) {
  x <- anosim(as.data.frame(mds_plot_selected)[,c(2:143)], as.data.frame(mds_plot_selected)[,i], 
              distance = "bray", permutations = 999)
  cat(i, "R=", x$statistic, "p=", x$signif, "\n")
  plot(x)
  rm(x)
}
par(opar)

# pairwise result of ANOSIM of trees
x <- vector("list")
for (i in c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level")) {
  set.seed(1234)
  a <- combn(levels(factor(as.data.frame(mds_plot_selected)[,i])), 2)
  z <- vector("list", 3)
  for (j in 1:ncol(a)) {
    b <- mds_plot_selected %>% as.data.frame() %>% 
      filter(mds_plot_selected[[i]] == a[1,j] | mds_plot_selected[[i]] == a[2,j]) 
    c <- anosim(b[2:143], b[,i])
    z[[1]] <- c(z[[1]], as.character(a[1,j]))
    z[[2]] <- c(z[[2]], as.character(a[2,j]))
    z[[3]] <- c(z[[3]], c$signif)
  }
  y <- data.frame(comp_1 = z[[1]], comp_2 = z[[2]], p = z[[3]]) %>% 
    mutate(comp_1 = factor(comp_1, levels = levels(factor(as.data.frame(mds_plot_selected)[,i]))), 
           comp_2 = factor(comp_2, levels(factor(as.data.frame(mds_plot_selected)[,i]))))
  print(y)
  x <- c(x, list(ggplot(y, aes(comp_1, comp_2)) + geom_tile(aes(fill = p)) + 
                   scale_fill_gradient(high = "orange", low = "red", limit = c(0, 0.01)) + 
                   theme(axis.text.x = element_text(angle = 90))))
}
Rmisc::multiplot(plotlist = x, layout = matrix(1:4, ncol = 2))
rm(x, y, z)

# and nmds plot for shrubs
mds_plot_ID <- shrub_diversity$Plot_ID[!(shrub_diversity$Plot_ID %in% c(269, 214, 75, 164, 244))]
mds_plot_selected <- shrub_diversity %>% filter(Plot_ID %in% mds_plot_ID)
mds_plot_data <- mds_plot_selected %>% 
  select(2:196) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
mds_plot_data$stress
stressplot(mds_plot_data)
mds_plot_point_data <- mds_plot_data$points %>%
  as.data.frame() %>%
  mutate(Plot_ID = mds_plot_ID) %>% 
  left_join(plot_info, by = "Plot_ID")
mds_plot <- mapply(.myplot, atr = c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level"), SIMPLIFY = FALSE)
Rmisc::multiplot(plotlist = mds_plot, layout = matrix(1:4, nrow = 2))

# general result of ANOSIM
par(mfrow = c(2,2))
for (i in c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level")) {
  x <- anosim(as.data.frame(mds_plot_selected)[,c(2:143)], as.data.frame(mds_plot_selected)[,i], 
              distance = "bray", permutations = 999)
  cat(i, "R=", x$statistic, "p=", x$signif, "\n")
  plot(x)
  rm(x)
}
par(opar)
# pairwise result of ANOSIM of shrubs
x <- vector("list")
for (i in c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level")) {
  set.seed(1234)
  a <- combn(levels(factor(as.data.frame(mds_plot_selected)[,i])), 2)
  z <- vector("list", 3)
  for (j in 1:ncol(a)) {
    b <- mds_plot_selected %>% as.data.frame() %>% 
      filter(mds_plot_selected[[i]] == a[1,j] | mds_plot_selected[[i]] == a[2,j]) 
    c <- anosim(b[2:196], b[,i])
    z[[1]] <- c(z[[1]], as.character(a[1,j]))
    z[[2]] <- c(z[[2]], as.character(a[2,j]))
    z[[3]] <- c(z[[3]], c$signif)
  }
  y <- data.frame(comp_1 = z[[1]], comp_2 = z[[2]], p = z[[3]]) %>% 
    mutate(comp_1 = factor(comp_1, levels = levels(factor(as.data.frame(mds_plot_selected)[,i]))), 
           comp_2 = factor(comp_2, levels(factor(as.data.frame(mds_plot_selected)[,i]))))
  print(y)
  x <- c(x, list(ggplot(y, aes(comp_1, comp_2)) + geom_tile(aes(fill = p)) + 
                   scale_fill_gradient(high = "orange", low = "red", limit = c(0, 0.01)) + 
                   theme(axis.text.x = element_text(angle = 90))))
}
Rmisc::multiplot(plotlist = x, layout = matrix(1:4, ncol = 2))
rm(x, y, z)

# cor among the indexes
chart.Correlation(subset(tree_diversity, select = c("Sum_stem", "Richness", "Shannon", "Simpson", "Evenness")))
chart.Correlation(subset(shrub_diversity, select = c("Sum_area", "Richness", "Shannon", "Simpson", "Evenness")))


# kruskal test & boxplot for trees
{
  pvalue_list <- vector("list", 3)
  for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
    for (j in c("Ward", "Landuse_class", "Plot_pub_pri_level")) {
      pvalue_list[[1]] <- c(pvalue_list[[1]], i)
      pvalue_list[[2]] <- c(pvalue_list[[2]], j)
      pvalue_list[[3]] <- c(pvalue_list[[3]], round(kruskal.test(as.data.frame(tree_diversity)[, i] ~ as.data.frame(tree_diversity)[, j])$p.value,digits = 3))
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

tree_diversity %>% subset(select = c("Sum_stem", "Richness", "Shannon", "Simpson", "Evenness", 
                                     "Ward", "Landuse_class", "Plot_pub_pri_level")) %>% 
  pivot_longer(cols = c("Sum_stem", "Richness", "Shannon", "Evenness"), 
               names_to = "index", values_to = "index_value") %>% 
  pivot_longer(cols = c("Ward", "Landuse_class", "Plot_pub_pri_level"), 
               names_to = "attr", values_to = "attr_value") %>% 
  mutate(index = factor(index, levels = c("Sum_stem", "Richness", "Shannon", "Evenness")), 
         attr = factor(attr, levels = c("Ward", "Landuse_class", "Plot_pub_pri_level")), 
         attr_value = factor(attr_value, levels = c(Ward_faclev, Landuse_class_faclev, Plot_pub_pri_level_faclev))) %>% 
  na.omit() %>%
  ggplot(aes(attr_value, index_value)) + geom_boxplot() + 
  facet_grid(index ~ attr, scales = "free", space = "free_x", switch = "both") + 
  scale_y_continuous(expand = expand_scale(mult = c(0.05,0.3))) +
  geom_text(data = pvalue, aes(x =Inf, y = Inf, label = label), size=3.5, hjust = 1.05, vjust = 1.5) +
  theme(axis.text = element_text(angle = 90))
# delete the vars
rm(pvalue, pvalue_list)

# dunn test & pairwise of trees
for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
  pvalue_pairwise_list <- vector("list", 2)
  for (j in c("Ward", "Landuse_class", "Plot_pub_pri_level")) {
    pvalue_pairwise_list[[1]] <- c(pvalue_pairwise_list[[1]], dunn.test(as.data.frame(tree_diversity)[, i], as.data.frame(tree_diversity)[, j])$comparisons)
    pvalue_pairwise_list[[2]] <- c(pvalue_pairwise_list[[2]], dunn.test(as.data.frame(tree_diversity)[, i], as.data.frame(tree_diversity)[, j])$P.adjusted)
  }
  pvalue_pairwise_df_up <- 
    data.frame(comparison = pvalue_pairwise_list[[1]], p = pvalue_pairwise_list[[2]]) %>% 
    separate(comparison, into = c("comparison_1", "comparison_2"), sep = " - ")
  pvalue_pairwise_df_down <- data.frame(comparison_1 = pvalue_pairwise_df_up$comparison_2, 
                                           comparison_2 = pvalue_pairwise_df_up$comparison_1, 
                                           p = pvalue_pairwise_df_up$p)
  pvalue_list_pairwise_df <- rbind(pvalue_pairwise_df_up, pvalue_pairwise_df_down) %>% 
    mutate(comparison_1 = factor(comparison_1, levels = c(Ward_faclev, Landuse_class_faclev, Plot_pub_pri_level_faclev)), 
           comparison_2 = factor(comparison_2, levels = c(Ward_faclev, Landuse_class_faclev, Plot_pub_pri_level_faclev)))
  print(ggplot(pvalue_list_pairwise_df, aes(comparison_1, comparison_2, fill = p)) + geom_tile() +
          geom_text(aes(label = round(p, digits = 2)), size = 3.5) +
          scale_fill_gradient2(high = "blue", low = "red", midpoint = 0.05, limits = c(0, 0.05)) + 
          theme(axis.text.x = element_text(angle = 90)) + 
          labs(title = i)) + xlab(NULL) + ylab(NULL)
}
# delete the vars
rm(pvalue_pairwise_list, pvalue_pairwise_df_up, pvalue_pairwise_df_down, pvalue_list_pairwise_df)

# kruskal test & boxplot for shrubs
{
  pvalue_list <- vector("list", 3)
  for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
    for (j in c("Ward", "Landuse_class", "Plot_pub_pri_level")) {
      pvalue_list[[1]] <- c(pvalue_list[[1]], i)
      pvalue_list[[2]] <- c(pvalue_list[[2]], j)
      pvalue_list[[3]] <- c(pvalue_list[[3]], round(kruskal.test(as.data.frame(shrub_diversity)[, i] ~ as.data.frame(shrub_diversity)[, j])$p.value,digits = 3))
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

shrub_diversity %>% subset(select = c("Sum_area", "Richness", "Shannon", "Evenness", 
                                      "Ward", "Landuse_class", "Plot_pub_pri_level")) %>% 
  pivot_longer(cols = c("Sum_area", "Richness", "Shannon", "Evenness"), 
               names_to = "index", values_to = "index_value") %>% 
  pivot_longer(cols = c("Ward", "Landuse_class", "Plot_pub_pri_level"), 
               names_to = "attr", values_to = "attr_value") %>% 
  mutate(index = factor(index, levels = c("Sum_area", "Richness", "Shannon","Evenness")), 
         attr = factor(attr, levels = c("Ward", "Landuse_class", "Plot_pub_pri_level")), 
         attr_value = factor(attr_value, levels = c(Ward_faclev, Landuse_class_faclev, Plot_pub_pri_level_faclev))) %>% 
  na.omit() %>%
  ggplot(aes(attr_value, index_value)) + geom_boxplot() + 
  facet_grid(index ~ attr, scales = "free", space = "free_x", switch = "both") + 
  scale_y_continuous(expand = expand_scale(mult = c(0.05,0.3))) +
  geom_text(data = pvalue, aes(x =Inf, y = Inf, label = label), size=3.5, hjust = 1.05, vjust = 1.5) +
  theme(axis.text = element_text(angle = 90))
# delete the vars
rm(pvalue, pvalue_list)

# dunn test & pairwise of shrubs
for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
  pvalue_pairwise_list <- vector("list", 2)
  for (j in c("Ward", "Landuse_class", "Plot_pub_pri_level")) {
    pvalue_pairwise_list[[1]] <- c(pvalue_pairwise_list[[1]], dunn.test(as.data.frame(shrub_diversity)[, i], as.data.frame(shrub_diversity)[, j])$comparisons)
    pvalue_pairwise_list[[2]] <- c(pvalue_pairwise_list[[2]], dunn.test(as.data.frame(shrub_diversity)[, i], as.data.frame(shrub_diversity)[, j])$P.adjusted)
  }
  pvalue_pairwise_df_up <- 
    data.frame(comparison = pvalue_pairwise_list[[1]], p = pvalue_pairwise_list[[2]]) %>% 
    separate(comparison, into = c("comparison_1", "comparison_2"), sep = " - ")
  pvalue_pairwise_df_down <- data.frame(comparison_1 = pvalue_pairwise_df_up$comparison_2, 
                                        comparison_2 = pvalue_pairwise_df_up$comparison_1, 
                                        p = pvalue_pairwise_df_up$p)
  pvalue_list_pairwise_df <- rbind(pvalue_pairwise_df_up, pvalue_pairwise_df_down) %>% 
    mutate(comparison_1 = factor(comparison_1, levels = c(Ward_faclev, Landuse_class_faclev, Plot_pub_pri_level_faclev)), 
           comparison_2 = factor(comparison_2, levels = c(Ward_faclev, Landuse_class_faclev, Plot_pub_pri_level_faclev)))
  print(ggplot(pvalue_list_pairwise_df, aes(comparison_1, comparison_2, fill = p)) + geom_tile() +
          geom_text(aes(label = round(p, digits = 2)), size = 3.5) +
          scale_fill_gradient2(high = "blue", low = "red", midpoint = 0.05, limits = c(0, 0.05)) + 
          theme(axis.text.x = element_text(angle = 90)) + 
          labs(title = i)) + xlab(NULL) + ylab(NULL)
}
# delete the vars
rm(pvalue_pairwise_list, pvalue_pairwise_df_up, pvalue_pairwise_df_down, pvalue_list_pairwise_df)


# test linear models for indexes ~ distance
# y ~ x for trees
for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
  print(i)
  print(round(summary(lm(tree_diversity[[i]] ~ tree_diversity$Dist))$r.squared*100, digits = 2))
}
# y ~ x + x^2 for trees
for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
  print(i)
  print(round(summary(lm(tree_diversity[[i]] ~ tree_diversity$Dist + I(tree_diversity$Dist^2)))$r.squared*100, digits = 2))
}
# y ~ x + x^2 + x^3 for trees
for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
  print(i)
  print(round(summary(lm(tree_diversity[[i]] ~ tree_diversity$Dist + I(tree_diversity$Dist^2) + I(tree_diversity$Dist^3)))$r.squared*100, digits = 2))
}

# y ~ x for shrubs
for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
  print(i)
  print(round(summary(lm(shrub_diversity[[i]] ~ shrub_diversity$Dist))$r.squared*100, digits = 2))
}
# y ~ x + x^2 for shrubs
for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
  print(i)
  print(round(summary(lm(shrub_diversity[[i]] ~ shrub_diversity$Dist + I(shrub_diversity$Dist^2)))$r.squared*100, digits = 2))
}
# y ~ x + x^2 + x^3 for shrubs
for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
  print(i)
  print(round(summary(lm(shrub_diversity[[i]] ~ shrub_diversity$Dist + I(shrub_diversity$Dist^2) + I(shrub_diversity$Dist^3)))$r.squared*100, digits = 2))
}

# percentage of exotic species ~ distance
tree_data_percex <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(percex = sum(Nt_ex == "ex")/n()) %>%
  left_join(plot_info, by = "Plot_ID") %>%
  as.data.frame()
ggplot(tree_data_percex, aes(Dist, percex)) + geom_point(aes(color = Landuse_class), size = 3) + geom_smooth(method = "lm")
ggplot(tree_data_percex, aes(Landuse_class, percex)) + geom_boxplot()
kruskal.test(tree_data_percex$percex ~ tree_data_percex$Landuse_class)
summary(lm(tree_data_percex$percex ~ tree_data_percex$Dist))

# dist
ggplot(subset(tree_diversity, Landuse_class == "Com"), aes(Dist)) + geom_histogram(binwidth = 1000) + xlim(-1000, 12000)
ggplot(subset(tree_diversity, Landuse_class == "R low"), aes(Dist)) + geom_histogram(binwidth = 1000)

ggplot(tree_diversity, aes(Dist)) + geom_histogram(binwidth = 1000) + facet_wrap( ~ Landuse_class)

# cor between distance and indexes
for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
  print(i)
  print(cor.test(as.data.frame(tree_diversity)[, i], tree_diversity$Dist))
}
for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
  print(i)
  print(cor.test(as.data.frame(shrub_diversity)[, i], shrub_diversity$Dist))
}

for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
  print(i)
  print(cor.test(as.data.frame(subset(tree_diversity, Dist < 6000))[, i], subset(tree_diversity, Dist < 6000)$Dist))
}
# change something 
