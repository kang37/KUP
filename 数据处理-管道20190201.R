setwd("C:/Rdata/KUP")
library(tidyverse)
library(Rmisc)
library(vegan)
library(digest)
opar <- par(no.readonly = TRUE)

# get the data
plot_info <- read.csv("In_Plot_info.csv", stringsAsFactors = FALSE) %>%
  mutate(Dist_level = ifelse(Dist < 3000, "Center", 
                              ifelse(Dist < 6000, "Cen-Mid", 
                                     ifelse(Dist < 9000, "Mid_Sub", "Suburban"))))
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
                                     ifelse(Plot_pub_pri ==0, "Private area", "Mix area")))
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
  ggplot(mds_plot_data, aes_string(x = "MDS1", y = "MDS2", color = atr)) + geom_point()
}
# for trees 
mds_plot_ID <- tree_diversity$Plot_ID[rowSums(tree_diversity[2:143])>2]
mds_plot_data <- tree_diversity %>% select(2:143) %>%
  subset(., rowSums(.)>2) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
mds_plot_data$stress
stressplot(mds_plot_data)
mds_plot_data <- mds_plot_data$points %>%
  as.data.frame() %>%
  mutate(Plot_ID = mds_plot_ID) %>% 
  left_join(plot_info, by = "Plot_ID")
mds_plot <- mapply(.myplot, atr = c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level"), SIMPLIFY = FALSE)
Rmisc::multiplot(plotlist = mds_plot, layout = matrix(1:4, nrow = 2))

# pairwise ANOSIM analysis data of trees for past
past <- tree_diversity %>% 
  select(2:143) %>% 
  subset(., rowSums(.)>2) %>% 
  mutate(Plot_ID = mds_plot_ID) %>%
  left_join(plot_info, by = "Plot_ID")
# general result of ANOSIM
par(mfrow = c(2,2))
for (i in c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level")) {
  x <- anosim(as.data.frame(past)[,c(1:142)], as.data.frame(past)[,i], 
              distance = "bray", permutations = 999)
  cat(i, "R=", x$statistic, "p=", x$signif, "\n")
  plot(x)
  rm(x)
}
par(opar)

# pairwise ANOSIM: the result almost = PAST 
for (i in c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level")) {
  x <- combn(unique(as.data.frame(past)[,i]), 2)
  for (j in 1:ncol(x)) {
    y <- past %>% as.data.frame() %>% filter(past[,i] == x[1,j] | past[,i] == x[2,j]) 
    z <- anosim(y[1:142], y[,i])
    if (z$signif < 0.05) cat(x[1,j], x[2,j], z$signif, "\n", sep = "  ")
  }
}

# and for shrubs
mds_plot_ID <- shrub_diversity$Plot_ID[rowSums(shrub_diversity[2:196]) > 0]
mds_plot_data <- shrub_diversity %>% select(2:196) %>%
  subset(., rowSums(.) > 0) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE)
mds_plot_data$stress
stressplot(mds_plot_data)
mds_plot_data <- mds_plot_data$points %>%
  as.data.frame() %>%
  mutate(Plot_ID = mds_plot_ID) %>% 
  left_join(plot_info, by = "Plot_ID")
mds_plot <- mapply(.myplot, atr = c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level"), SIMPLIFY = FALSE)
Rmisc::multiplot(plotlist = mds_plot, layout = matrix(1:4, nrow = 2))
# pairwise ANOSIM analysis data of shrubs for past
past <- shrub_diversity %>% 
  select(2:196) %>% 
  subset(., rowSums(.)>2) %>% 
  mutate(Plot_ID = mds_plot_ID) %>%
  left_join(plot_info, by = "Plot_ID")
# general result of ANOSIM
par(mfrow = c(2,2))
for (i in c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level")) {
  x <- anosim(as.data.frame(past)[,c(1:142)], as.data.frame(past)[,i], 
              distance = "bray", permutations = 999)
  cat(i, "R=", x$statistic, "p=", x$signif, "\n")
  plot(x)
  rm(x)
}
par(opar)













