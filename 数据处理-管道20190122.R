setwd("C:/Rdata/KUP")
library(tidyverse)
library(vegan)
library(digest)

plot_info <- read.csv("In_Plot_info.csv", stringsAsFactors = FALSE) %>%
  mutate(Dist_level = ifelse(Dist < 3000, "Center", 
                              ifelse(Dist < 6000, "Cen-Mid", 
                                     ifelse(Dist < 9000, "Mid_Sub", "Suburban"))))
plant_info <- read.csv("In_plant_info.csv", stringsAsFactors = FALSE)
plant_data <- read.csv("In_plant_data.csv", stringsAsFactors = FALSE) %>%
  filter(Plot_ID != 67) %>% 
  left_join(plant_info, by = "Species_CN") %>%
  left_join(plot_info,by = "Plot_ID")
plot_info <- plant_data %>% group_by(Plot_ID) %>% 
  summarise(Plot_pub_pri = sum(ifelse(Pub_pri == "Public", 1, 0)/n())) %>%
  left_join(plot_info, by = "Plot_ID") %>%
  mutate(Plot_pub_pri_level = ifelse(Plot_pub_pri ==1, "Public area", 
                                     ifelse(Plot_pub_pri ==0, "Private area", "Mix area")))
shrub_data <- plant_data %>% 
  filter(Tree_shrub == "Shrub") %>%
  mutate(Stem = NULL)
tree_diversity <- plant_data %>% 
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
shrub_diversity <- plant_data %>% 
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

# Top species by species and stem or area
plant_info %>% group_by(Family) %>%
  summarise(Number =n(), Prop = n()/nrow(plant_info)) %>% 
  arrange(desc(Prop))
# of trees
plant_data %>% subset(Tree_shrub == "Tree") %>%
  group_by(Family) %>%
  summarise(Number =n(), Prop = n()/nrow(plant_info)) %>% 
  arrange(desc(Prop))
# of shrubs
plant_data %>% subset(Tree_shrub == "Shrub") %>%
  group_by(Family) %>%
  summarise(Prop = sum(Area)/sum(subset(plant_data, Tree_shrub =="Shrub")$Area)) %>% 
  arrange(desc(Prop))

# mds graph for
# trees 
mds_plot_ID <- tree_diversity$Plot_ID[rowSums(tree_diversity[2:143])>2]
mds_plot_data <- tree_diversity %>% select(2:143) %>%
  subset(., rowSums(.)>2) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) %>%
  .$points %>%
  as.data.frame() %>%
  mutate(Plot_ID = mds_plot_ID) %>% 
  left_join(plot_info, by = "Plot_ID")
.myplot <- function(atr) {
  ggplot(mds_plot_data, aes_string(x = "MDS1", y = "MDS2", color = atr)) + geom_point()
}
mds_plot <- mapply(.myplot, atr = c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level"), SIMPLIFY = FALSE)
Rmisc::multiplot(plotlist = mds_plot, layout = matrix(1:4, nrow = 2))

for (i in c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level")) {
  x <- adonis(as.data.frame(tree_diversity)[,c(2:143)] ~ as.data.frame(tree_diversity)[,i], 
         method = "bray", permutations = 999)
  cat(i, x$aov.tab$`Pr(>F)`[1], "\n")
  rm(x)
}

# the data for PAST analysis
write.csv(tree_diversity %>% 
            select(2:143) %>% 
            subset(., rowSums(.)>2) %>% 
            mutate(Plot_ID = mds_plot_ID) %>%
            left_join(plot_info, by = "Plot_ID"), "past.csv")

# and for shrubs
mds_plot_ID <- shrub_diversity$Plot_ID[rowSums(shrub_diversity[2:196])>2]
mds_plot_data <- shrub_diversity %>% select(2:196) %>%
  subset(., rowSums(.)>2) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) %>%
  .$points %>%
  as.data.frame() %>%
  mutate(Plot_ID = mds_plot_ID) %>% 
  left_join(plot_info, by = "Plot_ID")
.myplot <- function(atr) {
  ggplot(mds_plot_data, aes_string(x = "MDS1", y = "MDS2", color = atr)) + geom_point()
}
mds_plot <- mapply(.myplot, atr = c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level"), SIMPLIFY = FALSE)
Rmisc::multiplot(plotlist = mds_plot, layout = matrix(1:4, nrow = 2))

for (i in c("Ward", "Dist_level", "Landuse_class", "Plot_pub_pri_level")) {
  x <- adonis(as.data.frame(shrub_diversity)[,c(2:196)] ~ as.data.frame(shrub_diversity)[,i], 
              method = "bray", permutations = 999)
  cat(i, x$aov.tab$`Pr(>F)`[1], "\n")
  rm(x)
}

# pic of the plants
par(mfrow= c(2,5))
for (i in c("Pla_spo", "Pot", "Pub_pri", "Street", "Nt_ex")) {
  barplot(table(subset(plant_data, Tree_shrub == "Tree")[,i]), 
          ylim = c(0, 1200), cex.axis = 3, cex.names = 3)
}

barplot(c("Planted" = sum(filter(shrub_data, shrub_data$Pla_spo ==  "Planted")$Area), 
        "Spontaneous" = sum(filter(shrub_data, shrub_data$Pla_spo ==  "Spontaneous")$Area)), 
        ylim = c(0, 1000), cex.axis = 3, cex.names = 3)
barplot(c("Non-Pot" = sum(filter(shrub_data, shrub_data$Pot ==  "Non_pot")$Area), 
          "Pot" = sum(filter(shrub_data, shrub_data$Pot ==  "Pot")$Area)), 
        ylim = c(0, 1000), cex.axis = 3, cex.names = 3)
barplot(c("Private" = sum(filter(shrub_data, shrub_data$Pub_pri ==  "Private")$Area), 
          "Public" = sum(filter(shrub_data, shrub_data$Pub_pri ==  "Public")$Area)), 
        ylim = c(0, 1000), cex.axis = 3, cex.names = 3)
barplot(c("Non-street" = sum(filter(shrub_data, shrub_data$Street ==  "Non_street")$Area), 
          "Street" = sum(filter(shrub_data, shrub_data$Street ==  "Street")$Area)), 
        ylim = c(0, 1000), cex.axis = 3, cex.names = 3) 
barplot(c("ex" = sum(filter(shrub_data, shrub_data$Nt_ex ==  "ex")$Area), 
          "nt" = sum(filter(shrub_data, shrub_data$Nt_ex ==  "nt")$Area)), 
        ylim = c(0, 1000), cex.axis = 3, cex.names = 3) 

# or use ggplot2
# p <- list()
# for (i in c("Pla_spo", "Pot", "Pub_pri")) {
#   p <- c(p, 
#          list(ggplot(shrub_data, aes_string(i)) + geom_bar() +
#                 theme(axis.title.x = element_text(size = 10))))
# }
# Rmisc::multiplot(plotlist = p, layout = matrix(1:4, ncol = 2))
# unsolved: for shrubs? 






