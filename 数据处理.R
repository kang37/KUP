setwd("C:/Users/kangj/Documents/R/KUP")
library(Rmisc)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(PerformanceAnalytics)
library(dunn.test)
library(BiodiversityR)
library(codyn)
library(ggpubr)

## get data
# define default parameters
opar <- par(no.readonly = TRUE)
Land_use_type_faclev <- c("Com", "Com-neigh", "R-low", "R-high", "R-other", "Ind")

# information of all the plots
all_plot_info <- read.csv("In_plot_info.csv", stringsAsFactors = FALSE) %>% 
  mutate(Land_use_type = factor(Land_use_type, levels = Land_use_type_faclev))

# information of all the plants species: provenance and taxonomy 
all_plant_info <- read.csv("In_plant_info.csv", stringsAsFactors = FALSE)

# data of all the plants: taxonomy, attributes, abundance, etc.
all_plant_data <- read.csv("In_plant_data.csv", stringsAsFactors = FALSE) %>%
  left_join(all_plant_info, by = "Species_LT") %>%
  left_join(all_plot_info,by = "Plot_ID")

# data of trees and shrubs and diversity matrix of trees
tree_data <- subset(all_plant_data, Tree_shrub == "tree") %>% mutate(Area = NULL)
tree_diversity <- subset(tree_data, select = c("Plot_ID", "Species_LT", "Stem")) %>%
  pivot_wider(names_from = Species_LT, values_from = Stem, 
              values_fn = list(Stem = sum), values_fill = list(Stem = 0)) %>% 
  mutate(Density = rowSums(.[2:ncol(.)]), 
         Richness = specnumber(.[2:ncol(.)], MARGIN = 1),
         Shannon = diversity(.[2:ncol(.)], index = "shannon"), 
         Simpson = diversity(.[2:ncol(.)], index = "simpson"),
         Evenness = Shannon / log(Richness)) %>%
  left_join(all_plot_info, by = "Plot_ID")
tree_diversity_perc_planted <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pla_spo == "planted", Stem, 0)/sum(Stem))) 
tree_diversity_perc_nonpot <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pot == "non_pot", Stem, 0)/sum(Stem)))
tree_diversity_perc_private <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pub_pri == "private", Stem, 0)/sum(Stem)))
tree_diversity_perc_nonstreet <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Street == "non_street", Stem, 0)/sum(Stem)))
tree_diversity_perc_native <- tree_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Nt_ex == "native", Stem, 0)/sum(Stem)))
tree_diversity <- tree_diversity %>% mutate(
  perc_planted = tree_diversity_perc_planted$perc, 
  perc_nonpot = tree_diversity_perc_nonpot$perc, 
  perc_private = tree_diversity_perc_private$perc, 
  perc_nonstreet = tree_diversity_perc_nonstreet$perc, 
  perc_native = tree_diversity_perc_native$perc
)
rm(tree_diversity_perc_planted, tree_diversity_perc_nonpot, tree_diversity_perc_private, tree_diversity_perc_nonstreet, tree_diversity_perc_native)

# data of shrubs and diversity matrix of shrubs
shrub_data <- subset(all_plant_data, Tree_shrub == "shrub") %>% mutate(Stem = NULL)
shrub_diversity <- subset(shrub_data, select = c("Plot_ID", "Species_LT", "Area")) %>%
  pivot_wider(names_from = Species_LT, values_from = Area, 
              values_fn = list(Area = sum), values_fill = list(Area = 0)) %>% 
  mutate(Density = rowSums(.[2:ncol(.)]), 
         Richness = apply(.[2:ncol(.)]>0, 1, sum),
         Shannon = diversity(.[2:ncol(.)], index = "shannon"), 
         Simpson = diversity(.[2:ncol(.)], index = "simpson"),
         Evenness = Shannon / log(Richness)) %>%
  left_join(all_plot_info, by = "Plot_ID") 
shrub_diversity_perc_planted <- shrub_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pla_spo == "planted", Area, 0)/sum(Area))) 
shrub_diversity_perc_nonpot <- shrub_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pot == "non_pot", Area, 0)/sum(Area)))
shrub_diversity_perc_private <- shrub_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Pub_pri == "private", Area, 0)/sum(Area)))
shrub_diversity_perc_nonstreet <- shrub_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Street == "non_street", Area, 0)/sum(Area)))
shrub_diversity_perc_native <- shrub_data %>% group_by(Plot_ID) %>% 
  dplyr::summarise(perc = sum(ifelse(Nt_ex == "native", Area, 0)/sum(Area)))
shrub_diversity <- shrub_diversity %>% mutate(
  perc_planted = shrub_diversity_perc_planted$perc, 
  perc_nonpot = shrub_diversity_perc_nonpot$perc, 
  perc_private = shrub_diversity_perc_private$perc, 
  perc_nonstreet = shrub_diversity_perc_nonstreet$perc, 
  perc_native = shrub_diversity_perc_native$perc
)
rm(shrub_diversity_perc_planted, shrub_diversity_perc_nonpot, shrub_diversity_perc_private, shrub_diversity_perc_nonstreet, shrub_diversity_perc_native)

# tree diversity longer and shrub diversity longer data set
tree_diversity_long <- 
  subset(tree_diversity, 
         select = c("Density", "Richness", "Shannon", "Evenness","Land_use_type")) %>% 
  pivot_longer(cols = c("Density", "Richness", "Shannon", "Evenness"), 
               names_to = "index", values_to = "index_value") %>% 
  mutate(index = factor(index, levels = c("Density", "Richness", "Shannon", "Evenness")), 
         Land_use_type = factor(Land_use_type, levels = c(Land_use_type_faclev)))

shrub_diversity_long <- 
  subset(shrub_diversity, 
         select = c("Density", "Richness", "Shannon", "Evenness", "Land_use_type")) %>% 
  pivot_longer(cols = c("Density", "Richness", "Shannon", "Evenness"), 
               names_to = "index", values_to = "index_value") %>% 
  mutate(index = factor(index, levels = c("Density", "Richness", "Shannon","Evenness")), 
         Land_use_type = factor(Land_use_type, levels = c(Land_use_type_faclev)))

# some other variables 
number_tree_species <- length(unique(tree_data$Species_LT))
number_shrub_species <- length(unique(shrub_data$Species_LT))



### analysis begins
## general description
# the number of species
cat("total species:", length(unique(all_plant_data$Species_LT)), "\n", 
    "total genera:", length(unique(all_plant_data$Genus)), "\n", 
    "total families:", length(unique(all_plant_data$Family)), "\n", "\n", 
    "total species of trees:", number_tree_species, "\n", 
    "total species of shrubs:", number_shrub_species, "\n", 
    "common species of trees and shrubs:", length(intersect(
      unique(tree_data$Species_LT), unique(shrub_data$Species_LT))), "\n", 
    "species solely for trees:", length(setdiff(unique(tree_data$Species_LT), 
                                               unique(shrub_data$Species_LT))), "\n", 
    "Species solely for shrubs:", length(setdiff(unique(shrub_data$Species_LT), 
                                                 unique(tree_data$Species_LT))))

# top species families of all plants by species number
all_plant_info %>% group_by(Family) %>% 
  dplyr::summarise(Number = n(), Prop = n()/nrow(all_plant_info)) %>% 
  arrange(desc(Prop))

# the number of trees or area of shrubs
cat("number of trees:", nrow(tree_data), "in", nrow(tree_diversity), "plot", "\n", 
    "area of shrubs:", sum(shrub_data$Area), "m2 in", nrow(shrub_diversity), "plot")

# top species families of trees by individual
tree_data %>% 
  group_by(Family) %>% 
  dplyr::summarise(Number = sum(Stem), Prop = Number/sum(tree_data$Stem)) %>% 
  arrange(desc(Prop))

# top species families of shrubs by area
shrub_data %>% 
  group_by(Family) %>% 
  dplyr::summarise(SArea = sum(Area), Prop = SArea/sum(shrub_data$Area)) %>% 
  arrange(desc(Prop))



## attributes of trees and shrubs
# the number of exotic vs. native species
table(all_plant_info$Nt_ex)/nrow(all_plant_info)

# the attributes of trees
for (i in c("Pla_spo", "Pub_pri", "Nt_ex")) {
  print(tapply(tree_data$Stem, tree_data[,i], sum)/sum(tree_data$Stem), digits = 3)
  cat("\n")
}

# the attributes of shrubs
for (i in c("Pla_spo", "Pub_pri", "Nt_ex")) {
  print(tapply(shrub_data$Area, shrub_data[,i], sum)/sum(shrub_data$Area), digits = 3)
  cat("\n")
}



## Species accumulation curve 
accum_list <- vector("list", 4)
for (i in Land_use_type_faclev) {
  accum_ori <- specaccum(subset(tree_diversity, 
                                Land_use_type == i)[,2:(number_tree_species+1)])
  accum_list[[1]] <- c(accum_list[[1]], accum_ori$sites)
  accum_list[[2]] <- c(accum_list[[2]], accum_ori$richness)
  accum_list[[3]] <- c(accum_list[[3]], rep(i, length(accum_ori$sites)))
  accum_list[[4]] <- c(accum_list[[4]], rep("tree", length(accum_ori$sites)))
}
for (i in Land_use_type_faclev) {
  accum_ori <- specaccum(subset(shrub_diversity, 
                                Land_use_type == i)[,2:(number_shrub_species+1)])
  accum_list[[1]] <- c(accum_list[[1]], accum_ori$sites)
  accum_list[[2]] <- c(accum_list[[2]], accum_ori$richness)
  accum_list[[3]] <- c(accum_list[[3]], rep(i, length(accum_ori$sites)))
  accum_list[[4]] <- c(accum_list[[4]], rep("shrub", length(accum_ori$sites)))
}
accum_df <- data.frame(
  number_of_plot = accum_list[[1]], 
  number_of_species = accum_list[[2]], 
  Land_use_type = accum_list[[3]], 
  tree_shrub = accum_list[[4]]
) %>% mutate(Land_use_type = factor(Land_use_type, levels = Land_use_type_faclev), 
             tree_shrub = factor(tree_shrub, levels = c("tree", "shrub")))
ggplot(accum_df, aes(x = number_of_plot, y = number_of_species, color = Land_use_type)) + 
  geom_line(size = 1.5) + facet_wrap(~tree_shrub) + 
  theme(strip.text = element_text(size = 15), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 12)) +
  labs(x = "Number of quadrats", y = "Accumulative number of species") +
  scale_color_discrete("Land use types")
rm(accum_ori, accum_list, accum_df)



## rank abundance plot
# the rank abundance plot by land use types of trees
tree_rank_list <- vector("list", 6)
for (i in Land_use_type_faclev) {
  tree_rank_ori <- 
    subset(tree_diversity, Land_use_type == i)[2:(number_tree_species+1)] %>% 
    subset(select = colSums(.) != 0) %>% 
    as.data.frame() %>%
    rankabundance() %>%
    as.data.frame()
  tree_rank_list[[1]] <- c(tree_rank_list[[1]], rownames(tree_rank_ori))
  tree_rank_list[[2]] <- c(tree_rank_list[[2]], tree_rank_ori$rank)
  tree_rank_list[[3]] <- c(tree_rank_list[[3]], tree_rank_ori$rankfreq)
  tree_rank_list[[4]] <- c(tree_rank_list[[4]], tree_rank_ori$abundance)
  tree_rank_list[[5]] <- c(tree_rank_list[[5]], tree_rank_ori$proportion)
  tree_rank_list[[6]] <- c(tree_rank_list[[6]], rep(i, nrow(tree_rank_ori)))
}
tree_rank_df <- data.frame(
  Species_LT = tree_rank_list[[1]], 
  Rank = tree_rank_list[[2]], 
  Scaled_rank = tree_rank_list[[3]], 
  Abundance = tree_rank_list[[4]], 
  Proportion = tree_rank_list[[5]], 
  Land_use_type = tree_rank_list[[6]]
) %>% 
  left_join(all_plant_info[, c("Species_LT", "Nt_ex")], by = "Species_LT") %>% 
  mutate(Land_use_type = factor(Land_use_type, levels = Land_use_type_faclev))

# the rank abundance plot by land use types of shrubs
shrub_rank_list <- vector("list", 6)
for (i in Land_use_type_faclev) {
  shrub_rank_ori <- 
    subset(shrub_diversity, Land_use_type == i)[2:(number_shrub_species+1)] %>% 
    subset(select = colSums(.) != 0) %>% 
    as.data.frame() %>%
    rankabundance() %>%
    as.data.frame()
  shrub_rank_list[[1]] <- c(shrub_rank_list[[1]], rownames(shrub_rank_ori))
  shrub_rank_list[[2]] <- c(shrub_rank_list[[2]], shrub_rank_ori$rank)
  shrub_rank_list[[3]] <- c(shrub_rank_list[[3]], shrub_rank_ori$rankfreq)
  shrub_rank_list[[4]] <- c(shrub_rank_list[[4]], shrub_rank_ori$abundance)
  shrub_rank_list[[5]] <- c(shrub_rank_list[[5]], shrub_rank_ori$proportion)
  shrub_rank_list[[6]] <- c(shrub_rank_list[[6]], rep(i, nrow(shrub_rank_ori)))
}
shrub_rank_df <- data.frame(
  Species_LT = shrub_rank_list[[1]], 
  Rank = shrub_rank_list[[2]], 
  Scaled_rank = shrub_rank_list[[3]], 
  Abundance = shrub_rank_list[[4]], 
  Proportion = shrub_rank_list[[5]], 
  Land_use_type = shrub_rank_list[[6]]
) %>% 
  left_join(all_plant_info[, c("Species_LT", "Nt_ex")], by = "Species_LT") %>% 
  mutate(Land_use_type = factor(Land_use_type, levels = Land_use_type_faclev))

# rearrange and plot
ggarrange(ggplot(tree_rank_df, aes(Scaled_rank, log(Abundance), label = Species_LT)) + 
            geom_line() + 
            geom_point(aes(color = Nt_ex), alpha = 0.3, size = 2) + 
            facet_wrap(~Land_use_type, nrow = 1) + labs(title = "(a)") + 
            labs(x = "Scaled rank of species", y = "Log (abundance)") +
            scale_color_discrete("Provenance") +
            theme(strip.text = element_text(size = 12),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 10),
                  legend.title = element_text(size = 12),
                  legend.text = element_text(size = 12)),
          ggplot(shrub_rank_df, aes(Scaled_rank, log(Abundance), label = Species_LT)) + 
            geom_line() + 
            geom_point(aes(color = Nt_ex), alpha = 0.3, size = 2) +
            facet_wrap(~Land_use_type, nrow = 1) + labs(title = "(b)") + 
            labs(x = "Scaled rank of species", y = "Log (abundance)") + 
            scale_color_discrete("Provenance") +
            theme(strip.text = element_text(size = 12),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 10), 
                  legend.title = element_text(size = 12),
                  legend.text = element_text(size = 12)), 
          nrow = 2, common.legend = T, legend = "bottom"
)

# the top 3 species regarding abundance
subset(tree_rank_df, Rank <= 3)
subset(shrub_rank_df, Rank <= 3)
# calculate the EQ evenness index and plot
community_structure(tree_rank_df, time.var = "Land_use_type", 
                               abundance.var = "Abundance", metric = "EQ") %>% 
  arrange(EQ)
community_structure(shrub_rank_df, time.var = "Land_use_type", 
                                abundance.var = "Abundance", metric = "EQ") %>%
  arrange(EQ)
rm(tree_rank_ori, tree_rank_list, tree_rank_df, 
   shrub_rank_ori, shrub_rank_list, shrub_rank_df)



## Non-metric multidimensional scaling analysis
set.seed(1234)

# nMDS calculation for tree
tree_mds_selected <- subset(tree_diversity, Density > 3)
tree_mds_meta <- tree_mds_selected %>% 
  select(2:(number_tree_species+1)) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
tree_mds_meta$stress
stressplot(tree_mds_meta)
tree_mds_selected <- cbind(tree_mds_selected, tree_mds_meta$points)

# nMDS calculation for shrub
shrub_mds_selected <- shrub_diversity %>% filter(Density > 5)
shrub_mds_meta <- shrub_mds_selected %>% 
  select(2:(number_shrub_species+1)) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
shrub_mds_meta$stress
stressplot(shrub_mds_meta)
shrub_mds_selected <- cbind(shrub_mds_selected, shrub_mds_meta$points)

# ANOSIM of trees and shrubs as labels for the nMDS plots
tree_anosim <- 
  anosim(tree_mds_selected[2:(number_tree_species+1)], tree_mds_selected$Land_use_type)
shrub_anosim <- 
  anosim(shrub_mds_selected[2:(number_shrub_species+1)], shrub_mds_selected$Land_use_type)

# get hull for nMDS plots of trees and shrubs
fun_find_hull <- function(x) {x[chull(x$MDS1, x$MDS2), ]}
tree_hulls <- ddply(tree_mds_selected, "Land_use_type", fun_find_hull)
shrub_hulls <- ddply(shrub_mds_selected, "Land_use_type", fun_find_hull)

# nMDS plots for trees and shrubs by land use types
ggarrange(ggplot(tree_mds_selected, aes(MDS1, MDS2, color = Land_use_type)) + 
            geom_point(size=3) +
            labs(title = "Tree", 
                 subtitle = paste("stress=", round(tree_mds_meta$stress, digits = 3),
                                  ", R=", round(tree_anosim$statistic, digits = 3),
                                  ", p=", round(tree_anosim$signif, digits = 3),sep = "")) +
            geom_polygon(data=tree_hulls, alpha = 0, aes(fill=Land_use_type), size=1) +
            theme(axis.text = element_text(size = 12), 
                  legend.title = element_text(size = 15), 
                  legend.text = element_text(size = 12)) +
            theme_bw(),
          ggplot(shrub_mds_selected, aes(MDS1, MDS2, color = Land_use_type)) + 
            geom_point(size=3) +
            labs(title = "Shrub", 
                 subtitle = paste("stress=", round(shrub_mds_meta$stress, digits = 3),
                                  ", R=", round(shrub_anosim$statistic, digits = 3),
                                  ", p=", round(shrub_anosim$signif, digits = 3), sep = "")) +
            geom_polygon(data=shrub_hulls, alpha = 0, aes(fill=Land_use_type), size=1) +
            theme(axis.text = element_text(size = 12), 
                  legend.title = element_text(size = 15), 
                  legend.text = element_text(size = 12)) +
            theme_bw(),
          common.legend = T, legend = "right"
)

# pairwise result of ANOSIM of trees by Land_use_type
tree_pair_anosim_list <- vector("list",3)
tree_pair_anosim_list[[1]] <- c(combn(levels(factor(tree_mds_selected$Land_use_type)),2)[1,])
tree_pair_anosim_list[[2]] <- c(combn(levels(factor(tree_mds_selected$Land_use_type)),2)[2,])
set.seed(1234)
for (i in 1:length(tree_pair_anosim_list[[1]])) {
  tree_mds_selected_sub <- subset(tree_mds_selected, Land_use_type == tree_pair_anosim_list[[1]][i] |
                                    Land_use_type == tree_pair_anosim_list[[2]][i])
  tree_pair_anosim_list[[3]] <- c(tree_pair_anosim_list[[3]], 
                                  anosim(tree_mds_selected_sub[2:(number_tree_species+1)], 
                                         tree_mds_selected_sub$Land_use_type)$signif)
}
tree_pair_anosim_df <- data.frame(comp_1 = tree_pair_anosim_list[[1]], 
                                  comp_2 = tree_pair_anosim_list[[2]], 
                                  p = tree_pair_anosim_list[[3]]) %>% subset(p < 0.05) %>% print()

# pairwise result of ANOSIM of shrubs by Land_use_type
shrub_pair_anosim_list <- vector("list",3)
shrub_pair_anosim_list[[1]] <- c(combn(levels(factor(shrub_mds_selected$Land_use_type)),2)[1,])
shrub_pair_anosim_list[[2]] <- c(combn(levels(factor(shrub_mds_selected$Land_use_type)),2)[2,])
set.seed(1234)
for (i in 1:length(shrub_pair_anosim_list[[1]])) {
  shrub_mds_selected_sub <- subset(shrub_mds_selected, Land_use_type == shrub_pair_anosim_list[[1]][i] |
                                     Land_use_type == shrub_pair_anosim_list[[2]][i])
  shrub_pair_anosim_list[[3]] <- c(shrub_pair_anosim_list[[3]], 
                                   anosim(shrub_mds_selected_sub[2:(number_shrub_species+1)], 
                                          shrub_mds_selected_sub$Land_use_type)$signif)
}
shrub_pair_anosim_df <- data.frame(comp_1 = shrub_pair_anosim_list[[1]], 
                                   comp_2 = shrub_pair_anosim_list[[2]], 
                                   p = shrub_pair_anosim_list[[3]]) %>% 
  subset(p < 0.05) %>% print()




## cor among the indexes
chart.Correlation(subset(tree_diversity, select = c("Density", "Richness", "Shannon", "Simpson", "Evenness")))
chart.Correlation(subset(shrub_diversity, select = c("Density", "Richness", "Shannon", "Simpson", "Evenness")))



## Kruskal-Wallis test & boxplot for trees 
# p-value matrix for tree 
set.seed(1234)
boxplot_list_index_attr <- vector("list", 2)
#
{
  pvalue_list <- vector("list", 2)
  for (i in c("Density", "Richness", "Shannon", "Evenness")) {
    pvalue_list[[1]] <- c(pvalue_list[[1]], i)
    pvalue_list[[2]] <- c(pvalue_list[[2]], 
                          round(kruskal.test(tree_diversity[, i] ~ tree_diversity$Land_use_type)$p.value,digits = 3))
  }
  pvalue <- data.frame(index = pvalue_list[[1]],
                       pvalue = pvalue_list[[2]]) %>% 
    mutate(index = factor(index, levels = c("Density", "Richness", "Shannon", "Evenness")))
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
  scale_y_continuous(expand = expansion(mult = c(0.05,0.3))) +
  geom_text(data = pvalue, aes(x =Inf, y = Inf, label = label), 
            size=3.5, hjust = 1.05, vjust = 1.5) +
  theme(axis.text.x = element_text(angle = 90, size = 12)) + 
  labs(title = "(a)", x = NULL, y = NULL)
#
# for shrub
{
  pvalue_list <- vector("list", 2)
  for (i in c("Density", "Richness", "Shannon", "Evenness")) {
    pvalue_list[[1]] <- c(pvalue_list[[1]], i)
    pvalue_list[[2]] <- c(pvalue_list[[2]], 
                          round(kruskal.test(shrub_diversity[, i] ~ shrub_diversity$Land_use_type)$p.value,digits = 3))
  }
  pvalue <- data.frame(index = pvalue_list[[1]],
                       pvalue = pvalue_list[[2]]) %>% 
    mutate(index = factor(index, levels = c("Density", "Richness", "Shannon", "Evenness")))
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
  scale_y_continuous(expand = expansion(mult = c(0.05,0.3))) +
  geom_text(data = pvalue, aes(x =Inf, y = Inf, label = label), 
            size=3.5, hjust = 1.05, vjust = 1.5) +
  theme(axis.text.x = element_text(angle = 90, size = 12)) + 
  labs(title = "(b)", x = NULL, y = NULL)
Rmisc::multiplot(plotlist = boxplot_list_index_attr, cols = 2)
# delete the vars
rm(pvalue, pvalue_list)



## pairwise dunn test of diversity ~ landuse class
pairwise_list <- vector("list", 4)
# list of pairwise test of diversity of tree
for (i in c("Density", "Richness", "Shannon", "Evenness")) {
  pairwise_list[[1]] <- c(pairwise_list[[1]], 
                          rep("tree", 15))
  pairwise_list[[2]] <- c(pairwise_list[[2]], 
                          rep(i, 15))
  pairwise_list[[3]] <- c(pairwise_list[[3]], 
                          dunn.test(tree_diversity[, i], tree_diversity$Land_use_type)$comparisons)
  pairwise_list[[4]] <- c(pairwise_list[[4]], 
                          dunn.test(tree_diversity[, i], tree_diversity$Land_use_type)$P.adjusted)
}
# list of pairwise test of diversity and attrs of shrub
for (i in c("Density", "Richness", "Shannon", "Evenness")) {
  pairwise_list[[1]] <- c(pairwise_list[[1]], 
                          rep("shrub", 15))
  pairwise_list[[2]] <- c(pairwise_list[[2]], 
                          rep(i, 15))
  pairwise_list[[3]] <- c(pairwise_list[[3]], 
                          dunn.test(shrub_diversity[, i], shrub_diversity$Land_use_type)$comparisons)
  pairwise_list[[4]] <- c(pairwise_list[[4]], 
                          dunn.test(shrub_diversity[, i], shrub_diversity$Land_use_type)$P.adjusted)
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
  mutate(index = factor(index, levels = c("Density", "Density", "Richness", "Shannon", "Evenness")), 
         comparison_1 = factor(comparison_1, levels = Land_use_type_faclev), 
         comparison_2 = factor(comparison_2, levels = Land_use_type_faclev))
# plot the pairwise test results
ggarrange(ggplot(subset(pairwise_df, taxa == "tree"), 
                 aes(comparison_1, comparison_2, fill = p))+
            geom_tile() + geom_text(aes(label = round(p*100)), size = 2.5) +
            scale_fill_gradient2(high = "blue", low = "red", 
                                 midpoint = 0.05, limits = c(0, 0.05)) + 
            theme(axis.text.x = element_text(angle = 90)) + 
            xlab(NULL) + ylab(NULL) + guides(fill = FALSE) + 
            facet_wrap(~ index, scales = "free") +
            labs(title = "Tree"),
          ggplot(subset(pairwise_df, taxa == "shrub"), 
                 aes(comparison_1, comparison_2, fill = p))+
            geom_tile() + geom_text(aes(label = round(p*100)), size = 2.5) +
            scale_fill_gradient2(high = "blue", low = "red", 
                                 midpoint = 0.05, limits = c(0, 0.05)) + 
            theme(axis.text.x = element_text(angle = 90)) + 
            xlab(NULL) + ylab(NULL) + guides(fill = FALSE) + 
            facet_wrap(~ index, scales = "free") + 
            labs(title = "Shrub")
)




