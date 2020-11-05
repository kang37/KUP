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
#some default parameters
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
rm(tree_diversity_perc_planted, tree_diversity_perc_nonpot, tree_diversity_perc_private,
   tree_diversity_perc_nonstreet, tree_diversity_perc_native)

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
rm(shrub_diversity_perc_planted, shrub_diversity_perc_nonpot, shrub_diversity_perc_private,
   shrub_diversity_perc_nonstreet, shrub_diversity_perc_native)

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

# top species families of trees and shrubs by abundance
tree_data %>% 
  group_by(Family) %>% 
  dplyr::summarise(Number = sum(Stem), Prop = Number/sum(tree_data$Stem)) %>% 
  arrange(desc(Prop))
shrub_data %>% 
  group_by(Family) %>% 
  dplyr::summarise(SArea = sum(Area), Prop = SArea/sum(shrub_data$Area)) %>% 
  arrange(desc(Prop))


## attributes of trees and shrubs
# the number of exotic vs. native species
table(all_plant_info$Nt_ex)/nrow(all_plant_info)

# the attributes of trees and shrubs
for (i in c("Pla_spo", "Pub_pri", "Nt_ex")) {
  print(tapply(tree_data$Stem, tree_data[,i], sum)/sum(tree_data$Stem), digits = 3)
  cat("\n")
}
for (i in c("Pla_spo", "Pub_pri", "Nt_ex")) {
  print(tapply(shrub_data$Area, shrub_data[,i], sum)/sum(shrub_data$Area), digits = 3)
  cat("\n")
}


## Species accumulation curve 
fun_accum <- function(x, y, z) {
  ori <- specaccum(x[, 2:(y+1)])
  data.frame("number_of_plot" = ori$sites,
             "number_of_species" = ori$richness, 
             "tree_shrub" = z)
}
tree_accum <- ddply(tree_diversity, .(Land_use_type), 
                    y = number_tree_species, z = "tree", fun_accum)
shrub_accum <- ddply(shrub_diversity, .(Land_use_type), 
                     y = number_shrub_species, z = "shrub", fun_accum)
accum_df <- rbind(tree_accum, shrub_accum) %>% 
  mutate(Land_use_type = factor(Land_use_type, levels = Land_use_type_faclev), 
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
rm(tree_accum, shrub_accum, fun_accum)


## rank abundance plot
fun_rank <- function(x, y) {
  as.data.frame(x[2:(y+1)]) %>%
    subset(select = colSums(.) != 0) %>%
    rankabundance() %>% 
    as.data.frame() %>%
    mutate(Species_LT = rownames(.))
}
tree_rank_df <- 
  ddply(tree_diversity, "Land_use_type", y = number_tree_species, fun_rank) %>% 
  left_join(select(all_plant_info, c("Species_LT", "Nt_ex")), by = "Species_LT")
shrub_rank_df <- 
  ddply(shrub_diversity, "Land_use_type", y = number_shrub_species, fun_rank) %>% 
  left_join(select(all_plant_info, c("Species_LT", "Nt_ex")), by = "Species_LT")
ggarrange(ggplot(tree_rank_df, aes(rankfreq, log(abundance))) + 
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
          ggplot(shrub_rank_df, aes(rankfreq, log(abundance))) + 
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
subset(tree_rank_df[,c("rank", "Species_LT", "Nt_ex")], rank <= 3)
subset(shrub_rank_df[,c("rank", "Species_LT", "Nt_ex")], rank <= 3)
# calculate the EQ evenness index and plot
community_structure(
  tree_rank_df, time.var = "Land_use_type", abundance.var = "abundance", metric = "EQ") %>% 
  arrange(EQ)
community_structure(
  shrub_rank_df, time.var = "Land_use_type", abundance.var = "abundance", metric = "EQ") %>%
  arrange(EQ)
rm(tree_rank_df, shrub_rank_df, fun_rank)


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
fun_nmds_plot <- function(mds_selected, hull, plot_title, mds_meta, anosim) {
  ggplot(mds_selected, aes(MDS1, MDS2, color = Land_use_type)) + 
    geom_point(size=3) +
    geom_polygon(data = hull, alpha = 0, aes(fill=Land_use_type), size=1) +
    labs(title = plot_title, 
         subtitle = paste("stress=", round(mds_meta$stress, digits = 3),
                          ", R=", round(anosim$statistic, digits = 3),
                          ", p=", round(anosim$signif, digits = 3),sep = "")) +
    theme(axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 12)) +
    theme_bw()
}
ggarrange(
  fun_nmds_plot(tree_mds_selected, tree_hulls, "tree", tree_mds_meta, tree_anosim), 
  fun_nmds_plot(shrub_mds_selected, shrub_hulls, "shrub", shrub_mds_meta, shrub_anosim),
  common.legend = T, legend = "right"
)

# pairwise result of ANOSIM of trees by Land_use_type
anosim_pairs <- combn(levels(factor(tree_mds_selected$Land_use_type)),2)
pair_anosim_list <- vector("list",4)
pair_anosim_list[[1]] <- rep(c("tree", "shrub"), each = ncol(anosim_pairs))
pair_anosim_list[[2]] <- rep(anosim_pairs[1,], 2)
pair_anosim_list[[3]] <- rep(anosim_pairs[2,], 2)
for (i in 1:ncol(anosim_pairs)) {
  set.seed(1234)
  tree_mds_selected_sub <- subset(
    tree_mds_selected, Land_use_type == pair_anosim_list[[2]][i] | 
      Land_use_type == pair_anosim_list[[3]][i])
  pair_anosim_list[[4]] <- c(pair_anosim_list[[4]], 
                             anosim(tree_mds_selected_sub[2:(number_tree_species+1)], 
                                    tree_mds_selected_sub$Land_use_type)$signif)
}
for (i in 1:ncol(anosim_pairs)) {
  set.seed(1234)
  shrub_mds_selected_sub <- subset(
    shrub_mds_selected, Land_use_type == pair_anosim_list[[2]][i] | 
      Land_use_type == pair_anosim_list[[3]][i])
  pair_anosim_list[[4]] <- c(pair_anosim_list[[4]], 
                             anosim(shrub_mds_selected_sub[2:(number_shrub_species+1)], 
                                    shrub_mds_selected_sub$Land_use_type)$signif)
}
data.frame(Tree_shrub = pair_anosim_list[[1]],
           Comp_1 = pair_anosim_list[[2]], 
           Comp_2 = pair_anosim_list[[3]], 
           p = pair_anosim_list[[4]]) %>% 
  subset(p < 0.05)
rm(tree_anosim, tree_hulls, tree_mds_meta, tree_mds_selected, tree_mds_selected_sub,
   shrub_anosim, shrub_hulls, shrub_mds_meta, shrub_mds_selected, shrub_mds_selected_sub, 
   anosim_pairs, pair_anosim_list, fun_find_hull)



## cor among the indexes
chart.Correlation(subset(tree_diversity, 
                         select = c("Density", "Richness", "Shannon", "Simpson", "Evenness")))
chart.Correlation(subset(shrub_diversity, 
                         select = c("Density", "Richness", "Shannon", "Simpson", "Evenness")))



## Kruskal-Wallis test & box plot for trees 
# tree diversity longer and shrub diversity longer data set
fun_cons_long <- function(x) {
  subset(x, select = c("Density", "Richness", "Shannon", "Evenness","Land_use_type")) %>% 
    pivot_longer(cols = c("Density", "Richness", "Shannon", "Evenness"), 
                 names_to = "Index", values_to = "Index_value") %>% 
    mutate(Index = factor(Index, levels = c("Density", "Richness", "Shannon", "Evenness")), 
           Land_use_type = factor(Land_use_type, levels = c(Land_use_type_faclev)), 
           Attr = c("Land use type")) %>% 
    na.omit()
}
tree_diversity_long <- fun_cons_long(tree_diversity)
shrub_diversity_long <- fun_cons_long(shrub_diversity)

# get p-values for box plots
fun_get_pvalue <- function(x) {
  y <- data.frame(Index = c("Density", "Richness", "Shannon", "Evenness"), 
             Pvalue = NA, Label = NA) %>% 
    mutate(Index = factor(Index, levels = c("Density", "Richness", "Shannon", "Evenness")))
  j <- 0
  for (i in c("Density", "Richness", "Shannon", "Evenness")) {
    j <- j+1
    y$Pvalue[j] <- round(kruskal.test(
      as.data.frame(x)[, i] ~ x$Land_use_type)$p.value,digits = 3)
  }
  y$Label <- case_when(
    y$Pvalue >= 0.05 ~ 
      paste("p=", sprintf("%.3f",y$Pvalue), sep = ""), 
    y$Pvalue < 0.05 & y$Pvalue >= 0.01 ~ 
      paste("p=", sprintf("%.3f",y$Pvalue), "*", sep = ""), 
    y$Pvalue < 0.01 & y$Pvalue >= 0.001 ~ 
      paste("p=", sprintf("%.3f",y$Pvalue), "**", sep = ""),
    y$Pvalue < 0.001 ~ 
      paste("p=", sprintf("%.3f",y$Pvalue), "***", sep = "")
  )
  data.frame(y)
}
tree_box_pvalue <- fun_get_pvalue(tree_diversity)
shrub_box_pvalue <- fun_get_pvalue(shrub_diversity)

# get box plots
fun_box_plot <- function(x, y, z) {
  ggplot(x, aes(Land_use_type, Index_value)) + 
    geom_boxplot() + 
    facet_grid(Index ~ Attr, scales = "free", space = "free_x", switch = "both") + 
    scale_y_continuous(expand = expansion(mult = c(0.05,0.3))) +
    geom_text(data = y, aes(x =Inf, y = Inf, label = Label), 
              size=3.5, hjust = 1.05, vjust = 1.5) +
    theme(axis.text.x = element_text(angle = 90, size = 12)) + 
    labs(title = z, x = NULL, y = NULL)
}
ggarrange(fun_box_plot(tree_diversity_long, tree_box_pvalue, "(a)"), 
          fun_box_plot(shrub_diversity_long, shrub_box_pvalue, "(b)"))
rm(tree_box_pvalue, tree_diversity_long, 
   shrub_box_pvalue, shrub_diversity_long, 
   fun_box_plot, fun_cons_long, fun_get_pvalue)



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




