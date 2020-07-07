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



## attributes of the tree and shrub
# the graph of attributes of the plants
par(mfrow= c(2,2), cex.axis = 1.5)
j <- 0
for (i in c("Pla_spo", "Pub_pri")) {
  j <- j + 1
  barplot(tapply(subset(all_plant_data, Tree_shrub == "tree")[, "Stem"], 
                 subset(all_plant_data, Tree_shrub == "tree")[, i], 
                 sum), ylim = c(0, 1400))
  title(main = paste("(", letters[j], ")"), adj = 0)
}
for (i in c("Pla_spo", "Pub_pri")) {
  j <- j + 1
  barplot(tapply(subset(all_plant_data, Tree_shrub == "shrub")[, "Area"], 
                 subset(all_plant_data, Tree_shrub == "shrub")[, i], 
                 sum), ylim = c(0, 1400))
  title(main = paste("(", letters[j], ")"), adj = 0)
}
par(opar)



## Species accumulation curve
par(mfrow = c(1,2))
# Species accumulation curve for trees
plot(specaccum(subset(tree_native_diversity, Landuse_class == "Com")[,2:(number_tree_native_species+1)]), col = "red", lty = 1, 
     ci.lty = 0, xlim = c(0, 40), ylim = c(0, 70))
plot(specaccum(subset(tree_native_diversity, Landuse_class == "Com neigh")[,2:(number_tree_native_species+1)]),col = "red", lty = 2, 
     ci.lty = 0, add = T)
plot(specaccum(subset(tree_native_diversity, Landuse_class == "R low")[,2:(number_tree_native_species+1)]), col = "blue", lty = 1, 
     ci.lty = 0, add = T)
plot(specaccum(subset(tree_native_diversity, Landuse_class == "R high")[,2:(number_tree_native_species+1)]), col = "blue", lty = 2, 
     ci.lty = 0, add = T)
plot(specaccum(subset(tree_native_diversity, Landuse_class == "R resi")[,2:(number_tree_native_species+1)]), col = "blue", lty = 3, 
     ci.lty = 0, add = T)
plot(specaccum(subset(tree_native_diversity, Landuse_class == "Ind")[,2:(number_tree_native_species+1)]), col = "black", lty = 1, 
     ci.lty = 0, add = T)
# Species accumulation curve for shrubs
plot(specaccum(subset(shrub_native_diversity, Landuse_class == "Com")[,2:(number_shrub_native_species+1)]), col = "red", lty = 1, 
     ci.lty = 0, xlim = c(0, 40), ylim = c(0, 70))
plot(specaccum(subset(shrub_native_diversity, Landuse_class == "Com neigh")[,2:(number_shrub_native_species+1)]),col = "red", lty = 2, 
     ci.lty = 0, add = T)
plot(specaccum(subset(shrub_native_diversity, Landuse_class == "R low")[,2:(number_shrub_native_species+1)]), col = "blue", lty = 1, 
     ci.lty = 0, add = T)
plot(specaccum(subset(shrub_native_diversity, Landuse_class == "R high")[,2:(number_shrub_native_species+1)]), col = "blue", lty = 2, 
     ci.lty = 0, add = T)
plot(specaccum(subset(shrub_native_diversity, Landuse_class == "R resi")[,2:(number_shrub_native_species+1)]), col = "blue", lty = 3, 
     ci.lty = 0, add = T)
plot(specaccum(subset(shrub_native_diversity, Landuse_class == "Ind")[,2:(number_shrub_native_species+1)]), col = "black", lty = 1, 
     ci.lty = 0, add = T)
par(opar)



## rank ahundance plot
# the rank abundance plot by land use types of trees: doesn't omit site 279
tree_native_rankabun_list <- vector("list", 5)
for (i in c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")) {
  tree_native_rankabun_ori <- as.data.frame(rankabundance(subset(tree_native_diversity, Landuse_class == i)[2:(number_tree_native_species+1)]))
  tree_native_rankabun_list[[1]] <- c(tree_native_rankabun_list[[1]], rownames(tree_native_rankabun_ori))
  tree_native_rankabun_list[[2]] <- c(tree_native_rankabun_list[[2]], tree_native_rankabun_ori$rank)
  tree_native_rankabun_list[[3]] <- c(tree_native_rankabun_list[[3]], tree_native_rankabun_ori$abundance)
  tree_native_rankabun_list[[4]] <- c(tree_native_rankabun_list[[4]], tree_native_rankabun_ori$proportion)
  tree_native_rankabun_list[[5]] <- c(tree_native_rankabun_list[[5]], rep(i, nrow(tree_native_rankabun_ori)))
}
tree_native_rankabun_df <- data.frame(
  Species_CN = tree_native_rankabun_list[[1]], 
  rank = tree_native_rankabun_list[[2]], 
  abundance = tree_native_rankabun_list[[3]], 
  proportion = tree_native_rankabun_list[[4]], 
  Landuse_class = tree_native_rankabun_list[[5]]
)
tree_native_rankabun_df <- tree_native_rankabun_df[which(tree_native_rankabun_df$abundance != 0),]

# the rank abundance plot by land use types of shrubs
shrub_native_rankabun_list <- vector("list", 5)
for (i in c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")) {
  shrub_native_rankabun_ori <- as.data.frame(rankabundance(subset(shrub_native_diversity, Landuse_class == i)[2:(number_shrub_native_species+1)]))
  shrub_native_rankabun_list[[1]] <- c(shrub_native_rankabun_list[[1]], rownames(shrub_native_rankabun_ori))
  shrub_native_rankabun_list[[2]] <- c(shrub_native_rankabun_list[[2]], shrub_native_rankabun_ori$rank)
  shrub_native_rankabun_list[[3]] <- c(shrub_native_rankabun_list[[3]], shrub_native_rankabun_ori$abundance)
  shrub_native_rankabun_list[[4]] <- c(shrub_native_rankabun_list[[4]], shrub_native_rankabun_ori$proportion)
  shrub_native_rankabun_list[[5]] <- c(shrub_native_rankabun_list[[5]], rep(i, nrow(shrub_native_rankabun_ori)))
}
shrub_native_rankabun_df <- data.frame(
  Species_CN = shrub_native_rankabun_list[[1]], 
  rank = shrub_native_rankabun_list[[2]], 
  abundance = shrub_native_rankabun_list[[3]], 
  proportion = shrub_native_rankabun_list[[4]], 
  Landuse_class = shrub_native_rankabun_list[[5]]
) 
shrub_native_rankabun_df <- shrub_native_rankabun_df[which(shrub_native_rankabun_df$abundance != 0),]

# rearrange and plot
tree_native_rankabun_df <- tree_native_rankabun_df %>% 
  mutate(Landuse_class = factor(Landuse_class, levels = c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")))
shrub_native_rankabun_df <- shrub_native_rankabun_df %>% 
  mutate(Landuse_class = factor(Landuse_class, levels = c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")))
Rmisc::multiplot(plotlist = c(
  list(ggplot(tree_native_rankabun_df, aes(rank, proportion, label = Species_CN)) + 
         geom_line() + 
         geom_point(alpha = 0.3, size = 2) + 
         geom_text(aes(label = ifelse(rank<4, as.character(Species_CN), "")), hjust = -0.5, vjust = 0) + 
         facet_wrap(~Landuse_class, nrow = 1) + labs(title = "(a)")), 
  list(ggplot(shrub_native_rankabun_df, aes(rank, proportion, label = Species_CN)) + 
         geom_line() + 
         geom_point(alpha = 0.3, size = 2) + 
         geom_text(aes(label = ifelse(rank<4, as.character(Species_CN), "")), hjust = -0.5, vjust = 0) + 
         facet_wrap(~Landuse_class, nrow = 1) + labs(title = "(b)"))
))
rm(tree_native_rankabun_list, tree_native_rankabun_ori, tree_native_rankabun_df, 
   shrub_native_rankabun_list, shrub_native_rankabun_ori, shrub_native_rankabun_df)



## mds analysis
set.seed(1234)
#
# nmds calculation for tree
tree_native_mds_selected_ID <- tree_native_diversity$Plot_ID[!(tree_native_diversity$Plot_ID %in% c(172))]
tree_native_mds_selected <- tree_native_diversity %>% filter(Plot_ID %in% tree_native_mds_selected_ID)
tree_native_mds_metaMDS <- tree_native_mds_selected %>% 
  select(2:(number_tree_native_species+1)) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
tree_native_mds_metaMDS$stress
stressplot(tree_native_mds_metaMDS)
tree_native_mds_selected <- cbind(tree_native_mds_selected, tree_native_mds_metaMDS$points)
#
# nmds calculation for shrub
shrub_native_mds_selected_ID <- shrub_native_diversity$Plot_ID[!(shrub_native_diversity$Plot_ID %in% c(65, 274, 244, 164, 75, 172))]
shrub_native_mds_selected <- shrub_native_diversity %>% filter(Plot_ID %in% shrub_native_mds_selected_ID)
shrub_native_mds_metaMDS <- shrub_native_mds_selected %>% 
  select(2:(number_shrub_native_species+1)) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
shrub_native_mds_metaMDS$stress
stressplot(shrub_native_mds_metaMDS)
shrub_native_mds_selected <- cbind(shrub_native_mds_selected, shrub_native_mds_metaMDS$points)
#
# mds plot for trees and shrubs
# general Anosim of trees and shrubs
tree_landuseclass_anosim_result <- anosim(tree_native_mds_selected[2:(number_tree_native_species+1)], tree_native_mds_selected$Landuse_class)
shrub_landuseclass_anosim_result <- anosim(shrub_native_mds_selected[2:(number_shrub_native_species+1)], shrub_native_mds_selected$Landuse_class)
# get statistic results as labels for the mds plots
tree_landuseclass_mds_lab <- paste("stress=", round(tree_native_mds_metaMDS$stress, digits = 3), 
                                   ", R=", round(tree_landuseclass_anosim_result$statistic, digits = 3), 
                                   ", p=", round(tree_landuseclass_anosim_result$signif, digits = 3), 
                                   sep = "")
shrub_landuseclass_mds_lab <- paste("stress=", round(shrub_native_mds_metaMDS$stress, digits = 3), 
                                    ", R=", round(shrub_landuseclass_anosim_result$statistic, digits = 3), 
                                    ", p=", round(shrub_landuseclass_anosim_result$signif, digits = 3), 
                                    sep = "")
# mds plots for trees and shrubs by land use types and land ownership 
Rmisc::multiplot(plotlist = list(
  ggplot(tree_native_mds_selected, aes(MDS1, MDS2, color = Landuse_class)) + geom_point(alpha = 0.7) + 
    labs(title = "Tree - Land use type", subtitle = tree_landuseclass_mds_lab), 
  ggplot(shrub_native_mds_selected, aes(MDS1, MDS2, color = Landuse_class)) + geom_point(alpha = 0.7) + 
    labs(title = "Shrub - Land use type", subtitle = shrub_landuseclass_mds_lab)
), layout = matrix(1:2, nrow = 1, byrow = T))
#
# pairwise result of ANOSIM of trees by landuse_class
tree_pair_anosim_list <- vector("list",3)
tree_pair_anosim_list[[1]] <- c(combn(levels(factor(tree_native_mds_selected$Landuse_class)),2)[1,])
tree_pair_anosim_list[[2]] <- c(combn(levels(factor(tree_native_mds_selected$Landuse_class)),2)[2,])
set.seed(1234)
for (i in 1:length(tree_pair_anosim_list[[1]])) {
  tree_native_mds_selected_sub <- subset(tree_native_mds_selected, Landuse_class == tree_pair_anosim_list[[1]][i] |
                                           Landuse_class == tree_pair_anosim_list[[2]][i])
  tree_pair_anosim_list[[3]] <- c(tree_pair_anosim_list[[3]], 
                                  anosim(tree_native_mds_selected_sub[2:(number_tree_native_species+1)], 
                                         tree_native_mds_selected_sub$Landuse_class)$signif)
}
tree_pair_anosim_df <- data.frame(comp_1 = tree_pair_anosim_list[[1]], 
                                  comp_2 = tree_pair_anosim_list[[2]], 
                                  p = tree_pair_anosim_list[[3]]) %>% subset(p < 0.05) %>% print()
#
# pairwise result of ANOSIM of shrubs by landuse_class
shrub_pair_anosim_list <- vector("list",3)
shrub_pair_anosim_list[[1]] <- c(combn(levels(factor(shrub_native_mds_selected$Landuse_class)),2)[1,])
shrub_pair_anosim_list[[2]] <- c(combn(levels(factor(shrub_native_mds_selected$Landuse_class)),2)[2,])
set.seed(1234)
for (i in 1:length(shrub_pair_anosim_list[[1]])) {
  shrub_native_mds_selected_sub <- subset(shrub_native_mds_selected, Landuse_class == shrub_pair_anosim_list[[1]][i] |
                                            Landuse_class == shrub_pair_anosim_list[[2]][i])
  shrub_pair_anosim_list[[3]] <- c(shrub_pair_anosim_list[[3]], 
                                   anosim(shrub_native_mds_selected_sub[2:(number_shrub_native_species+1)], 
                                          shrub_native_mds_selected_sub$Landuse_class)$signif)
}
shrub_pair_anosim_df <- data.frame(comp_1 = shrub_pair_anosim_list[[1]], 
                                   comp_2 = shrub_pair_anosim_list[[2]], 
                                   p = shrub_pair_anosim_list[[3]]) %>% subset(p < 0.05) %>% print()

## cor among the indexes
chart.Correlation(subset(tree_native_diversity, select = c("Sum_stem", "Richness", "Shannon", "Simpson", "Evenness")))
chart.Correlation(subset(shrub_native_diversity, select = c("Sum_area", "Richness", "Shannon", "Simpson", "Evenness")))



## kruskal test & boxplot for trees
set.seed(1234) 
boxplot_list_index_attr <- vector("list", 2)
# for tree_native 
{
  pvalue_list <- vector("list", 2)
  for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
    pvalue_list[[1]] <- c(pvalue_list[[1]], i)
    pvalue_list[[2]] <- c(pvalue_list[[2]], 
                          round(kruskal.test(tree_native_diversity[, i] ~ tree_native_diversity$Landuse_class)$p.value,digits = 3))
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
boxplot_list_index_attr[[1]] <- tree_native_diversity_long %>% 
  na.omit() %>%
  ggplot(aes(attr_value, index_value)) + geom_boxplot() + 
  facet_grid(index ~ attr, scales = "free", space = "free_x", switch = "both") + 
  scale_y_continuous(expand = expand_scale(mult = c(0.05,0.3))) +
  geom_text(data = pvalue, aes(x =Inf, y = Inf, label = label), size=3.5, hjust = 1.05, vjust = 1.5) +
  theme(axis.text = element_text(angle = 90)) + 
  labs(title = "(a)", x = NULL, y = NULL)
#
# for shrub_native
{
  pvalue_list <- vector("list", 2)
  for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
    pvalue_list[[1]] <- c(pvalue_list[[1]], i)
    pvalue_list[[2]] <- c(pvalue_list[[2]], 
                          round(kruskal.test(shrub_native_diversity[, i] ~ shrub_native_diversity$Landuse_class)$p.value,digits = 3))
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
boxplot_list_index_attr[[2]] <- shrub_native_diversity_long %>% 
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
# list of pairwise test of diversity of tree_native
for (i in c("Sum_stem", "Richness", "Shannon", "Evenness")) {
  pairwise_list[[1]] <- c(pairwise_list[[1]], 
                          rep("tree", 15))
  pairwise_list[[2]] <- c(pairwise_list[[2]], 
                          rep(i, 15))
  pairwise_list[[3]] <- c(pairwise_list[[3]], 
                          dunn.test(tree_native_diversity[, i], tree_native_diversity$Landuse_class)$comparisons)
  pairwise_list[[4]] <- c(pairwise_list[[4]], 
                          dunn.test(tree_native_diversity[, i], tree_native_diversity$Landuse_class)$P.adjusted)
}
# list of pairwise test of diversity and attrs of shrub
for (i in c("Sum_area", "Richness", "Shannon", "Evenness")) {
  pairwise_list[[1]] <- c(pairwise_list[[1]], 
                          rep("shrub", 15))
  pairwise_list[[2]] <- c(pairwise_list[[2]], 
                          rep(i, 15))
  pairwise_list[[3]] <- c(pairwise_list[[3]], 
                          dunn.test(shrub_native_diversity[, i], shrub_native_diversity$Landuse_class)$comparisons)
  pairwise_list[[4]] <- c(pairwise_list[[4]], 
                          dunn.test(shrub_native_diversity[, i], shrub_native_diversity$Landuse_class)$P.adjusted)
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