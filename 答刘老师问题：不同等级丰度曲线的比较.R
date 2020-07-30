tree_rankabun_list <- vector("list", 6)
for (i in c("Com", "Com neigh", "R low", "R high", "R resi", "Ind")) {
  tree_diversity_landuse <- subset(tree_diversity, Landuse_class == i)[2:(number_tree_species+1)]
  tree_diversity_landuse <- subset(tree_diversity_landuse, select = (colSums(tree_diversity_landuse) != 0))
  tree_rankabun_ori <- as.data.frame(rankabundance(tree_diversity_landuse))
  tree_rankabun_list[[1]] <- c(tree_rankabun_list[[1]], rownames(tree_rankabun_ori))
  tree_rankabun_list[[2]] <- c(tree_rankabun_list[[2]], tree_rankabun_ori$rank)
  tree_rankabun_list[[3]] <- c(tree_rankabun_list[[3]], tree_rankabun_ori$rankfreq)
  tree_rankabun_list[[4]] <- c(tree_rankabun_list[[4]], tree_rankabun_ori$abundance)
  tree_rankabun_list[[5]] <- c(tree_rankabun_list[[5]], tree_rankabun_ori$proportion)
  tree_rankabun_list[[6]] <- c(tree_rankabun_list[[6]], rep(i, nrow(tree_rankabun_ori)))
}
tree_rankabun_df <- data.frame(
  Species_CN = tree_rankabun_list[[1]], 
  rank = tree_rankabun_list[[2]],
  rankfreq = tree_rankabun_list[[3]], 
  abundance = tree_rankabun_list[[4]], 
  proportion = tree_rankabun_list[[5]], 
  Landuse_class = tree_rankabun_list[[6]]
) %>% 
  left_join(all_plant_info[, c("Species_CN", "Nt_ex")], by = "Species_CN") %>% 
  mutate(Landuse_class = factor(Landuse_class, levels = Landuse_class_faclev))
# plot
ggarrange(ggplot(tree_rankabun_df, aes(rank, abundance, label = Species_CN)) + 
            geom_line() + 
            geom_point(alpha = 0.3, size = 2) + 
            facet_wrap(~Landuse_class, nrow = 1) + labs(title = "(rank~abundance)"),
          ggplot(tree_rankabun_df, aes(rank, proportion, label = Species_CN)) + 
            geom_line() + 
            geom_point(alpha = 0.3, size = 2) + 
            facet_wrap(~Landuse_class, nrow = 1) + labs(title = "(rank~proportion)"),
          ggplot(tree_rankabun_df, aes(rankfreq, abundance, label = Species_CN)) + 
            geom_line() + 
            geom_point(alpha = 0.3, size = 2) + 
            facet_wrap(~Landuse_class, nrow = 1) + labs(title = "(rankfreq~abundance)"),
          ggplot(tree_rankabun_df, aes(rankfreq, proportion, label = Species_CN)) + 
            geom_line() + 
            geom_point(alpha = 0.3, size = 2) + 
            facet_wrap(~Landuse_class, nrow = 1) + labs(title = "(rankfreq~proportion)"),
          nrow = 4
)
# log Y scale
ggarrange(ggplot(tree_rankabun_df, aes(rank, log(abundance), label = Species_CN)) + 
            geom_line() + 
            geom_point(alpha = 0.3, size = 2) + 
            facet_wrap(~Landuse_class, nrow = 1) + labs(title = "(rank~log abundance)"),
          ggplot(tree_rankabun_df, aes(rank, log(proportion), label = Species_CN)) + 
            geom_line() + 
            geom_point(alpha = 0.3, size = 2) + 
            facet_wrap(~Landuse_class, nrow = 1) + labs(title = "(rank~log proportion)"),
          ggplot(tree_rankabun_df, aes(rankfreq, log(abundance), label = Species_CN)) + 
            geom_line() + 
            geom_point(alpha = 0.3, size = 2) + 
            facet_wrap(~Landuse_class, nrow = 1) + labs(title = "(rankfreq~log abundance)"),
          ggplot(tree_rankabun_df, aes(rankfreq, log(proportion), label = Species_CN)) + 
            geom_line() + 
            geom_point(alpha = 0.3, size = 2) + 
            facet_wrap(~Landuse_class, nrow = 1) + labs(title = "(rankfreq~log proportion)"),
          nrow = 4
)
