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

# Get data ----
# some default factors
opar <- par(no.readonly = TRUE)
Land_use_type_faclev <- c("Com", "Com-neigh", "R-low", "R-high", "R-other", "Ind")
index_faclev <- c("Density", "Richness", "Shannon", "Simpson", "Evenness")

# information of all the plots
all_plot_info <- read.csv("In_plot_info.csv", stringsAsFactors = FALSE) %>% 
  mutate(Land_use_type = factor(Land_use_type, levels = Land_use_type_faclev))

# information of all the plants species: provenance and taxonomy 
all_plant_info <- read.csv("In_plant_info.csv", stringsAsFactors = FALSE)

# data of all the plants: taxonomy, attributes, abundance, etc.
all_plant_data <- read.csv("In_plant_data.csv", stringsAsFactors = FALSE) %>%
  left_join(all_plant_info, by = "Species_LT") %>%
  left_join(all_plot_info,by = "Plot_ID")

# data of trees and shrubs and diversity table of trees and shrubs 
tree_data <- subset(all_plant_data, Tree_shrub == "tree") %>% mutate(Area = NULL)
shrub_data <- subset(all_plant_data, Tree_shrub == "shrub") %>% mutate(Stem = NULL) %>% 
  mutate(Area = Area/10000)
fun_div <- function(x, y, z) {
  perc_planted <- x %>% group_by(Plot_ID) %>% 
    summarise(perc = sum(ifelse(Pla_spo == "planted", x[, z], 0)/sum(x[, z])))
  perc_nonpot <- x %>% group_by(Plot_ID) %>% 
    summarise(perc = sum(ifelse(Pot == "non_pot", x[, z], 0)/sum(x[, z])))
  perc_private <- x %>% group_by(Plot_ID) %>% 
    summarise(perc = sum(ifelse(Pub_pri == "private", x[, z], 0)/sum(x[, z])))
  perc_nonstreet <- x %>% group_by(Plot_ID) %>% 
    summarise(perc = sum(ifelse(Street == "non_street", x[, z], 0)/sum(x[, z])))
  perc_native <- x %>% group_by(Plot_ID) %>% 
    summarise(perc = sum(ifelse(Nt_ex == "native", x[, z], 0)/sum(x[, z])))
  subset(x, select = c("Plot_ID", "Species_LT", z)) %>%
    pivot_wider(names_from = Species_LT, values_from = {{y}}, 
                values_fn = sum, values_fill = 0) %>% 
    mutate(Density = rowSums(.[2:ncol(.)]), 
           Richness = apply(.[2:ncol(.)]>0, 1, sum),
           Shannon = diversity(.[2:ncol(.)], index = "shannon"), 
           Simpson = diversity(.[2:ncol(.)], index = "simpson"),
           Evenness = Shannon / log(Richness), 
           perc_planted = perc_planted$perc, 
           perc_nonpot = perc_nonpot$perc, 
           perc_private = perc_private$perc, 
           perc_nonstreet = perc_nonstreet$perc, 
           perc_native = perc_native$perc) %>% 
    left_join(all_plot_info, by = "Plot_ID") %>% 
    as.data.frame()
}
tree_diversity <- fun_div(tree_data, Stem, "Stem")
shrub_diversity <- fun_div(shrub_data, Area, "Area")
rm(fun_div)

# some other variables 
number_tree_species <- length(unique(tree_data$Species_LT))
number_shrub_species <- length(unique(shrub_data$Species_LT))


# City level ----
## Number of species ----
cat("\n", "total species:", length(unique(all_plant_data$Species_LT)), "\n", 
    "total genera:", length(unique(all_plant_data$Genus)), "\n", 
    "total families:", length(unique(all_plant_data$Family)), "\n", "\n")

# species accumulation curve 
# func to plot species accumulation curve and extrapolation 
# func parameters: x, raw data; y, size of extrapolation; z, extent of x axis; method, plot the curve separately by land use (method = "landuse") or not (method = "all")
fun_accum <- function(x, y, z, method) {
  # func to derive list of incidence data of certain land use
  # func para: xtemp, raw data; ytemp, name of land use
  fun_temp <- function(xtemp, ytemp) {
    xtemp %>% filter(Land_use_type %in% ytemp) %>% 
      group_by(Species_LT) %>% 
      summarise(n = length(unique(Plot_ID))) %>% 
      arrange(desc(n)) %>% 
      .$n %>% 
      c(length(unique(xtemp[which(xtemp[, "Land_use_type"] %in% ytemp), "Plot_ID"])), .)
  }
  if (method == "landuse") {
    incidence <- list(
      fun_temp(x, "Com"), 
      fun_temp(x, "Com-neigh"),
      fun_temp(x, "Ind"), 
      fun_temp(x, "R-other"), 
      fun_temp(x, "R-high"), 
      fun_temp(x, "R-low")
    )
    names(incidence) <- c("Com", "Com-neigh", "Ind", "R-other", "R-high", "R-low")
  } else {
    incidence <- list(fun_temp(x, Land_use_type_faclev))
    names(incidence) <- "all"
  }
  
  accum <- iNEXT(incidence, q = 0, datatype = "incidence_freq", 
                 size = seq(1:y), se = FALSE)
  if (method == "landuse") {
    accum$iNextEst$Com$land_use <- "Com"
    accum$iNextEst$"Com-neigh"$land_use <- "Com-neigh"
    accum$iNextEst$Ind$land_use <- "Ind"
    accum$iNextEst$"R-other"$land_use <- "R-other"
    accum$iNextEst$"R-high"$land_use <- "R-high"
    accum$iNextEst$"R-low"$land_use <- "R-low"
    accum <- Reduce(rbind, accum$iNextEst)
  } else {
    accum$iNextEst$all$land_use <- "all"
    accum <- accum$iNextEst$all
  }
  accum$method[accum$method == "interpolated"] <- "observed"
  
  hline <- accum %>% group_by(land_use) %>% 
    summarise(asymtote = max(qD))
  
  if (method == "all") {
    plotdata <- ggplot(accum) + 
      geom_line(aes(t, qD, linetype = method), size = 1) +
      geom_hline(data = hline, aes(yintercept = asymtote), 
                 linetype = 2, size = 1) + 
      scale_linetype_discrete(limits = c("observed", "extrapolated"))
  } else {
    ggplot(accum[which(accum$method == "observed"), ]) + 
      geom_line(aes(t, qD, color = land_use, linetype = method), size = 1)+
      geom_hline(data = hline, aes(yintercept = asymtote, color = land_use), 
                 linetype = 2, size = 1)
  }
  plotdata + 
    labs(x = "Number of quadrats", y = "Nuber of species") + 
    scale_x_continuous(limits = c(0, z)) + 
    theme_bw()
}
fun_accum(all_plant_data, 600, 600, method = "all") + 
  theme(legend.position = "none")

## Top taxa ----
# top species families of all plants by species number
all_plant_info %>% group_by(Family) %>% 
  dplyr::summarise(Number = n(), Prop = n()/nrow(all_plant_info)) %>% 
  arrange(desc(Prop))

# abundance of trees and shrubs
cat("number of trees:", nrow(tree_data), "in", nrow(tree_diversity), "plot", "\n", 
    "area of shrubs:", sum(shrub_data$Area), "m2 in", nrow(shrub_diversity), "plot")

# top species families of trees and shrubs by abundance
# func to generate top species with abundance data
# func para: x, row data; y, level of classification; z, column of abundance data; k, calculation of total abundance; n, number of top "n"
fun_top <- function(x, y, z, k, n) {
  top_plant <- x %>%¡¡group_by(get(y)) %>% 
    dplyr::summarise(Number = sum(get(z)), Prop = Number/k) %>% 
    arrange(desc(Number)) %>% 
    head(n)
  names(top_plant)[1] = c(y)
  print(top_plant)
}
tree_top_species <- fun_top(tree_data, "Species_LT", "Stem", sum(tree_data$Stem), 10)
tree_top_family <- fun_top(tree_data, "Family", "Stem", sum(tree_data$Stem), 10)
shrub_top_species <- fun_top(shrub_data, "Species_LT", "Area", sum(shrub_data$Area), 10)
shrub_top_family <- fun_top(shrub_data, "Family", "Area", sum(shrub_data$Area), 10)
intersect(tree_top_family$Family, shrub_top_family$Family)

# func to test if top species are contained in top families
# func para: x, raw data of species; y, raw data of families 
fun_contain <- function(x, y) {
  merge_data <- merge(x, all_plant_info, by = "Species_LT")
  data.frame("Family" = merge_data$Family, 
             "Speices_LT" = merge_data$Species_LT, 
             "Contain" = merge_data$Family %in% y$Family
  )
}
fun_contain(tree_top_species, tree_top_family)
fun_contain(shrub_top_species, shrub_top_family)
rm(tree_top_species, tree_top_family, 
   shrub_top_species, shrub_top_family, 
   fun_top, fun_contain)

# attributes of trees and shrubs
# the number of exotic vs. native by number of species
table(all_plant_info$Nt_ex)/nrow(all_plant_info)

# the attributes of trees and shrubs
for (i in c("Pla_spo", "Pub_pri", "Nt_ex")) {
  print(tapply(tree_data$Stem, tree_data[,i], sum)/sum(tree_data$Stem), digits = 2)
  cat("\n")
}
for (i in c("Pla_spo", "Pub_pri", "Nt_ex")) {
  print(tapply(shrub_data$Area, shrub_data[,i], sum)/sum(shrub_data$Area), digits = 2)
  cat("\n")
}


## species accumulation curve
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
tapply(accum_df$number_of_species, accum_df[, c("Land_use_type", "tree_shrub")], max)
rm(tree_accum, shrub_accum, accum_df,fun_accum)


## rank abundance plot
fun_rank <- function(x, y) {
  x[2:(y+1)] %>%
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
fun_rank_plot <- function(x, title) {
  ggplot(x, aes(rankfreq, log(abundance))) + 
    geom_line() + 
    geom_point(aes(color = Nt_ex), alpha = 0.3, size = 2) + 
    facet_wrap(~Land_use_type, nrow = 1) + labs(title = title) + 
    labs(x = "Scaled rank of species", y = "Log (abundance)") +
    scale_color_discrete("Provenance") +
    theme(strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
}
ggarrange(fun_rank_plot(tree_rank_df, "(a)"),
          fun_rank_plot(shrub_rank_df, "(b)"), 
          nrow = 2, common.legend = TRUE, legend = "bottom")
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
rm(tree_rank_df, shrub_rank_df, fun_rank, fun_rank_plot)


## Non-metric multidimensional scaling analysis
set.seed(1234)
# nMDS calculation for tree
tree_mds_selected <- subset(tree_diversity, Density > 1)
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
         subtitle = paste("stress=", sprintf("%.3f",mds_meta$stress),
                          ", R=", sprintf("%.3f",anosim$statistic),
                          ", p=", sprintf("%.3f",anosim$signif),sep = "")) +
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
anosim_pairs <- combn(Land_use_type_faclev, 2)
fun_anosim_pairs <- function(x, y) {
  result <- NULL
  for (i in 1:ncol(anosim_pairs)) {
    set.seed(1234)
    mds_selected_sub <- subset(x, Land_use_type %in% c(anosim_pairs[1,i], anosim_pairs[2,i]))
    result <- c(result, anosim(mds_selected_sub[2:(y+1)], 
                               mds_selected_sub$Land_use_type)$signif)
  }
  result
}
data.frame("Comp_1" = anosim_pairs[1,], "Comp_2" = anosim_pairs[2,], 
  "p" = fun_anosim_pairs(tree_mds_selected, number_tree_species)
) %>% subset(p < 0.05)
data.frame("Comp_1" = anosim_pairs[1,], "Comp_2" = anosim_pairs[2,], 
  "p" = fun_anosim_pairs(shrub_mds_selected, number_shrub_species)
) %>% subset(p < 0.05)
rm(tree_anosim, tree_hulls, tree_mds_meta, tree_mds_selected, 
   shrub_anosim, shrub_hulls, shrub_mds_meta, shrub_mds_selected, 
   anosim_pairs, fun_nmds_plot, fun_find_hull, fun_anosim_pairs)


## plant occupancy of species for different land use types
fun_occup_rate <- function(x,y) {
  apply(x[,2:(y+1)], 2, function(k)sum(ifelse(k>0,1,0))/length(k))
}
fun_occup_df <- function(x){
  data.frame("Com" = names(head(sort(x[,"Com"],decreasing = TRUE),10)), 
             "Com-neigh" = names(head(sort(x[,"Com-neigh"],decreasing = TRUE),10)),
             "R-low" = names(head(sort(x[,"R-low"],decreasing = TRUE),10)),
             "R-high" = names(head(sort(x[,"R-high"],decreasing = TRUE),10)), 
             "R-other" = names(head(sort(x[,"R-other"],decreasing = TRUE),10)), 
             "Ind" = names(head(sort(x[,"Ind"],decreasing = TRUE),10)))
}

tree_occup <- ddply(
  tree_diversity, .(Land_use_type), y = number_tree_species, fun_occup_rate) %>% 
  .[,-1] %>% t() 
colnames(tree_occup) = Land_use_type_faclev
tree_occup_top <- fun_occup_df(tree_occup)
write.csv(tree_occup_top, "C:/Users/kangj/Documents/R/KUP/occup.csv", row.names = FALSE)

shrub_occup <- ddply(
  shrub_diversity, .(Land_use_type), y = number_shrub_species, fun_occup_rate) %>% 
  .[,-1] %>% t() 
colnames(shrub_occup) = Land_use_type_faclev
shrub_occup_top <- fun_occup_df(shrub_occup)
write.table(shrub_occup_top, "C:/Users/kangj/Documents/R/KUP/occup.csv", 
              row.names = FALSE, append = TRUE, sep = ",")

occup_top <- rbind(tree_occup_top, shrub_occup_top)

# shared species over land use types
# shared species for both trees and shrubs over land use types
Reduce(intersect, list(tree_occup_top[,1], tree_occup_top[,2], tree_occup_top[,3], 
                       tree_occup_top[,4], tree_occup_top[,5], tree_occup_top[,6], 
                       shrub_occup_top[,1], shrub_occup_top[,2], shrub_occup_top[,3], 
                       shrub_occup_top[,4], shrub_occup_top[,5], shrub_occup_top[,6]))

# shared species for trees over land use types
Reduce(intersect, list(tree_occup_top[,1], tree_occup_top[,2], tree_occup_top[,3], 
                       tree_occup_top[,4], tree_occup_top[,5], tree_occup_top[,6]))
# shared species for shrubs over land use types
Reduce(intersect, list(shrub_occup_top[,1], shrub_occup_top[,2], shrub_occup_top[,3], 
                       shrub_occup_top[,4], shrub_occup_top[,5], shrub_occup_top[,6]))

# shared species for trees between land use types
fun_share_prop <- function(occup_top_data, plot_title) {
  share_prop <- as.data.frame(matrix(numeric(0),ncol=3, nrow = 36))
  colnames(share_prop) <- c("land_use_1", "land_use_2", "prop")
  k <- 0
  for (i in c(1:6)) {
    for (j in c(1:6)) {
      k <- k+1
      share_prop$land_use_1[k] <- Land_use_type_faclev[i]
      share_prop$land_use_2[k] <- Land_use_type_faclev[j]
      share_prop$prop[k] <- (length(intersect(occup_top_data[,i],occup_top_data[,j])))/
        (length(union(occup_top_data[,i],occup_top_data[,j])))
    }
  }
  ggplot(share_prop, aes(land_use_1, land_use_2, fill = prop)) + 
    geom_tile() + geom_text(aes(label = round(prop*100))) +
    scale_fill_gradient2(high = "red", low = "blue", midpoint = 0.4, limits = c(0,0.8)) +
    labs(title = plot_title) + theme(axis.text.x = element_text(angle = 90))
}
ggarrange(fun_share_prop(tree_occup_top, "tree"), 
          fun_share_prop(shrub_occup_top, "shrub"), 
          fun_share_prop(occup_top, "common"), 
          nrow = 1, common.legend = TRUE, legend = "right")

# unique ubiquitous species in certain land use types
occup_top_long <- occup_top %>% pivot_longer(
  cols = c("Com", "Com.neigh", "R.low", "R.high", "R.other",  "Ind"),  
  names_to = "Land_use_type", values_to = "Species") %>% unique()
Species_unique <- occup_top_long %>% group_by(Species) %>% 
  dplyr::summarise(Number = n()) %>% arrange(Number) %>% 
  subset(Number == 1) %>% .$Species 
occup_top_long[which(occup_top_long$Species %in% Species_unique),] %>% arrange(Land_use_type)

rm(tree_occup, shrub_occup, tree_occup_top, shrub_occup_top, occup_top, occup_top_long, 
  Species_unique, fun_occup_rate,fun_occup_df, fun_share_prop)


## cor among the indexes
chart.Correlation(subset(tree_diversity, select = index_faclev))
chart.Correlation(subset(shrub_diversity, select = index_faclev))


## Kruskal-Wallis test & box plot for trees 
# tree diversity longer and shrub diversity longer data set
fun_cons_long <- function(x) {
  subset(x, select = c("Density", "Richness", "Shannon", "Evenness", "Land_use_type")) %>% 
    pivot_longer(cols = c("Density", "Richness", "Shannon", "Evenness"), 
                 names_to = "Index", values_to = "Index_value") %>% 
    mutate(Index = factor(Index, levels = c("Density", "Richness", "Shannon", "Evenness")), 
           Land_use_type = factor(Land_use_type, levels = Land_use_type_faclev), 
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
      x[, i] ~ x$Land_use_type)$p.value,digits = 3)
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


## pairwise dunn test of diversity ~ land use type
fun_dunn <- function(x, taxa, index) {
  dunn_result <- dunn.test(x[ , index], x$Land_use_type, table = FALSE, kw = FALSE)
  x <- data.frame(
    "taxa" = taxa, 
    "index" = index, 
    "comparison" = dunn_result$comparisons, 
    "p" = dunn_result$P.adjusted
  ) %>% 
    separate(comparison, into = c("comparison_1", "comparison_2"), sep = " - ")
}
dunn_df_1 <- rbind(fun_dunn(tree_diversity, "tree", "Density"), 
                   fun_dunn(tree_diversity, "tree", "Richness"),
                   fun_dunn(tree_diversity, "tree", "Shannon"),
                   fun_dunn(tree_diversity, "tree", "Evenness"),
                   fun_dunn(shrub_diversity, "shrub", "Density"), 
                   fun_dunn(shrub_diversity, "shrub", "Richness"),
                   fun_dunn(shrub_diversity, "shrub", "Shannon"),
                   fun_dunn(shrub_diversity, "shrub", "Evenness")
) 
dunn_df_2 <- dunn_df_1[, c("taxa", "index", "comparison_2", "comparison_1", "p")]
names(dunn_df_2) <- c("taxa", "index", "comparison_1", "comparison_2", "p")
dunn_df <- rbind(dunn_df_1, dunn_df_2)%>% 
  mutate(index = factor(index, levels = index_faclev), 
         comparison_1 = factor(comparison_1, Land_use_type_faclev), 
         comparison_2 = factor(comparison_2, Land_use_type_faclev))
# plot the pairwise test results
fun_dunn_plot <- function(x, title) {
  ggplot(x, aes(comparison_1, comparison_2, fill = p)) +
    geom_tile() + 
    geom_text(aes(label = round(p*100)), size = 2.5) +
    scale_fill_gradient2(high = "blue", low = "red", 
                         midpoint = 0.05, limits = c(0, 0.05)) + 
    theme(axis.text.x = element_text(angle = 90)) + 
    xlab(NULL) + ylab(NULL) + guides(fill = FALSE) + 
    facet_wrap(~ index, scales = "free", nrow = 1) +
    labs(title = title)
}
ggarrange(fun_dunn_plot(subset(dunn_df, taxa == "tree"), "Tree"), 
          fun_dunn_plot(subset(dunn_df, taxa == "shrub"), "Shrub"), 
          nrow = 2)
rm(dunn_df_1, dunn_df_2, dunn_df, 
   fun_dunn, fun_dunn_plot)

## data for discussion
# means of quadrat density and richness for trees
tree_diversity %>% group_by(Land_use_type) %>% 
  dplyr::summarise(Density = mean(Density), Richness = mean(Richness))

