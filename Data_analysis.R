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
library(iNEXT)

# Functions ----
fun_div <- function(x, y, z, k, m, method) {
  perc_planted <- x %>% group_by({{k}}) %>% 
    summarise(perc = sum(ifelse(Pla_spo == "planted", x[, z], 0)/sum(x[, z])))
  perc_nonpot <- x %>% group_by({{k}}) %>% 
    summarise(perc = sum(ifelse(Pot == "non_pot", x[, z], 0)/sum(x[, z])))
  perc_private <- x %>% group_by({{k}}) %>% 
    summarise(perc = sum(ifelse(Pub_pri == "private", x[, z], 0)/sum(x[, z])))
  perc_nonstreet <- x %>% group_by({{k}}) %>% 
    summarise(perc = sum(ifelse(Street == "non_street", x[, z], 0)/sum(x[, z])))
  perc_native <- x %>% group_by({{k}}) %>% 
    summarise(perc = sum(ifelse(Nt_ex == "native", x[, z], 0)/sum(x[, z])))
  
  newdf <- subset(x, select = c(m, "Species_LT", z)) %>% 
    pivot_wider(names_from = Species_LT, values_from = {{y}}, 
                values_fn = sum, values_fill = 0) %>% 
    mutate(Abundance = rowSums(.[2:ncol(.)]), 
           Richness = apply(.[2:ncol(.)]>0, 1, sum),
           Shannon = diversity(.[2:ncol(.)], index = "shannon"), 
           Simpson = diversity(.[2:ncol(.)], index = "simpson"),
           Evenness = Shannon / log(Richness), 
           perc_planted = perc_planted$perc, 
           perc_nonpot = perc_nonpot$perc, 
           perc_private = perc_private$perc, 
           perc_nonstreet = perc_nonstreet$perc, 
           perc_native = perc_native$perc)
  if (method == "land_use") {
    as.data.frame(newdf)
  } else {
    newdf %>% left_join(all_plot_info, by = m) %>% 
      as.data.frame()
  }
}

# func to plot species accumulation curve and extrapolation 
# func parameters: x, raw data; y, size of extrapolation; z, extent of x axis; method, plot the curve separately by land use (method = "land_use") or not (method = "city")
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
  if (method == "city") {
    incidence <- list(fun_temp(x, Land_use_type_faclev))
    names(incidence) <- "city"
  } else {
    incidence <- list(
      fun_temp(x, "Com"), 
      fun_temp(x, "ComNbr"),
      fun_temp(x, "Ind"), 
      fun_temp(x, "ResOther"), 
      fun_temp(x, "ResHigh"), 
      fun_temp(x, "ResLow")
    )
    names(incidence) <- c("Com", "ComNbr", "Ind", "ResOther", "ResHigh", "ResLow")
  }
  
  accum <- iNEXT(incidence, q = 0, datatype = "incidence_freq", 
                 size = seq(1:y), se = FALSE)
  if (method == "city") {
    accum$iNextEst$city$land_use <- "city"
    accum <- accum$iNextEst$city
  } else {
    accum$iNextEst$Com$land_use <- "Com"
    accum$iNextEst$"ComNbr"$land_use <- "ComNbr"
    accum$iNextEst$Ind$land_use <- "Ind"
    accum$iNextEst$"ResOther"$land_use <- "ResOther"
    accum$iNextEst$"ResHigh"$land_use <- "ResHigh"
    accum$iNextEst$"ResLow"$land_use <- "ResLow"
    accum <- Reduce(rbind, accum$iNextEst)
    accum$land_use <- factor(
      accum$land_use, 
      levels = c("Com", "ComNbr", "Ind", "ResOther", "ResHigh", "ResLow"))
  }
  accum$method[accum$method == "interpolated"] <- "observed"
  
  hline <- accum %>% group_by(land_use) %>% 
    summarise(asymtote = max(qD))
  print(arrange(hline, desc(asymtote)))
  
  if (method == "city") {
    plotdata <- ggplot(accum) + 
      geom_line(aes(t, qD, color = land_use, linetype = method), size = 1) +
      geom_hline(data = hline, aes(yintercept = asymtote, color = land_use), 
                 linetype = 2, size = 1) + 
      scale_linetype_discrete(limits = c("observed", "extrapolated"))
  } else {
    plotdata <- ggplot(accum[which(accum$method == "observed"), ]) + 
      geom_line(aes(t, qD, color = land_use), size = 1) +
      geom_hline(data = hline, aes(yintercept = asymtote, color = land_use), 
                 linetype = 2, size = 1)
  }
  plotdata + 
    labs(x = "Number of quadrats", y = "Nuber of species") + 
    scale_x_continuous(limits = c(0, z)) + 
    theme_bw()
}

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

# func to test if top species are contained in top families
# func para: x, raw data of species; y, raw data of families 
fun_contain <- function(x, y) {
  merge_data <- merge(x, all_plant_info, by = "Species_LT")
  data.frame("Family" = merge_data$Family, 
             "Speices_LT" = merge_data$Species_LT, 
             "Contain" = merge_data$Family %in% y$Family
  )
}

fun_rank_data <- function(x, y) {
  x[2:(y+1)] %>%
    subset(select = colSums(.) != 0) %>%
    rankabundance() %>% 
    as.data.frame() %>%
    mutate(Species_LT = rownames(.))
}

# func to plot rank-abundance curve
# func para: x, raw data; title, plot title; method, plot at city level (method = "city") or at land use level (method = "land_use)
fun_rank_plot <- function(x, title, method) {
  if (method == "city") {
    plotdata <- ggplot(x, aes(rankfreq, proportion)) + 
      geom_line() + labs(title = title)
  } else {
    plotdata <- ggplot(x, aes(rankfreq, proportion)) + 
      geom_line() + 
      facet_wrap(~Land_use_type, nrow = 1) + labs(title = title)
  }
  plotdata + 
    geom_point(aes(color = Nt_ex), alpha = 0.3, size = 2) + 
    labs(x = "Scaled rank of species", y = "Relative abundance") +
    scale_color_discrete("Provenance") +
    theme(strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
}

# Get data ----
# some default factors
Land_use_type_faclev <- c("Com", "ComNbr", "Ind", "ResOther", "ResHigh", "ResLow")
index_faclev <- c("Abundance", "Richness", "Shannon", "Simpson", "Evenness")

# information of all the plots
all_plot_info <- read.csv("In_plot_info.csv", stringsAsFactors = FALSE) %>% 
  mutate(Land_use_type = factor(Land_use_type, levels = Land_use_type_faclev))

# information of all the plants species: provenance and taxonomy 
all_plant_info <- read.csv("In_plant_info.csv", stringsAsFactors = FALSE)

# data of all the plants: taxonomy, attributes, abundance, etc.
all_plant_data <- read.csv("In_plant_data.csv", stringsAsFactors = FALSE) %>%
  left_join(all_plant_info, by = "Species_LT") %>%
  left_join(all_plot_info,by = "Plot_ID")
qua_plant_div <- all_plant_data %>% 
  mutate(presence = 1) %>% 
  subset(select = c("Plot_ID", "Species_LT", "presence")) %>% 
  unique() %>% 
  pivot_wider(names_from = Species_LT, values_from = presence, 
              values_fn = sum, values_fill = 0) %>% 
  mutate(Richness = apply(.[2:ncol(.)]>0, 1, sum)) %>% 
  left_join(all_plot_info, by = "Plot_ID") %>% 
  as.data.frame()

# data of trees and shrubs and diversity table of trees and shrubs 
tree_data <- subset(all_plant_data, Tree_shrub == "tree") %>% mutate(Area = NULL)
shrub_data <- subset(all_plant_data, Tree_shrub == "shrub") %>% mutate(Stem = NULL) %>% 
  mutate(Area = Area/10000)

lu_tree_div <- fun_div(tree_data, Stem, "Stem", 
                       Land_use_type, "Land_use_type", method = "land_use")
lu_shrub_div <- fun_div(shrub_data, Area, "Area", 
                        Land_use_type, "Land_use_type", method = "land_use")
qua_tree_div <- fun_div(tree_data, Stem, "Stem", 
                        Plot_ID, "Plot_ID", method = "quadrat")
qua_shrub_div <- fun_div(shrub_data, Area, "Area", 
                         Plot_ID, "Plot_ID", method = "quadrat")


# some other variables 
number_plant_species <- length(unique(all_plant_data$Species_LT))
number_tree_species <- length(unique(tree_data$Species_LT))
number_shrub_species <- length(unique(shrub_data$Species_LT))


# City level ----
## Number of species ----
cat("\n", "total species:", length(unique(all_plant_data$Species_LT)), "\n", 
    "total genera:", length(unique(all_plant_data$Genus)), "\n", 
    "total families:", length(unique(all_plant_data$Family)), "\n", "\n")

# species accumulation curve 
# extrapolation up to double the reference sample size
fun_accum(all_plant_data, 348, 350, method = "city") + 
  theme(legend.position = "none")

## Top taxa ----
# top species families of all plants by species number
all_plant_info %>% group_by(Family) %>% 
  dplyr::summarise(Num_spe = n(), Prop = n()/nrow(all_plant_info)) %>% 
  arrange(desc(Prop))

# abundance of trees and shrubs
cat("number of trees:", nrow(tree_data), "in", 
    nrow(qua_tree_div), "plot", "\n", 
    "area of shrubs:", sum(shrub_data$Area), "m2 in", 
    nrow(qua_shrub_div), "plot")

# top species families of trees and shrubs by abundance
tree_top_species <- fun_top(
  tree_data, "Species_LT", "Stem", sum(tree_data$Stem), 10)
tree_top_family <- fun_top(
  tree_data, "Family", "Stem", sum(tree_data$Stem), 10)
fun_contain(tree_top_species, tree_top_family)

shrub_top_species <- fun_top(
  shrub_data, "Species_LT", "Area", sum(shrub_data$Area), 10)
shrub_top_family <- fun_top(
  shrub_data, "Family", "Area", sum(shrub_data$Area), 10)
intersect(tree_top_family$Family, shrub_top_family$Family)
fun_contain(shrub_top_species, shrub_top_family)

rm(tree_top_species, tree_top_family, 
   shrub_top_species, shrub_top_family, 
   fun_top, fun_contain)

# attributes of trees and shrubs
# the number of exotic vs. native by number of species
table(all_plant_info$Nt_ex)/nrow(all_plant_info)

# the attributes of trees and shrubs
for (i in c("Pla_spo", "Pub_pri", "Nt_ex")) {
  print(tapply(tree_data$Stem, tree_data[,i], sum)/sum(tree_data$Stem), 
        digits = 2)
  cat("\n")
}
for (i in c("Pla_spo", "Pub_pri", "Nt_ex")) {
  print(tapply(shrub_data$Area, shrub_data[,i], sum)/sum(shrub_data$Area), 
        digits = 2)
  cat("\n")
}

## Distribution of species abundance ----
# rank abundance curves for trees and shrubs
city_tree_rank <- fun_rank_data(qua_tree_div, number_tree_species) %>% 
  left_join(select(all_plant_info, c("Species_LT", "Nt_ex")), by = "Species_LT")
city_shrub_rank <- fun_rank_data(qua_shrub_div, number_tree_species) %>% 
  left_join(select(all_plant_info, c("Species_LT", "Nt_ex")), by = "Species_LT")

ggarrange(plotlist = list(
  fun_rank_plot(city_tree_rank, title = "(a)", method = "city"), 
  fun_rank_plot(city_shrub_rank, title = "(b)", method = "city")
), nrow = 1, common.legend = TRUE)

# the top 3 species regarding abundance
subset(city_tree_rank[,c("rank", "Species_LT", "Nt_ex")], rank <= 3)
subset(city_shrub_rank[,c("rank", "Species_LT", "Nt_ex")], rank <= 3)

# EQ evenness index and plot
city_tree_rank %>% 
  community_structure(abundance.var = "abundance", metric = "EQ") %>% 
  arrange(EQ)
city_shrub_rank %>% 
  community_structure(abundance.var = "abundance", metric = "EQ") %>%
  arrange(EQ)

# Land use level ----
## Species accumulation curve ----
fun_accum(all_plant_data, 348, 50, method = "land_use") + 
  labs(color = "Land use") + 
  scale_color_manual(values = c("#FF0000", "#FF7800", "#DF73FF", 
                                "#BFBF30", "#6BE400", "#00733E"))

## Distribution of species abundance ----
lu_tree_rank <- 
  ddply(qua_tree_div, "Land_use_type", 
        y = number_tree_species, fun_rank_data) %>% 
  left_join(select(all_plant_info, c("Species_LT", "Nt_ex")), by = "Species_LT")
lu_shrub_rank <- 
  ddply(qua_shrub_div, "Land_use_type", 
        y = number_shrub_species, fun_rank_data) %>% 
  left_join(select(all_plant_info, c("Species_LT", "Nt_ex")), by = "Species_LT")
ggarrange(fun_rank_plot(lu_tree_rank, "(a)", method = "land_use"),
          fun_rank_plot(lu_shrub_rank, "(b)", method = "land_use"), 
          nrow = 2, common.legend = TRUE, legend = "bottom")
# the top 3 species regarding abundance
subset(lu_tree_rank[,c("rank", "Species_LT", "Nt_ex")], rank <= 3)
subset(lu_shrub_rank[,c("rank", "Species_LT", "Nt_ex")], rank <= 3)
# calculate the EQ evenness index and plot
community_structure(
  lu_tree_rank, time.var = "Land_use_type", 
  abundance.var = "abundance", metric = "EQ") %>% 
  arrange(desc(EQ))
community_structure(
  lu_shrub_rank, time.var = "Land_use_type", 
  abundance.var = "abundance", metric = "EQ") %>%
  arrange(desc(EQ))

## Species composition ----
# Bray-Curtis dissimilarity of pairs of land use
lu_dissim <- vegdist(lu_tree_div[2: 133]) %>% as.matrix() %>% round(digits = 2)
rownames(lu_dissim) <- Land_use_type_faclev
colnames(lu_dissim) <- Land_use_type_faclev
lu_dissim

# plant occupancy of species for different land use types
fun_occup_rate <- function(x,y) {
  apply(x[,2:(y+1)], 2, function(k)sum(ifelse(k>0,1,0))/length(k))
}
fun_occup_df <- function(x){
  data.frame("Com" = names(head(sort(x[,"Com"],decreasing = TRUE),10)), 
             "ComNbr" = names(head(sort(x[,"ComNbr"],decreasing = TRUE),10)),
             "ResLow" = names(head(sort(x[,"ResLow"],decreasing = TRUE),10)),
             "ResHigh" = names(head(sort(x[,"ResHigh"],decreasing = TRUE),10)), 
             "ResOther" = names(head(sort(x[,"ResOther"],decreasing = TRUE),10)), 
             "Ind" = names(head(sort(x[,"Ind"],decreasing = TRUE),10)))
}

plant_occup <- ddply(
  qua_plant_div, .(Land_use_type), y = number_plant_species, fun_occup_rate) %>% 
  .[,-1] %>% t() 
colnames(plant_occup) = Land_use_type_faclev
(plant_occup_top <- fun_occup_df(plant_occup))

# shared species over land use types
Reduce(intersect, list(plant_occup_top[,1], plant_occup_top[,2], 
                       plant_occup_top[,3], plant_occup_top[,4], 
                       plant_occup_top[,5], plant_occup_top[,6]))

# shared top species between land use types
fun_share_prop <- function(occup_top_data) {
  share_prop <- as.data.frame(matrix(numeric(0),ncol=3, nrow = 36))
  colnames(share_prop) <- c("land_use_1", "land_use_2", "prop")
  k <- 0
  for (i in c(1:6)) {
    for (j in c(1:6)) {
      k <- k+1
      share_prop$land_use_1[k] <- Land_use_type_faclev[i]
      share_prop$land_use_2[k] <- Land_use_type_faclev[j]
      share_prop$prop[k] <- 
        (length(intersect(occup_top_data[,i],occup_top_data[,j])))/
        (length(union(occup_top_data[,i],occup_top_data[,j])))
    }
  }
  ggplot(share_prop, aes(land_use_1, land_use_2, fill = prop)) + 
    geom_tile() + geom_text(aes(label = round(prop*100))) +
    scale_fill_gradient2(high = "red", low = "blue", 
                         midpoint = 0.4, limits = c(0,0.8)) +
    theme(axis.text.x = element_text(angle = 90))
}
fun_share_prop(plant_occup_top)

# unique ubiquitous species in certain land use types
plant_occup_top %>% pivot_longer(
  cols = all_of(Land_use_type_faclev),  
  names_to = "Land_use_type", values_to = "Species") %>%
  .[which(
  occup_top_long$Species %in% 
    # unique ubiquitous species - present once only
    names(table(as.character(as.matrix(plant_occup_top)))[
      table(as.character(as.matrix(plant_occup_top))) == 1])),] %>% 
  arrange(Land_use_type)

# Quadrat level ----
## Richness ~ land use for all plants ----
ggplot(qua_plant_div) + 
  geom_boxplot(aes(Land_use_type, Richness)) + 
  labs(x = "Land use zone", y = "Quadrat richness") + 
  geom_text(data = data.frame(
    Land_use_type = Land_use_type_faclev, 
    Label = c(rep("a", 4), "ab", "b")
  ), aes(x = Land_use_type, y = Inf, label = Label), vjust = 2) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme_bw()
kruskal.test(qua_plant_div$Richness ~ qua_plant_div$Land_use_type)
dunn.test(x = qua_plant_div$Richness, g = qua_plant_div$Land_use_type)

## Indexes ~ land use for trees and shrubs ----
# Kruskal-Wallis test & box plot for trees 
# tree diversity longer and shrub diversity longer data set
fun_cons_long <- function(x) {
  subset(x, select = c("Abundance", "Richness", "Evenness", "Land_use_type")) %>% 
    pivot_longer(cols = c("Abundance", "Richness", "Evenness"), 
                 names_to = "Index", values_to = "Index_value") %>% 
    mutate(Index = factor(Index, levels = c("Abundance", "Richness", "Evenness")), 
           Land_use_type = factor(Land_use_type, levels = Land_use_type_faclev), 
           Attr = c("Land use type")) %>% 
    na.omit()
}
qua_tree_div_long <- fun_cons_long(qua_tree_div)
qua_shrub_div_long <- fun_cons_long(qua_shrub_div)

# get p-values for box plots
fun_get_pvalue <- function(x) {
  y <- data.frame(Index = c("Abundance", "Richness", "Evenness"), 
                  Pvalue = NA, Label = NA) %>% 
    mutate(Index = factor(Index, levels = c("Abundance", "Richness", "Evenness")))
  j <- 0
  for (i in c("Abundance", "Richness", "Evenness")) {
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
tree_box_pvalue <- fun_get_pvalue(qua_tree_div)
shrub_box_pvalue <- fun_get_pvalue(qua_shrub_div)

# get box plots
fun_box_plot <- function(x, y, z) {
  ggplot(x, aes(Land_use_type, Index_value)) + 
    geom_boxplot() + 
    facet_grid(Index ~ Attr, scales = "free", 
               space = "free_x", switch = "both") + 
    scale_y_continuous(expand = expansion(mult = c(0.05,0.3))) +
    geom_text(data = y, aes(x =Inf, y = Inf, label = Label), 
              size=3.5, hjust = 1.05, vjust = 1.5) +
    theme(axis.text.x = element_text(angle = 90, size = 12)) + 
    labs(title = z, x = NULL, y = NULL) + 
    theme_bw()
}
ggarrange(fun_box_plot(qua_tree_div_long, tree_box_pvalue, "(a)"), 
          fun_box_plot(qua_shrub_div_long, shrub_box_pvalue, "(b)"))
rm(tree_box_pvalue, qua_tree_div_long, 
   shrub_box_pvalue, qua_shrub_div_long, 
   fun_box_plot, fun_cons_long, fun_get_pvalue)

# pairwise dunn test of indexes ~ land use type
fun_dunn <- function(x, taxa, index) {
  dunn_result <- dunn.test(x[ , index], x$Land_use_type, 
                           table = FALSE, kw = FALSE)
  x <- data.frame(
    "taxa" = taxa, 
    "index" = index, 
    "comparison" = dunn_result$comparisons, 
    "p" = dunn_result$P.adjusted
  ) %>% 
    separate(comparison, into = c("comparison_1", "comparison_2"), sep = " - ")
}
dunn_df_1 <- rbind(fun_dunn(qua_tree_div, "tree", "Abundance"), 
                   fun_dunn(qua_tree_div, "tree", "Richness"),
                   fun_dunn(qua_tree_div, "tree", "Evenness"),
                   fun_dunn(qua_shrub_div, "shrub", "Abundance"), 
                   fun_dunn(qua_shrub_div, "shrub", "Richness"),
                   fun_dunn(qua_shrub_div, "shrub", "Evenness")
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
    xlab(NULL) + ylab(NULL) + guides(fill = "none") + 
    facet_wrap(~ index, scales = "free", nrow = 1) +
    labs(title = title)
}
ggarrange(fun_dunn_plot(subset(dunn_df, taxa == "tree"), "Tree"), 
          fun_dunn_plot(subset(dunn_df, taxa == "shrub"), "Shrub"), 
          nrow = 2)
rm(dunn_df_1, dunn_df_2, dunn_df, 
   fun_dunn, fun_dunn_plot)

## Species composition ----
# use non-metric multidimensional scaling
set.seed(1234)
# nMDS calculation for tree
tree_mds_selected <- subset(qua_tree_div, Abundance > 1)
tree_mds_meta <- tree_mds_selected %>% 
  select(2:(number_tree_species+1)) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
tree_mds_meta$stress
stressplot(tree_mds_meta)
tree_mds_selected <- cbind(tree_mds_selected, tree_mds_meta$points)

# nMDS calculation for shrub
shrub_mds_selected <- qua_shrub_div %>% filter(Abundance > 5)
shrub_mds_meta <- shrub_mds_selected %>% 
  select(2:(number_shrub_species+1)) %>%
  metaMDS(distance = "bray", trace = FALSE, autotransform = FALSE) 
shrub_mds_meta$stress
stressplot(shrub_mds_meta)
shrub_mds_selected <- cbind(shrub_mds_selected, shrub_mds_meta$points)

# ANOSIM of trees and shrubs as labels for the nMDS plots
tree_anosim <- 
  anosim(tree_mds_selected[2:(number_tree_species+1)], 
         tree_mds_selected$Land_use_type)
shrub_anosim <- 
  anosim(shrub_mds_selected[2:(number_shrub_species+1)], 
         shrub_mds_selected$Land_use_type)

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
  fun_nmds_plot(tree_mds_selected, tree_hulls, "tree", 
                tree_mds_meta, tree_anosim), 
  fun_nmds_plot(shrub_mds_selected, shrub_hulls, "shrub", 
                shrub_mds_meta, shrub_anosim),
  common.legend = T, legend = "right"
)

# pairwise result of ANOSIM of trees by Land_use_type
anosim_pairs <- combn(Land_use_type_faclev, 2)
fun_anosim_pairs <- function(x, y) {
  result <- NULL
  for (i in 1:ncol(anosim_pairs)) {
    set.seed(1234)
    mds_selected_sub <- 
      subset(x, Land_use_type %in% c(anosim_pairs[1,i], anosim_pairs[2,i]))
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

# Cor among the indexes ----
chart.Correlation(subset(qua_tree_div, select = index_faclev))
chart.Correlation(subset(qua_shrub_div, select = index_faclev))


# Data for discussion ----
# means of quadrat Abundance and richness for trees
qua_tree_div %>% group_by(Land_use_type) %>% 
  dplyr::summarise(Abundance = mean(Abundance), Richness = mean(Richness))

