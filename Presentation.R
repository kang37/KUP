library(openxlsx)
library(ggplot2)
library(ggpubr)

citycomp <- read.xlsx("Compare_cities_diversity.xlsx", "plot")
citycomp_richness <- citycomp[c("City", "Land.use.class", "Richness")]
citycomp_richness$Richness <- as.numeric(citycomp_richness$Richness)
citycomp_richness <- citycomp_richness[which(
  citycomp_richness$Land.use.class %in% 
    c("residential", "neigh-commercial", "commercial")), ]


citycomp_density <- citycomp[c("City", "Land.use.class", "Density")]
citycomp_density$Density <- as.numeric(citycomp_density$Density)
citycomp_density <- citycomp_density[which(
  citycomp_density$Land.use.class %in% 
    c("residential", "neigh-commercial", "commercial")), ]


ggarrange(plotlist = list(ggplot(citycomp_richness, aes(x = Land.use.class, y = Richness)) + geom_col(aes(fill = City), width = 0.5, position = "dodge"), 
                       ggplot(citycomp_density, aes(x = Land.use.class, y = Density)) + geom_col(aes(fill = City), width = 0.5, position = "dodge")), 
          common.legend = TRUE, ncol = 1)
