setwd('C:/Rdata/KUP')
opar <- par(no.readonly = TRUE)

# 载入需要的R包
library(reshape)
library(vegan)
library(ggplot2)
library(ggpubr)


### 构造数据集
# 读取数据
In.plotinfo <- read.csv("In_plot_info.csv")
In.tree.data <- read.csv("In_tree_data.csv")
#? 删除异常点
In.tree.data <- subset(In.tree.data, Plot_ID!=67)
# 将输入数据变成物种-数量矩阵
In1 <- data.frame(In.tree.data[,c("Plot_ID", "Plant_ID", "Species_CN")], value=1)
names(In1)[3] <- c("variable")
Tree.comm <- cast(In1,Plot_ID~variable, sum)
rm(In1)
# 添加生物多样性指数
Tree.comm$Species_sum <- rowSums(Tree.comm)
Tree.comm$Richness <- apply(Tree.comm[,c(2:143)]>0,1,sum)
Tree.comm$Shannon <- diversity(Tree.comm[,c(2:143)], index = "shannon")
Tree.comm$Simpson <- diversity(Tree.comm[,c(2:143)], index = "simpson")
Tree.comm$Evenness <- diversity(Tree.comm[, c(2:143)], index = "shannon")/log(Tree.comm$Richness)
# 添加样地信息
Tree.comm <- merge(Tree.comm, In.plotinfo, by="Plot_ID")
# 完善样地信息：将Dist分级
Tree.comm$Dist_level <- NA
Tree.comm$Dist_level[Tree.comm$Dist<3000] <- "Center"
Tree.comm$Dist_level[and(Tree.comm$Dist>=3000, Tree.comm$Dist<6000)] <- "Cen-Mid"
Tree.comm$Dist_level[and(Tree.comm$Dist>=6000, Tree.comm$Dist<9000)] <- "Mid-Sub"
Tree.comm$Dist_level[Tree.comm$Dist>9000] <- "Suburban"
# 完善样地信息：种的是Pub树还是Pri树或者都有？
In1 <- In.tree.data[,c("Plot_ID","Pub_Pri")]
In1$Pub_Pri_cat[In1$Pub_Pri=="N"] <- 0
In1$Pub_Pri_cat[In1$Pub_Pri=="Y"] <- 1
In2 <- aggregate(In1[,c(1,3)], by=list(In1$Plot_ID),FUN=mean)
In2$Pub_Pri_new[In2$Pub_Pri_cat==1] <- "Y"
In2$Pub_Pri_new[In2$Pub_Pri_cat==0] <- "N"
In2$Pub_Pri_new[In2$Pub_Pri_cat<1 &
                  In2$Pub_Pri_cat>0] <- NA
Tree.comm <- merge(Tree.comm, In2, by="Plot_ID")
Tree.comm$Group.1 <- NULL
Tree.comm$Pub_Pri_cat <- NULL
rm(In1, In2)


### 分析开始
### 一般描述
cat("乔木物种数为：",length(unique(In.tree.data$Species_CN)),"\n",
    "乔木植株数量为：", nrow(In.tree.data), "\n",
    "具体结构见右图：","\n")
# 并且作图
par(mfrow=c(2,2))
barplot(table(In.tree.data$Pla_Spo), main = "种植或自生")
barplot(table(In.tree.data$Pot), main = "是否长盆里")
barplot(table(In.tree.data$Pub_Pri), main = "是否公共")
barplot(table(In.tree.data$Street), main = "是否行道树")
par(opar)

### 各种指数之间的关系。
library(car)
library(corrplot)
scatterplotMatrix(Tree.comm[,c("Species_sum", "Richness", "Shannon", "Simpson", "Evenness")])
cor(Tree.comm[,c("Species_sum", "Richness", "Shannon", "Simpson", "Evenness")], use = "pairwise.complete.obs")
corrplot(cor(Tree.comm[,c("Species_sum", "Richness", "Shannon", "Simpson", "Evenness")], use = "pairwise.complete.obs"), method = "color", type="upper", addCoef.col = "black")

### MDS
# 计算相异度矩阵
# 整合数据，构建包含样地属性的群落数据
Tree.mds <- Tree.comm[,c(1:143)][which(rowSums(Tree.comm[,c(2:143)])>1),]
Tree.mds.result <- metaMDS(comm = Tree.mds[,c(2:143)], distance = "bray", trace = FALSE, autotransform = FALSE)
Tree.mds.result <- cbind(Tree.mds, Tree.mds.result$points)
Tree.mds.result <- merge(Tree.mds.result, In.plotinfo, by = "Plot_ID")
Tree.mds.result <- merge(Tree.mds.result, Tree.comm[,c("Plot_ID","Pub_Pri_new","Dist_level")], by = "Plot_ID")
rm(Tree.mds)
# 作图
P1 <- ggplot(Tree.mds.result, aes(MDS1, MDS2, color = Tree.mds.result$Ward)) + 
  geom_point() + theme_bw()+
  stat_ellipse(aes(x=MDS1,y=MDS2,color= Tree.mds.result$Ward))
P2 <- ggplot(Tree.mds.result, aes(MDS1, MDS2, color = Tree.mds.result$Dist_level)) + 
  geom_point() + theme_bw()+
  stat_ellipse(aes(x=MDS1,y=MDS2,color= Tree.mds.result$Dist_level))
P3 <- ggplot(Tree.mds.result, aes(MDS1, MDS2, color = Tree.mds.result$Landuse_class)) + 
  geom_point() + theme_bw()+
  stat_ellipse(aes(x=MDS1,y=MDS2,color= Tree.mds.result$Landuse_class))
P4 <- ggplot(Tree.mds.result, aes(MDS1, MDS2, color = Tree.mds.result$Landuse_agg)) + 
  geom_point() + theme_bw()+
  stat_ellipse(aes(x=MDS1,y=MDS2,color= Tree.mds.result$Landuse_agg))
P5 <- ggplot(na.omit(Tree.mds.result), aes(MDS1, MDS2, color = na.omit(Tree.mds.result)$Pub_Pri_new)) + 
  geom_point() + theme_bw()+
  stat_ellipse(aes(x=MDS1,y=MDS2,color= na.omit(Tree.mds.result)$Pub_Pri_new))
ggarrange(P1, P2, P3, P4, P5, ncol = 3, nrow = 2)
rm(P1, P2, P3, P4, P5)
# 定量分析质心是否有差异
for (i in c("Ward", "Dist_level", "Landuse_class", "Landuse_agg", "Pub_Pri_new")) {
  print(c(i))
  if (sum(is.na(Tree.mds.result[,i])>0)) {
    In <- adonis(na.omit(Tree.mds.result)[,c(2:143)]~na.omit(Tree.mds.result)[,i], method = "bray", permutations = 999)
  } else {
    In <- adonis(Tree.mds.result[,c(2:143)]~Tree.mds.result[,i], method = "bray", permutations = 999)
  }
  print(In$aov.tab$`Pr(>F)`)
}
library(RVAideMemoire)
rm(In)
# 每次运行结果不同且无法通过set.seed解决，只能报告显著性






### 物种等级丰富度曲线：对上个分析中质心差异显著的组合？
# 提取物种列和Pub_Pri列，分为Pub和Pri两个部分，转化为物种频度表并排序
In.tree.data.pub <- In.tree.data[which(In.tree.data$Pub_Pri=="Y"), c("Plot_ID","Species_CN")]
In.tree.data.pri <- In.tree.data[which(In.tree.data$Pub_Pri=="N"), c("Plot_ID","Species_CN")]
# 将Pub树原始记录表转化成频数表并排序
In.tree.data.pub.freq <- as.data.frame(table(In.tree.data.pub[,"Species_CN"]))
In.tree.data.pub.freq <- In.tree.data.pub.freq[order(-In.tree.data.pub.freq$Freq),]
In.tree.data.pub.freq$Prop <- In.tree.data.pub.freq$Freq/nrow(In.tree.data.pub)
In.tree.data.pub.freq <- In.tree.data.pub.freq[which(In.tree.data.pub.freq$Prop>0),]
# 将Pri树原始记录表转化成频数表并排序
In.tree.data.pri.freq <- as.data.frame(table(In.tree.data.pri[,"Species_CN"]))
In.tree.data.pri.freq <- In.tree.data.pri.freq[order(-In.tree.data.pri.freq$Freq),]
In.tree.data.pri.freq$Prop <- In.tree.data.pri.freq$Freq/nrow(In.tree.data.pri)
In.tree.data.pri.freq <- In.tree.data.pri.freq[which(In.tree.data.pri.freq$Prop>0),]
# 作图
par(mfrow=c(1,2))
barplot(In.tree.data.pub.freq$Freq, xlab = "Rank of public trees", 
        ylim = c(0.8,100), xlim = c(0,130), space = 0,
        ylab = "Relative abundance", log = "y")
barplot(In.tree.data.pri.freq$Freq,
        ylim = c(0.8,100), xlim = c(0,130), space = 0,
        xlab = "Rank of private trees", log = "y")
par(opar)
rm(In.tree.data.pub, In.tree.data.pub.freq, In.tree.data.pri, In.tree.data.pri.freq)
# 可见Pri植株多，且更加均匀。

# 可以用方差分析吗？先看正态性
for (k in c("Species_sum","Shannon","Simpson")) {
  print(k)
  for (j in c("Ward", "Landuse_class", "Landuse_agg", "Pub_Pri_new")) {
    print(j)
    inter1 <- NULL
    Inter1 <- tapply(Tree.comm[,k],Tree.comm[,j], shapiro.test)
    Inter2 <- NULL
    for (i in c(1:length(Inter1))) {
      Inter2 <- rbind(Inter2, Inter1[[i]]$p.value>0.05)
    }
    print(table(Inter2))
  }
}

# 几乎都不符合正态分布，所以检查各项的log变换是否符合正态分布。

Inter1 <- tapply(log(Tree.comm$Species_sum),Tree.comm$Ward, shapiro.test)
Inter2 <- NULL
for (i in c(1:length(Inter1))) {
  Inter2 <- rbind(Inter2, Inter1[[i]]$p.value>0.05)
}
table(Inter2)

Inter1 <- tapply(log(Tree.comm$Species_sum),Tree.comm$Landuse_subclass, shapiro.test)
Inter2 <- NULL
for (i in c(1:length(Inter1))) {
  Inter2 <- rbind(Inter2, Inter1[[i]]$p.value>0.05)
}
table(Inter2)

Inter1 <- tapply(log(Tree.comm$Species_sum),Tree.comm$Landuse_class, shapiro.test)
Inter2 <- NULL
for (i in c(1:length(Inter1))) {
  Inter2 <- rbind(Inter2, Inter1[[i]]$p.value>0.05)
}
table(Inter2)

Inter1 <- tapply(log(Tree.comm$Species_sum),Tree.comm$Landuse_agg, shapiro.test)
Inter2 <- NULL
for (i in c(1:length(Inter1))) {
  Inter2 <- rbind(Inter2, Inter1[[i]]$p.value>0.05)
}
table(Inter2)

Inter1 <- tapply(log(Tree.comm$Species_sum),Tree.comm$Pub_Pri_new, shapiro.test)
Inter2 <- NULL
for (i in c(1:length(Inter1))) {
  Inter2 <- rbind(Inter2, Inter1[[i]]$p.value>0.05)
}
table(Inter2)

# 删除无关变量。
rm(Inter1, Inter2, i)

# 结果参差不齐，用非参数检验吧。
library(dunn.test)
for (j in c("Species_sum","Shannon","Simpson", "Evenness")) {
  print(j)
  par(mfrow=c(2,3))
  for (i in c("Ward", "Dist_level", "Landuse_class", "Landuse_agg", "Pub_Pri_new")) {
    print(i)
    Inter1 <- kruskal.test(Tree.comm[,j] ~ Tree.comm[,c(i)], data = Tree.comm)
    print(Inter1)
    Inter2 <- pairwise.wilcox.test(Tree.comm[,j], Tree.comm[,i], 
                                   p.adjust.method = "BH", exact=FALSE)
    print(Inter2)
    boxplot(Tree.comm[,j]~Tree.comm[,i], ylab = j, xlab = i)
  }
  par(opar)
}




