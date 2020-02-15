setwd('C:/Rdata/KUP')
opar <- par(no.readonly = TRUE)

### 载入需要的R包
library(reshape)
library(vegan)
library(ggplot2)
library(ggpubr)



### 构造乔木数据集
# 读取数据
In.plotinfo <- read.csv("In_plot_info.csv")
In.tree.data <- read.csv("In_tree_data.csv")
In.shrub.data <- read.csv("In_shrub_data.csv")
In.plantinfo <- read.csv("In_plant_info.csv")
#? 删除异常点
In.tree.data <- subset(In.tree.data, Plot_ID!=67)
In.shrub.data <- subset(In.shrub.data, Plot_ID!=67)
# 编入植物分类数据
In.tree.data <- merge(In.tree.data, In.plantinfo, by="Species_CN")
In.shrub.data <- merge(In.shrub.data, In.plantinfo, by="Species_CN")

## 乔木
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

## 灌木
# 将输入数据变成物种-数量矩阵
In1 <- data.frame(In.shrub.data[,c("Plot_ID", "Plant_ID", "Species_CN")], value=In.shrub.data$Area)
names(In1)[3] <- c("variable")
Shrub.comm <- cast(In1,Plot_ID~variable, sum)
rm(In1)
# 添加生物多样性指数
Shrub.comm$Species_sum <- rowSums(Shrub.comm)
Shrub.comm$Richness <- apply(Shrub.comm[,c(2:ncol(Shrub.comm))]>0,1,sum)
Shrub.comm$Shannon <- diversity(Shrub.comm[,c(2:ncol(Shrub.comm))], index = "shannon")
Shrub.comm$Simpson <- diversity(Shrub.comm[,c(2:ncol(Shrub.comm))], index = "simpson")
Shrub.comm$Evenness <- diversity(Shrub.comm[, c(2:ncol(Shrub.comm))], index = "shannon")/log(Shrub.comm$Richness)
# 添加样地信息
Shrub.comm <- merge(Shrub.comm, In.plotinfo, by="Plot_ID")
# 完善样地信息：将Dist分级
Shrub.comm$Dist_level <- NA
Shrub.comm$Dist_level[Shrub.comm$Dist<3000] <- "Center"
Shrub.comm$Dist_level[and(Shrub.comm$Dist>=3000, Shrub.comm$Dist<6000)] <- "Cen-Mid"
Shrub.comm$Dist_level[and(Shrub.comm$Dist>=6000, Shrub.comm$Dist<9000)] <- "Mid-Sub"
Shrub.comm$Dist_level[Shrub.comm$Dist>9000] <- "Suburban"
# 完善样地信息：种的是Pub树还是Pri树或者都有？
In1 <- In.shrub.data[,c("Plot_ID","Pub_Pri")]
In1$Pub_Pri_cat[In1$Pub_Pri=="N"] <- 0
In1$Pub_Pri_cat[In1$Pub_Pri=="Y"] <- 1
In2 <- aggregate(In1[,c(1,3)], by=list(In1$Plot_ID),FUN=mean)
In2$Pub_Pri_new[In2$Pub_Pri_cat==1] <- "Y"
In2$Pub_Pri_new[In2$Pub_Pri_cat==0] <- "N"
In2$Pub_Pri_new[In2$Pub_Pri_cat<1 &
                  In2$Pub_Pri_cat>0] <- NA
Shrub.comm <- merge(Shrub.comm, In2, by="Plot_ID")
Shrub.comm$Group.1 <- NULL
Shrub.comm$Pub_Pri_cat <- NULL
rm(In1, In2)



### 分析开始
### 乔灌一般描述
cat("乔灌物种数为：",length(unique(union(In.tree.data$Species_CN,In.shrub.data$Species_CN))),
    "乔木特有物种数为：", length(setdiff(unique(In.tree.data$Species_CN), unique(In.shrub.data$Species_CN))),
    "灌木特有物种数为：", length(setdiff(unique(In.shrub.data$Species_CN), unique(In.tree.data$Species_CN))),
    "乔灌共有物种数为：", length(intersect(unique(In.tree.data$Species_CN), unique(In.shrub.data$Species_CN))))
cat("从数量来看乔木物种数为：",length(unique(In.tree.data$Species_CN)),"\n",
    "乔木植株数量为：", nrow(In.tree.data),"\n",
    "灌木物种数为：",length(unique(In.shrub.data$Species_CN)),"\n",
    "灌木面积为：", sum(In.shrub.data$Area))

cat("物种最多的科是" )
In1 <- as.data.frame(table(In.plantinfo$Family))
head(In1[order(-In1$Freq),])
cat("但是数量最多的科是：")
In1 <- as.data.frame(table(In.tree.data$Family))
head(In1[order(-In1$Freq),])

cat("从物种数角度看外来植物占：")
In1 <- as.data.frame(table(In.plantinfo$Nt.Ex))
In1
In1$Freq[In1$Var1=="ex"]/sum(In1$Freq)
cat("其中，乔木物种中外来种占：")
In1 <- as.data.frame(table(unique(In.tree.data[,c("Species_CN", "Nt.Ex")])$Nt.Ex))
In1
In1$Freq[In1$Var1=="ex"]/sum(In1$Freq)
cat("而在灌木物种中外来种占：")
In1 <- as.data.frame(table(unique(In.shrub.data[,c("Species_CN", "Nt.Ex")])$Nt.Ex))
In1
In1$Freq[In1$Var1=="ex"]/sum(In1$Freq)
cat("但是从植株数量角度看，外来植物在乔木中占：")
In1 <- as.data.frame(table(In.tree.data$Nt.Ex))
In1
In1$Freq[In1$Var1=="ex"]/sum(In1$Freq)
cat("外来植物在灌木中面积占比：")
In1 <- data.frame("Var1"=c("ex","nt"), "Freq"=c(sum(In.shrub.data$Area[which(In.shrub.data$Nt.Ex=="ex")]), sum(In.shrub.data$Area[which(In.shrub.data$Nt.Ex=="nt")])))
In1
In1$Freq[In1$Var1=="ex"]/sum(In1$Freq)

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


### 乔灌植物性质作图
par(mfrow=c(4,2))
barplot(c("Planted"=nrow(In.tree.data[which(In.tree.data$Pla_Spo=="P"),]), 
          "Spontaneous"=nrow(In.tree.data[which(In.tree.data$Pla_Spo=="S"),])), main = "Tree")
barplot(c("Planted"=sum(In.shrub.data$Area[In.shrub.data$Pla_Spo=="P"]), 
          "Spontaneous"=sum(In.shrub.data$Area[In.shrub.data$Pla_Spo=="S"])), main = "Shrub")

barplot(c("Non-pot"=nrow(In.tree.data[which(In.tree.data$Pot=="N"),]), 
          "Pot"=nrow(In.tree.data[which(In.tree.data$Pot=="Y"),])))
barplot(c("Non-pot"=sum(In.shrub.data$Area[In.shrub.data$Pot=="N"]), 
          "pot"=sum(In.shrub.data$Area[In.shrub.data$Pot=="Y"])))

barplot(c("Private"=nrow(In.tree.data[which(In.tree.data$Pub_Pri=="N"),]), 
          "Public"=nrow(In.tree.data[which(In.tree.data$Pub_Pri=="Y"),])))
barplot(c("Privat"=sum(In.shrub.data$Area[In.shrub.data$Pub_Pri=="N"]), 
          "Public"=sum(In.shrub.data$Area[In.shrub.data$Pub_Pri=="Y"])))

barplot(c("Non street"=nrow(In.tree.data[which(In.tree.data$Street=="N"),]), 
          "Street"=nrow(In.tree.data[which(In.tree.data$Street=="Y"),])))
barplot(c("Non street"=sum(In.shrub.data$Area[In.shrub.data$Street=="N"]),
          "Street"=sum(In.shrub.data$Area[In.shrub.data$Street=="Y"])))
par(opar)



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



### 各种指数之间的关系。
library(car)
library(corrplot)

scatterplotMatrix(Tree.comm[,c("Species_sum", "Richness", "Shannon", "Simpson", "Evenness")])
scatterplotMatrix(Shrub.comm[,c("Species_sum", "Richness", "Shannon", "Simpson", "Evenness")])

# 乔灌各种指数相关系数
cor(Tree.comm[,c("Species_sum", "Richness", "Shannon", "Simpson", "Evenness")], use = "pairwise.complete.obs")
cor(Shrub.comm[,c("Species_sum", "Richness", "Shannon", "Simpson", "Evenness")], use = "pairwise.complete.obs")

# 乔灌相关系数热图
library(PerformanceAnalytics)
par(mfrow=c(2,2))
corrplot(cor(Tree.comm[,c("Species_sum", "Richness", "Shannon", "Simpson", "Evenness")], use = "pairwise.complete.obs"), method = "color", type="upper", addCoef.col = "black", sig.level = 0.01, insig = "label_sig")
corrplot(cor(Shrub.comm[,c("Species_sum", "Richness", "Shannon", "Simpson", "Evenness")], use = "pairwise.complete.obs"), method = "color", type="upper", addCoef.col = "black", insig = "label_sig")
chart.Correlation(Tree.comm[,c("Species_sum", "Richness", "Shannon", "Simpson", "Evenness")])
chart.Correlation(Shrub.comm[,c("Species_sum", "Richness", "Shannon", "Simpson", "Evenness")])
par(opar)
#? 分格对后面两图不适用

### 指数~属性
# 正态性检验
# 乔木
for (i in c("Species_sum", "Richness", "Shannon","Simpson", "Evenness")) {
  print(i)
  In2 <- NULL
  for (j in c("Ward", "Dist_level", "Landuse_class", "Pub_Pri_new")) {
    In1 <- data.frame("Attr"=NA, "Norm"=0, "Non-norm"=0)
    In3 <- NULL
    In3 <- tapply(Tree.comm[,i],Tree.comm[,j], shapiro.test)
    for (k in c(1:length(In3))) {
      In1$Attr <- j
      if (In3[[k]]$p.value>0.05) {
        In1$Norm <- In1$Norm+1
      } else {
        In1$Non.norm <- In1$Non.norm+1
      }
    }
    In2 <- rbind(In2,In1)
  }
  print(In2)
  rm(In1, In2, In3)
}
# 灌木
#? 出错提示

# 非参数检验
# 总体比较
# 乔木
In4 <- data.frame("attr"=c("Ward", "Dist_level", "Landuse_class", "Pub_Pri_new"))
for (i in c("Species_sum", "Richness", "Shannon","Simpson", "Evenness")) {
  In3 <- data.frame("attr"=as.character(), "p"=as.numeric())
  
  for (j in c("Ward", "Dist_level", "Landuse_class", "Pub_Pri_new")) {
    In1 <- kruskal.test(Tree.comm[,i] ~ Tree.comm[,c(j)], data = Tree.comm)
    In2 <- data.frame("attr"=j,"p"=In1$p.value)
    In3 <- rbind(In3,In2)
  }
  In3$sig <- NA
  In3$sig[In3$p<0.001] <- "***"
  In3$sig[and(In3$p>0.001, In3$p<0.01)] <- "**"
  In3$sig[and(In3$p>0.01, In3$p<0.05)] <- "*"
  In3$sig[In3$p>0.05] <- "."
  In4 <- cbind(In4,In3[,c(2:3)])
  colnames(In4)[ncol(In4)-1] <- i
}
print(In4)

# 灌木
for (i in c("Species_sum","Shannon","Simpson", "Evenness")) {
  print(i)
  In3 <- data.frame("attr"=as.character(), "p"=as.numeric())
  for (j in c("Ward", "Dist_level", "Landuse_class", "Pub_Pri_new")) {
    In1 <- kruskal.test(Shrub.comm[,i] ~ Shrub.comm[,c(j)], data = Shrub.comm)
    In2 <- data.frame("attr"=j,"p"=In1$p.value)
    In3 <- rbind(In3,In2)
  }
  In3$sig <- NA
  In3$sig[In3$p<0.001] <- "***"
  In3$sig[and(In3$p>0.001, In3$p<0.01)] <- "**"
  In3$sig[and(In3$p>0.01, In3$p<0.05)] <- "*"
  In3$sig[In3$p>0.05] <- ""
  print(In3)
}


# 两两比较
# 乔木
library(pheatmap)
for (i in c("Species_sum", "Richness", "Shannon","Simpson", "Evenness")) {
  print(i)
  for (j in c("Ward", "Dist_level", "Landuse_class", "Pub_Pri_new")) {
    print(j)
    In1 <- pairwise.wilcox.test(Tree.comm[,i], Tree.comm[,j], p.adjust.method = "BH", exact=FALSE)
    In2 <- as.data.frame(In1$p.value)
    In2[is.na(In2)] <- ""
    In2[In2>0.05] <- "."
    In2[In2>0.01] <- "*"
    In2[In2>0.001] <- "**"
    In2[In2>0] <- "***"
    print(In2)
    write.table(cbind(c("", rownames(In2)), rbind(colnames(In2), In2)), "in2.csv",sep = c(","), append = TRUE, 
                row.names = FALSE, col.names = FALSE)
  }
}

# 灌木
for (j in c("Species_sum","Shannon","Simpson", "Evenness")) {
  print(j)
  for (i in c("Ward", "Dist_level", "Landuse_class", "Pub_Pri_new")) {
    print(i)
    Inter2 <- pairwise.wilcox.test(Shrub.comm[,j], Shrub.comm[,i], 
                                   p.adjust.method = "BH", exact=FALSE)
    print(Inter2)
    boxplot(Shrub.comm[,j]~Shrub.comm[,i], ylab = j, xlab = i)
  }
}
par(opar)


par(mfrow=c(4,4))
for (j in c("Species_sum","Shannon","Simpson", "Evenness")) {
  for (i in c("Ward", "Dist_level", "Landuse_class", "Pub_Pri_new")) {
    boxplot(Tree.comm[,j]~Tree.comm[,i], ylab = j, xlab = i, par(las=2))
  }
}

par(mfrow=c(4,4))
for (j in c("Species_sum","Shannon","Simpson", "Evenness")) {
  for (i in c("Ward", "Dist_level", "Landuse_class", "Pub_Pri_new")) {
    boxplot(Shrub.comm[,j]~Shrub.comm[,i], ylab = j, xlab = NULL, par(las=2), cex.axis=1.5)
  }
}






