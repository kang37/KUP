setwd('C:/Rdata/KUP')
opar <- par(no.readonly = TRUE)

# 载入需要的R包
library(reshape)
library(vegan)
library(ggplot2)

# 读取数据，构造数据集
In.tree.data <- read.csv("In_tree_data.csv")
In.tree.data.inter1 <- data.frame(In.tree.data[,c(1:3)], value=1)
names(In.tree.data.inter1)[3] <- c("variable")
Tree.data <- cast(In.tree.data.inter1,Plot_ID~variable, sum)
Tree.data$Species_sum <- rowSums(Tree.data)
In.plotinfo <- read.csv("In_plot_info.csv")
Tree.data <- merge(Tree.data, In.plotinfo, by="Plot_ID")[,c(-153)]

In.tree.data.inter2 <- In.tree.data[,c("Plot_ID","Pub_Pri")]
In.tree.data.inter2$Pub_Pri_cat[In.tree.data.inter2$Pub_Pri=="N"] <- 0
In.tree.data.inter2$Pub_Pri_cat[In.tree.data.inter2$Pub_Pri=="Y"] <- 1
In.tree.data.inter3 <- aggregate(In.tree.data.inter2[,c(1,3)], by=list(In.tree.data.inter2$Plot_ID),FUN=mean)
In.tree.data.inter3$Pub_Pri_new[In.tree.data.inter3$Pub_Pri_cat==1] <- "Y"
In.tree.data.inter3$Pub_Pri_new[In.tree.data.inter3$Pub_Pri_cat==0] <- "N"
In.tree.data.inter3$Pub_Pri_new[In.tree.data.inter3$Pub_Pri_cat<1 &
                                  In.tree.data.inter3$Pub_Pri_cat>0] <- "NA"
Tree.data <- merge(Tree.data, In.tree.data.inter3, by="Plot_ID")[,c(-153,-154)]
rm(In.tree.data.inter1, In.tree.data.inter2, In.tree.data.inter3)

# 做物种等级丰富度曲线

## 提取物种列和Pub_Pri列，分为Pub和Pri两个部分，转化为物种频度表并排序
Species.pub <- In.tree.data[which(In.tree.data$Pub_Pri=="Y"), c("Plot_ID","Species_CN")]
Species.pri <- In.tree.data[which(In.tree.data$Pub_Pri=="N"), c("Plot_ID","Species_CN")]

## 将Pub树原始记录表转化成频数表并排序
Species.pub.freq <- as.data.frame(table(Species.pub[,"Species_CN"]))
Species.pub.freq <- Species.pub.freq[order(-Species.pub.freq$Freq),]
Species.pub.freq$Prop <- Species.pub.freq$Freq/nrow(Species.pub)
Species.pub.freq <- Species.pub.freq[which(Species.pub.freq$Prop>0),]

## 将Pri树原始记录表转化成频数表并排序
Species.pri.freq <- as.data.frame(table(Species.pri[,"Species_CN"]))
Species.pri.freq <- Species.pri.freq[order(-Species.pri.freq$Freq),]
Species.pri.freq$Prop <- Species.pri.freq$Freq/nrow(Species.pri)
Species.pri.freq <- Species.pri.freq[which(Species.pri.freq$Prop>0),]

## 删除无用变量
rm(Species.street, Species.nonstreet)

## 输出一般描述
cat("乔木物种数为：",length(unique(Species.tot$Species_CN)),"\n",
    "其中行道树物种数为：", length(unique(Species.pub$Species_CN)), "\n",
    "非行道树物种数为：", length(unique(Species.pri$Species_CN)))
cat("乔木植株数量为：", nrow(In.tree.data), "\n",
    "具体结构见右图：","\n")

par(mfrow=c(2,2))
barplot(table(In.tree.data$Pla_Spo), main = "种植或自生")
barplot(table(In.tree.data$Pot), main = "是否长盆里")
barplot(table(In.tree.data$Pub_Pri), main = "是否公共")
barplot(table(In.tree.data$Street), main = "是否行道树")
par(opar)

## 作图
par(mfrow=c(1,2))
barplot(Species.pub.freq$Freq, xlab = "Rank of public trees", 
        ylim = c(0.8,100), xlim = c(0,130), space = 0,
        ylab = "Relative abundance", log = "y")
barplot(Species.pri.freq$Freq,
        ylim = c(0.8,100), xlim = c(0,130), space = 0,
        xlab = "Rank of private trees", log = "y")
par(opar)

# 作MDS分析
library(vegan)
library(ggplot2)
library(reshape)

## 计算相异度矩阵
### 整合数据，构建包含样地属性的群落数据
In.tree.datapreana <- In.tree.data[,c(1:3)]
In.tree.datapreana$Value <- 1
names(In.tree.datapreana)[3:4] <- c("variable","value")
Tree.comm.attr <- cast(In.tree.datapreana, Plot_ID ~ variable, sum)

### 筛选出符合要求的数据并构建用于MDS计算的群落数据
Tree.comm.attr$Sum <- rowSums(Tree.comm.attr[,c(2:143)])

In.Plotinfo <- read.csv("In_plot_info.csv")
Tree.comm.attr <- merge(Tree.comm.attr, In.Plotinfo, by="Plot_ID")
Tree.comm <- Tree.comm.attr[which(Tree.comm.attr$Sum>1),]

Tree.mds <- metaMDS(comm = Tree.comm[,c(2:143)], distance = "bray", trace = FALSE, autotransform = FALSE)
Tree.mds.xy <- Tree.mds$points
plot(Tree.mds.xy)

Tree.mds.xy <- as.data.frame(Tree.mds.xy)
Tree.mds.xy$Plot_ID <- Tree.comm$Plot_ID
Tree.mds.xy$Landuse <- Tree.comm$Landuse_class
Tree.mds.xy$Landuse_agg <- "NA"
  Tree.mds.xy$Landuse_agg[Tree.mds.xy$Landuse=="R low"] <- "R"
  Tree.mds.xy$Landuse_agg[Tree.mds.xy$Landuse=="R high"] <- "R"
  Tree.mds.xy$Landuse_agg[Tree.mds.xy$Landuse=="R resi"] <- "R"
  Tree.mds.xy$Landuse_agg[Tree.mds.xy$Landuse=="Ind"] <- "Ind"
  Tree.mds.xy$Landuse_agg[Tree.mds.xy$Landuse=="Com neigh"] <- "Com"
  Tree.mds.xy$Landuse_agg[Tree.mds.xy$Landuse=="Com"] <- "Com"
Tree.mds.xy$Ward <- Tree.comm$Ward_EN
tree.mds.xy <- na.omit(merge(Tree.mds.xy, Tree.data, by="Plot_ID"))

## 作图
ggplot(Tree.mds.xy, aes(MDS1, MDS2, color = Tree.mds.xy$Landuse)) + geom_point() + theme_bw()
ggplot(Tree.mds.xy, aes(MDS1, MDS2, color = Tree.mds.xy$Ward)) + geom_point() + theme_bw()
ggplot(Tree.mds.xy, aes(MDS1, MDS2, color = Tree.mds.xy$Landuse_agg)) + geom_point() + theme_bw()
ggplot(Tree.mds.xy, aes(MDS1, MDS2, color = tree.mds.xy$Pub_Pri_new )) + geom_point() + theme_bw()
### 均无明显区分？















