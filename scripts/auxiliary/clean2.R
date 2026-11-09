######################Friends分析######################==
rm(list = ls())   #清除环境变量
options(stringsAsFactors = FALSE)   #禁止chr转成factor

#1. R包和数据准备
##1.1 R包安装与加载
###R包安装
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("GOSemSim")
#install.packages("ggplot2")
BiocManager::install("org.Hs.eg.db")

###R包加载
library(org.Hs.eg.db)
library(GOSemSim)
library(reshape2)
library(ggplot2)
# 安装 clusterProfiler 包（如果未安装）
install.packages("BiocManager")
BiocManager::install("clusterProfiler")

# 加载 clusterProfiler 包
library(clusterProfiler)


##1.2 数据读取
rt <- read.table("gene.txt", sep="\t", header = T, check.names = F)
str(rt)
 
getwd()
setwd("/mnt/data/")



###字符串转换
rt$ENTREZID <- as.character(rt$ENTREZID)
head(rt)

#2. 计算相似性
bp <- godata('org.Hs.eg.db', ont="BP", computeIC = FALSE)
cc <- godata('org.Hs.eg.db', ont="CC", computeIC = FALSE)
mf <- godata('org.Hs.eg.db', ont="MF", computeIC = FALSE)

simbp <- mgeneSim(rt$ENTREZID,
                  semData = bp,
                  measure = "Wang",
                  drop = NULL,
                  combine = "BMA")
simcc <- mgeneSim(rt$ENTREZID,
                  semData = cc,
                  measure = "Wang",
                  drop = NULL,
                  combine = "BMA")
simmf <- mgeneSim(rt$ENTREZID,
                  semData = mf,
                  measure = "Wang",
                  drop = NULL,
                  combine = "BMA")

#计算基因的几何平均值
fsim <- (simmf * simcc * simbp)^(1/3)

#将基因的名字由ENTREZID改为gene symbol
colnames(fsim) = rt$SYMBOL
rownames(fsim) = rt$SYMBOL

#将基因自己和自己的相似度设为NA
for (i in 1:ncol(fsim)){
  fsim[i,i] <- NA
}

dat <- melt(fsim)   #把宽格式数据转化成长格式
dat <- dat[!is.na(dat$value),] #删掉带NA的行
dat <- dat[,c(1,3)]
head(dat)

#3. 绘制boxplot图
#计算每个基因跟其他基因相似度的平均值
dat.mean <- aggregate(value~Var1, dat, mean) 

m <- dat.mean$value
names(m) <- dat.mean$Var1
#按平均值给基因名排序
dat$Var1 <- factor(dat$Var1,
                   levels=names(sort(m)))
str(dat)

##绘图
ggplot(dat,
       aes(x = Var1, y = value, fill = factor(Var1))) + 
  scale_fill_brewer(palette="Set2") +   #配色选择
  geom_boxplot() +
  coord_flip() +   #坐标轴互换
  xlab("") + ylab("") +
  theme_bw()
