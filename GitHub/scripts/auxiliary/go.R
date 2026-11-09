###############GO分析#################################
rm(list=ls())   #清空环境变量
options(stringsAsFactors = F)
load("TCGA-STADDEG_symbol.Rdata")

#1. 添加ENTREZID列
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
gene_up = DEG[DEG$change == 'UP','ENTREZID'] 
gene_down = DEG[DEG$change == 'DOWN','ENTREZID'] 
gene_diff = c(gene_up,gene_down)

#2. GO分析
#2.1 GO分析
ego <- enrichGO(gene = gene_diff,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
write.table(ego,file="GO.txt",sep="\t",quote=F,row.names = F)
##2.2 细胞组分
ego_CC <- enrichGO(gene = gene_diff,
                   OrgDb= org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
write.table(ego_CC,file="GO_CC.txt",sep="\t",quote=F,row.names = F)
##2.3 生物过程
ego_BP <- enrichGO(gene = gene_diff,
                   OrgDb= org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
write.table(ego_BP,file="GO_BP.txt",sep="\t",quote=F,row.names = F)
##2.4 分子功能
ego_MF <- enrichGO(gene = gene_diff,
                   OrgDb= org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
write.table(ego_MF,file="GO_MF.txt",sep="\t",quote=F,row.names = F)
save(ego,ego_CC,ego_BP,ego_MF,file = "ego_TCGA-STAD.Rdata")

#3. 可视化
#load("ego_TCGA-STAD.Rdata")
##3.1 柱状图
barplot(ego, showCategory = 20,color = "pvalue")
##3.2 气泡图
dotplot(ego, showCategory = 20)
##3.3 分类展示
barplot(ego, drop = TRUE, showCategory =10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
dotplot(ego,showCategory = 10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
##3.4 网络图
List = DEG$log2FoldChange
names(List)= DEG$ENTREZID
List = sort(List,decreasing = T)
cnetplot(ego, categorySize="pvalue", foldChange = List,colorEdge = TRUE)
##3.5 网络图2
cnetplot(ego, foldChange = List, circular = TRUE, colorEdge = TRUE)
##3.6 富集分布图
emapplot(ego)
##3.7 热图
heatplot(ego,foldChange = List)
##3.8 upset图

