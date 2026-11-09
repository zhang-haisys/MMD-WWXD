if(class(exp)[1]!="matrix")exp=as.matrix(exp)
if(length(unique(group))!=2)warning('分组不是2组，请检查')
if(!(case%in%group&control%in%group))stop('case和control与group中的信息未对应上')
if(!all(colnames(exp)==names(group)))stop('样本名未对应上')
if(sum(is.na(exp))>=1)stop('表达矩阵有空值')
if(is.null(colors))stop("缺少颜色参数")

####差异分析####
dir.create("C:/Users/18025/Desktop/EGFR/2.DEG/DEG", recursive = TRUE)
setwd("C:/Users/18025/Desktop/EGFR/2.DEG/DEG/")

if(is_count){
  cat("你选择了使用count矩阵计算差异基因，这会非常耗时，另外请检查数据是否为count矩阵")
  if(max(table(group))>15){  
    ####edgeR####
    library(edgeR)
    
    group_deg <- factor(group,levels = c(control,case))
    design <- model.matrix(~0+group_deg)
    row.names(design) <- colnames(exp)
    
    exp[exp<0] <- 0
    DGElist <- DGEList(counts = exp,group = group_deg)
    keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2
    
    DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
    DGElist <- calcNormFactors( DGElist )
    DGElist <- estimateGLMCommonDisp(DGElist, design)
    DGElist <- estimateGLMTrendedDisp(DGElist, design)
    DGElist <- estimateGLMTagwiseDisp(DGElist, design)
    
    fit <- glmFit(DGElist, design)
    results <- glmLRT(fit, contrast = c(-1, 1)) 
    nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
    nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
    
    dif <- na.omit(nrDEG_edgeR)
    dif$gene_symbol <- row.names(dif)
    library(dplyr)
    
    if(P=='p'){
      dif.up <- dif %>%
        filter(logFC > logFC_value & PValue < P_value)
      dif.down <- dif %>%
        filter(logFC < (-logFC_value) & PValue < P_value)
    }else{
      dif.up <- dif %>%
        filter(logFC > logFC_value & FDR < P_value)
      dif.down <- dif %>%
        filter(logFC < (-logFC_value) & FDR < P_value)
    }
    
    dif <- dif %>% mutate(group = case_when(
      gene_symbol %in% dif.up$gene_symbol ~ 'up',
      gene_symbol %in% dif.down$gene_symbol ~ 'down',
      TRUE ~ "no"
    ))
    
    DEG=data.frame(name=dif$gene_symbol,
                   logFC=dif$logFC,
                   P.value=dif$PValue,
                   P.adj=dif$FDR,
                   group=dif$group)
    
    ####输出数据####
    dir.create("C:/Users/18025/Desktop/EGFR/2.DEG/data", recursive = TRUE)
    dir.create("C:/Users/18025/Desktop/EGFR/2.DEG/Table", recursive = TRUE)
    print(table(DEG$group))
    
    DEG_list<-list(method="edgeR",
                   case=case,
                   control=control,
                   DEG=DEG,P=P,
                   logFC_value=logFC_value,P_value=P_value,
                   feature_type=feature_type,colors=colors)
    
    save(DEG_list,file = "C:/Users/18025/Desktop/EGFR/2.DEG/data/DEG_list.rda")
    write.csv(DEG,'C:/Users/18025/Desktop/EGFR/2.DEG/Table/DEG.csv')
    
    exp_plot=log2(cpm(DGElist)+1)
  }
  
  ####绘图####
  dir.create("C:/Users/18025/Desktop/EGFR/2.DEG/plott", recursive = TRUE)
  
  top10_padj = DEG[DEG$group!='no',] %>% group_by(group) %>% slice_max(n = 5, order_by = -log10(P.adj), with_ties = F)
  top10_p = DEG[DEG$group!='no',] %>% group_by(group) %>% slice_max(n=5, order_by = -log10(P.value), with_ties = F)
  
  ####火山图####
  library(ggplot2)
  if(P=="p") label = top10_p else label = top10_padj
  DEG$Label = ifelse(DEG$name %in% label$name, DEG$name, NA)
  Significant = factor(DEG$group, levels = c('up','down',"no"))
  dp <- ifelse(P=='p', 3, 4)
  
  library(ggrepel)
  p1 <- ggplot(DEG, aes(logFC, -log10(DEG[,dp]))) +
    geom_point(aes(col=Significant)) +
    scale_color_manual(values=c(colors[2], colors[1], "grey")) +
    labs(title = " ") +
    geom_vline(xintercept=c(-logFC_value,logFC_value), colour="black", linetype="dashed") +
    geom_hline(yintercept = -log10(0.05), colour="black", linetype="dashed") +
    labs(x="log2(FoldChange)", y="-log10(Pvalue)") +
    theme(plot.title = element_blank(), plot.subtitle = element_blank(),
          plot.background = element_blank(), plot.margin = margin(t=1,r=1,b=1,l=1,unit="pt"),
          panel.border = element_rect(size = 0.6, fill = NA, colour = "black"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(size=6, colour = "black"),
          axis.title = element_text(size=7),
          legend.background = element_blank(), legend.key = element_blank(),
          legend.title = element_text(size = 7), legend.text = element_text(size = 7),
          legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"), legend.box.spacing = unit(3,"pt")) +
    geom_text_repel(aes(label=Label), size=2.5)
  
  while (!is.null(dev.list())) dev.off()
  pdf(paste0('C:/Users/18025/Desktop/EGFR/2.DEG/plot/zhf',case,'_vs_',control,'_volplotzhf.pdf'), width = 8/2.54, height = 8/2.54)
  print(p1)
  dev.off()
  