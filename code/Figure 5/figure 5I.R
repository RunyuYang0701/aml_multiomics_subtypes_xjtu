library("DESeq2")
coldata=group
colnames(coldata)="condition"
counts=expr_gse132511
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData,
                              design = ~ condition)

colData <- data.frame(row.names = c("prim-AML-1-ROSlow","prim-AML-2-ROSlow","prim-AML-3-ROSlow","prim-AML-4-ROSlow","prim-AML-5-ROSlow","prim-AML-6-ROSlow","prim-AML-7-ROSlow","mono-AML-1-ROSlow","mono-AML-2-ROSlow","mono-AML-3-ROSlow","mono-AML-4-ROSlow","mono-AML-5-ROSlow"),
                      condition =
                        factor(c("prim","prim","prim","prim","prim","prim","prim","mono","mono","mono","mono","mono"),
                               levels = c("prim","mono")))
# 函数分析差异
dds <- DESeq(dds)
# 计算标准化因子
sizeFactors(dds)
#提取差异表达结果
res <- results(dds)
class(res)
res <- as.data.frame(res)
res=na.omit(res)

geneList <- res$log2FoldChange               # 获取GeneList
names(geneList) <-gene_set$entrez # 对GeneList命名
geneList <- sort(geneList, decreasing = T)  # 从高到低排序
head(geneList)
library(clusterProfiler)
#####富集分析######################
set.seed(123456)
geneSet_hall<- read.gmt("./h.all.v7.2.symbols.gmt")
geneSet_go<- read.gmt("./c5.go.v7.2.entrez.gmt")
egmt_hall<- GSEA(geneList, TERM2GENE=geneSet_go, verbose=FALSE, 
                     nPerm = 10000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
gsea_results <- egmt_hall@result
write.csv(gsea_results,file="gse_results.csv")

library(ggpubr)
gse_result$NES=-gse_result$NES
ggdotchart(gse_result, x = "Description", y = "NES", 
           color = "qvalue",
           sorting = "descending",
           add = "segments",
           rotate = TRUE,
           dot.size = 8,
           gradient.cols = c("#E71D36", "grey"),
           font.label = list(color = "white", size = 8, vjust = 0.5),              
           ggtheme = theme_pubr(),
           ylab = F)+
  geom_hline(yintercept = 0, linetype=2, linewidth=0.5)

##########亚型的biomarker富集结果###################
geneset=marker.up$templates
geneset=geneset[,c(2,1)]
colnames(geneset)=c("term","gene")
geneset$term=paste0(geneset$term,"suntypes")
egmt <- GSEA(geneList, TERM2GENE=geneset, verbose=FALSE, 
             nPerm = 10000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
gsea_results <- egmt@result

library(GseaVis)
gseaNb(object = egmt,
       geneSetID = 'CS3suntypes',
       newGsea = T,
       addPval = T)
geneSetID = c('CS3suntypes',
              'CS2suntypes')

# all plot
# sub plot
gseaNb(object = egmt,
       geneSetID = geneSetID,
       subPlot = 2,
       curveCol = c("#E71D36","#FF9F1C") )

