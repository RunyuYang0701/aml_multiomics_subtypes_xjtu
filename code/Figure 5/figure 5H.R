
###差异基因&富集分析
library(limma)
data=expr_gse132511
group_list <- factor(group[,1],ordered = F)
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(data)
fit <- lmFit(data,design)
cont.matrix <- makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2) # 贝叶斯检验
tmpOut <- topTable(fit2,coef = 1,n = Inf,adjust.method = "BH",sort.by = "logFC") 
limma.na <- na.omit(tmpOut)
dif <- limma.na[limma.na$P.Value <= 0.05 & abs(limma.na$logFC) > log2(2),]
colnames(expr_gse132511)
geneset=marker.up$templates
geneset=geneset[,c(1,2)]
colnames(geneset)=c("gene","term")
geneset$term=paste0(geneset$term,"suntypes")
geneList <- limma.na$logFC
names(geneList) <- rownames(limma.na)
geneList <- sort(geneList, decreasing = T)
geneList <- geneList[geneList != 0]
head(geneList)
geneset=geneset[,-1]
library(clusterProfiler)
set.seed(123456)
egmt <- GSEA(geneList, TERM2GENE=geneset, verbose=FALSE, 
             nPerm = 10000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
gsea_results <- egmt@result
gsea_results2 <- gsea_results[order(gsea_results$enrichmentScore,decreasing=F),]
head(geneset)



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

