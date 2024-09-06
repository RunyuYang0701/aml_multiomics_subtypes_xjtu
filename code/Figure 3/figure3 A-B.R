###数据读入
aml1012=read.table("GSM3587923_AML1012-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml1012_metadat=read.table("GSM3587924_AML1012-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml210=read.table("GSM3587925_AML210A-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml210_metadat=read.table("GSM3587926_AML210A-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml314=read.table("GSM3587927_AML314-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml314_metadat=read.table("GSM3587928_AML314-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml328=read.table("GSM3587931_AML328-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml328_metadat=read.table("GSM3587932_AML328-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml329=read.table("GSM3587940_AML329-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml329_metadat=read.table("GSM3587941_AML329-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml371=read.table("GSM3587946_AML371-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml371_metadat=read.table("GSM3587947_AML371-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml419=read.table("GSM3587950_AML419A-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml419_metadat=read.table("GSM3587951_AML419A-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml420=read.table("GSM3587953_AML420B-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml420_metadat=read.table("GSM3587954_AML420B-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml475=read.table("GSM3587959_AML475-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml475_metadat=read.table("GSM3587960_AML475-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml556=read.table("GSM3587963_AML556-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml556_metadat=read.table("GSM3587964_AML556-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml707=read.table("GSM3587969_AML707B-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml707_metadat=read.table("GSM3587970_AML707B-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml722=read.table("GSM3587980_AML722B-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml722_metadat=read.table("GSM3587981_AML722B-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml870=read.table("GSM3587984_AML870-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml870_metadat=read.table("GSM3587985_AML870-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml916=read.table("GSM3587988_AML916-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml916_metadat=read.table("GSM3587989_AML916-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml921=read.table("GSM3587990_AML921A-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml921_metadat=read.table("GSM3587991_AML921A-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
aml997=read.table("GSM3587992_AML997-D0.dem.txt", header=TRUE, sep="\t", row.names=1)
aml997_metadat=read.table("GSM3587993_AML997-D0.anno.txt", header=TRUE, sep="\t", row.names=1)
BM1=read.table("GSM3587996_BM1.dem.txt", header=TRUE, sep="\t", row.names=1)
BM1_metadat=read.table("GSM3587996_BM1.anno.txt", header=TRUE, sep="\t", row.names=1)
BM2=read.table("GSM3587997_BM2.dem.txt", header=TRUE, sep="\t", row.names=1)
BM2_metadat=read.table("GSM3587997_BM2.anno.txt", header=TRUE, sep="\t", row.names=1)
BM3=read.table("GSM3587998_BM3.dem.txt", header=TRUE, sep="\t", row.names=1)
BM3_metadat=read.table("GSM3587999_BM3.anno.txt", header=TRUE, sep="\t", row.names=1)
BM4=read.table("GSM3588000_BM4.dem.txt", header=TRUE, sep="\t", row.names=1)
BM4_metadat=read.table("GSM3588001_BM4.anno.txt", header=TRUE, sep="\t", row.names=1)
####s4文件构建
library(Seurat)
scaml210 = CreateSeuratObject(aml210,project="aml210")
aml210_metadat$orig.ident=scaml210@meta.data$orig.ident
scaml210@meta.data=aml210_metadat

scaml102 = CreateSeuratObject(aml1012,project="aml102")
aml1012_metadat$orig.ident=scaml102@meta.data$orig.ident
scaml102@meta.data=aml1012_metadat

scaml314 = CreateSeuratObject(aml314,project="aml314")
aml314_metadat$orig.ident=scaml314@meta.data$orig.ident
scaml314@meta.data=aml314_metadat

scaml328 = CreateSeuratObject(aml328,project="aml328")
aml328_metadat$orig.ident=scaml328@meta.data$orig.ident
scaml328@meta.data=aml328_metadat

scaml329 = CreateSeuratObject(aml329,project="aml329")
aml329_metadat$orig.ident=scaml329@meta.data$orig.ident
scaml329@meta.data=aml329_metadat

scaml371 = CreateSeuratObject(aml371,project="aml371")
aml371_metadat$orig.ident=scaml371@meta.data$orig.ident
scaml371@meta.data=aml371_metadat

scaml419 = CreateSeuratObject(aml419,project="aml419")
aml419_metadat$orig.ident=scaml419@meta.data$orig.ident
scaml419@meta.data=aml419_metadat

scaml420 = CreateSeuratObject(aml420,project="aml420")
aml420_metadat$orig.ident=scaml420@meta.data$orig.ident
scaml420@meta.data=aml420_metadat

scaml475 = CreateSeuratObject(aml475,project="aml475")
aml475_metadat$orig.ident=scaml475@meta.data$orig.ident
scaml475@meta.data=aml475_metadat

scaml556 = CreateSeuratObject(aml556,project="aml556")
aml556_metadat$orig.ident=scaml556@meta.data$orig.ident
scaml556@meta.data=aml556_metadat

scaml707 = CreateSeuratObject(aml707,project="aml707")
aml707_metadat$orig.ident=scaml707@meta.data$orig.ident
scaml707@meta.data=aml707_metadat

scaml722 = CreateSeuratObject(aml722,project="aml722")
aml722_metadat$orig.ident=scaml722@meta.data$orig.ident
scaml722@meta.data=aml722_metadat

scaml870 = CreateSeuratObject(aml870,project="aml870")
aml870_metadat$orig.ident=scaml870@meta.data$orig.ident
scaml870@meta.data=aml870_metadat

scaml916 = CreateSeuratObject(aml916,project="aml916")
aml916_metadat$orig.ident=scaml916@meta.data$orig.ident
scaml916@meta.data=aml916_metadat


scaml921 = CreateSeuratObject(aml921,project="aml921")
aml921_metadat$orig.ident=scaml921@meta.data$orig.ident
scaml921@meta.data=aml921_metadat

scaml997 = CreateSeuratObject(aml997,project="aml997")
aml997_metadat$orig.ident=scaml997@meta.data$orig.ident
scaml997@meta.data=aml997_metadat

scbm1=CreateSeuratObject(BM1,project="BM1")
aml997_metadat$orig.ident=scaml997@meta.data$orig.ident
scaml997@meta.data=aml997_metadat
########
library(Seurat)
########1.5对数据进行均一化，使用NormalizeData这个函数。
scRNA1 <- NormalizeData(scRNA_aml, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst", nfeatures = 3000) 
scale.genes <-  rownames(scRNA1)
scRNA1 <- ScaleData(scRNA1, features = scale.genes)
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1)) 
pc.num=1:30
scRNA1 <- FindNeighbors(scRNA1, dims = pc.num) 

scRNA1 <- FindClusters(scRNA1, resolution = 1.2)
scRNA1 <- RunUMAP(scRNA1, dims = pc.num)

plot2 = DimPlot(scRNA1, reduction = "umap",group.by = "celltype") 
 plot2
ggsave("UMAP.pdf", plot = plot2, width = 8, height = 7)
#合并tSNE与UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
plotc
ggsave("tSNE_UMAP.pdf", plot = plotc, width = 10, height = 5)
##保存数据这节课的数据
saveRDS(scRNA1, file="scRNA1.rds")
###
scRNA2=subset(scRNA1, subset = cel=="malignant")
plot2 = DimPlot(lymphoma, reduction = "umap",group.by = "celltype") 
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)

lymphoma <- NormalizeData(scRNA_aml) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
system.time({lymphoma <- RunHarmony(lymphoma, group.by.vars = "orig.ident")})
lymphoma <- FindNeighbors(lymphoma, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 1)
plot2 <- ElbowPlot(lymphoma, ndims=50, reduction="pca") 
lymphoma <- RunUMAP(lymphoma, reduction = "harmony", dims = 1:30)
plot2 = DimPlot(scRNA2, reduction = "umap")
scRNA2=subset(lymphoma, subset = cel=="malignant")
library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(Seurat)
library(SeuratData)
library(msigdbr)
library(patchwork)
cells_rankings <- AUCell_buildRankings(scRNA2@assays$RNA@data) 
geneset=marker.up$templates
genesets = split(geneset$probe,geneset$class)
cells_AUC <- AUCell_calcAUC(genesets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
geneSet <- "CS1"  
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])  #提取这个通路在每一个细胞的得分
scRNA2$AUC <- aucs  #将得分添加入scRNA（seruat）对象

#选择细胞展示的维度 
df<- data.frame(scRNA2@meta.data, scRNA2@reductions$umap@cell.embeddings)  #选择用UMAP维度看   也可以选择TSNE
head(df)
class_avg <- df %>%
  group_by(celltype) %>%        #这里可以改成cluster  seurat_clusters/或者其他的annotation
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

#通过ggplot画图
ggplot(df, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="A") +
  ggrepel::geom_label_repel(aes(label = celltype),
                            data = class_avg,
                            size = 4,
                            label.size = 0,
                            segment.color = NA
  )+   theme(legend.position = "none") + theme_bw()
aa=scRNA2@meta.data
write.csv(aa,file="aa.csv")
