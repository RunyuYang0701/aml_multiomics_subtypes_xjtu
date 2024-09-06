###单细胞亚组分型及预测
##作者：yry
##时间：2024.6.29


library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)






scRNA_aml1=subset(scRNA_aml,subset=cel=="malignant")
lymphoma <- NormalizeData(scRNA_aml1) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
system.time({lymphoma <- RunHarmony(lymphoma, group.by.vars = "orig.ident")})
lymphoma <- FindNeighbors(lymphoma, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 2)
lymphoma <- RunUMAP(lymphoma, reduction = "harmony", dims = 1:30)
scRNA_aml1=lymphoma
DimPlot(scRNA_aml1, reduction = "umap",group.by = "celltype")

library(MOVICS)
GSE1236.ntp.pred <- runNTP(expr       = exprSet,
                             templates  = marker.up$templates, # the template has been already prepared in runMarker()
                             scaleFlag  = TRUE, # scale input data (by default)
                             centerFlag = TRUE, # center input data (by default)
                             doPlot     = TRUE, # to generate heatmap
                             fig.name   = "NTP HEATMAP FOR 1236") 
####绘图#####
library(Seurat)
library(ggplot2)
library(ggrastr)
library(tidydr)
library(dplyr)
library(ggrepel)

adj_scRNA=scRNA_aml1
df <- adj_scRNA@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(cell_type = adj_scRNA@meta.data$celltype)

label <- df %>%group_by(cell_type) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))%>%
  as.data.frame()
rownames(label) <- label$cell_type

label$number <- seq(1:6)

cols= c('#7F3C8D' ,'#11A579', '#3969AC',
        '#E73F74', '#80BA5A', '#E68310'
        )

#ggplot作图
p = ggplot()+
  geom_point_rast(data=df, aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type),size = 1,shape=16) +
  scale_color_manual(values = alpha(cols,0.3))+ #设置下透明度
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none')+
  geom_point(data = label, aes(x= UMAP_1 , y = UMAP_2), size=6, color='white', alpha=0.9)+
  geom_point(data = label, aes(x= UMAP_1 , y = UMAP_2), size=6, color='black', shape=21)+
  geom_text(data = label,
            mapping = aes(x= UMAP_1 , y = UMAP_2, label = number),
            color='black')
library(AUCell)
library(clusterProfiler)
cells_rankings <- AUCell_buildRankings(scRNA_aml1@assays$RNA@data) 
geneset=marker.up$templates
genesets = split(geneset$probe,geneset$class)
cells_AUC <- AUCell_calcAUC(genesets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

geneSet <- "CS1"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
scRNA_aml1$AUC  <- aucs

library(ggraph)
ggplot(data.frame(scRNA_aml1@meta.data, scRNA_aml1@reductions$umap@cell.embeddings), aes(UMAP_1, UMAP_2, color=AUC)
) + geom_point( size=1.5
) + scale_color_viridis(option="A")  + theme_light(base_size = 15)+labs(title = "TNFA_SIGNALING_VIA_NFKB")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))

######clust1
library(plotrix)
cellratio <- as.data.frame(table(scRNA_aml1$clust, scRNA_aml1$celltype))
clust1<- subset(cellratio, Var1=='1')
clust2<- subset(cellratio, Var1=='2')
clust3<- subset(cellratio, Var1=='3')
A = pie3D(x=clust1$Freq,
          radius=1,
          height=0.1,
          theta=pi/6,
          explode=0,
          main="BM",
          col=c('#CDBAD4' ,'#C5DCCD', '#C0C6DF',
                '#F3C7CD', '#D6E4C6', '#F2D3AF'),
          border = "black",
          shade = 0.5,
          labels=paste0(c(clust1$Var2),
                        "\n",
                        round(clust1$Freq/sum(clust1$Freq) * 100,2), "%"),
          mar=c(2,2,2,3),
          labelcol = "black",
          labelcex = 0.8
)

pie3D(x=clust1$Freq,
      labelpos=A,
      radius=1,
      height=0.1,
      theta=pi/6,
      explode=0.1,
      main="BM",
      col=c('#CDBAD4' ,'#C5DCCD', '#C0C6DF',
            '#F3C7CD', '#D6E4C6', '#F2D3AF'),
      border = "black",
      shade = 0.5,
      labels=paste0(c(clust1$Var2),
                    "\n",
                    round(clust1$Freq/sum(clust1$Freq) * 100,2), "%"),
      mar=c(2,2,2,3),
      labelcol = "black",
      labelcex = 0.8
)
#######clust2
library(plotrix)
cellratio <- as.data.frame(table(scRNA_aml1$clust, scRNA_aml1$celltype))
clust1<- subset(cellratio, Var1=='1')
clust2<- subset(cellratio, Var1=='2')
clust3<- subset(cellratio, Var1=='3')
A = pie3D(x=clust2$Freq,
          radius=1,
          height=0.1,
          theta=pi/6,
          explode=0,
          main="BM",
          col=c('#7F3C8D' ,'#11A579', '#3969AC',
                '#E73F74', '#80BA5A', '#E68310'),
          border = "black",
          shade = 0.5,
          labels=paste0(c(clust2$Var2),
                        "\n",
                        round(clust2$Freq/sum(clust2$Freq) * 100,2), "%"),
          mar=c(2,2,2,3),
          labelcol = "black",
          labelcex = 0.8
)

pie3D(x=clust2$Freq,
      labelpos=A,
      radius=1,
      height=0.1,
      theta=pi/6,
      explode=0.1,
      main="BM",
      col=c('#CDBAD4' ,'#C5DCCD', '#C0C6DF',
            '#F3C7CD', '#D6E4C6', '#F2D3AF'),
      border = "black",
      shade = 0.5,
      labels=paste0(c(clust2$Var2),
                    "\n",
                    round(clust2$Freq/sum(clust2$Freq) * 100,2), "%"),
      mar=c(2,2,2,3),
      labelcol = "black",
      labelcex = 0.8
)
#######clust3
library(plotrix)
cellratio <- as.data.frame(table(scRNA_aml1$clust, scRNA_aml1$celltype))
clust1<- subset(cellratio, Var1=='1')
clust2<- subset(cellratio, Var1=='2')
clust3<- subset(cellratio, Var1=='3')
A = pie3D(x=clust3$Freq,
          radius=1,
          height=0.1,
          theta=pi/6,
          explode=0,
          main="BM",
          col=c('#7F3C8D' ,'#11A579', '#3969AC',
                '#E73F74', '#80BA5A', '#E68310'),
          border = "black",
          shade = 0.5,
          labels=paste0(c(clust3$Var2),
                        "\n",
                        round(clust3$Freq/sum(clust3$Freq) * 100,2), "%"),
          mar=c(2,2,2,3),
          labelcol = "black",
          labelcex = 0.8
)
pie3D(x=clust3$Freq,
      labelpos=A,
      radius=1,
      height=0.1,
      theta=pi/6,
      explode=0.1,
      main="BM",
      col=c('#CDBAD4' ,'#C5DCCD', '#C0C6DF',
            '#F3C7CD', '#D6E4C6', '#F2D3AF'),
      border = "black",
      shade = 0.5,
      labels=paste0(c(clust3$Var2),
                    "\n",
                    round(clust3$Freq/sum(clust3$Freq) * 100,2), "%"),
      mar=c(2,2,2,3),
      labelcol = "black",
      labelcex = 0.8
)
