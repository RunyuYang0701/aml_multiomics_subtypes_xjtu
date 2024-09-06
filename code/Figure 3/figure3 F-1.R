####拟时序分析 c3亚型的分析
##作者：yry
##时间：2024.6.24
library(Seurat)
library(monocle)
library(dplyr)
##自定义函数 
seurat_to_monocle <- function(otherCDS, assay, slot, lowerDetectionLimit = 0, import_all = FALSE) {
  if(class(otherCDS)[1] == 'Seurat') {
    requireNamespace("Seurat")
    data <- GetAssayData(otherCDS, assay = assay, slot = slot)
    data <- data[rowSums(as.matrix(data)) != 0,]
    pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    valid_data <- data[, row.names(pd)]
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
      } else {
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    } 
  } 
  return(monocle_cds)
}
##选取c3亚型的数据
scRNA_aml_cdc=subset(scRNA_aml,subset=celltype=="cDC")
scRNA_aml_cdc$clust="cdc"
scRNA_aml_mono=subset(scRNA_aml,subset=celltype=="Mono")
scRNA_aml_mono$clust="mono"

scRNA_aml1_mono=subset(scRNA_aml1,subset=celltype=="Mono-like")

scRNA_aml1_mono_c3=subset(scRNA_aml1_mono,subset=clust=="3")

scRNA_aml1_mono_c3$clust="mono-like-c3"

scRNA_aml1_cdc=subset(scRNA_aml1,subset=celltype=="cDC-like")

scRNA_aml1_cdc_c3=subset(scRNA_aml1_cdc,subset=clust=="3")


scRNA_aml1_cdc_c3$clust="cdc-like-c3"

scrna_c3=merge(scRNA_aml_cdc,c(scRNA_aml_mono,scRNA_aml1_mono_c3,scRNA_aml1_cdc_c3))
##monocle2 分析
aml_monocle_fun <- seurat_to_monocle(scrna_c3, assay = "RNA", slot = "counts")
aml_monocle <- estimateSizeFactors(aml_monocle_fun)
aml_monocle <- estimateDispersions(aml_monocle)
##数据质控
aml_monocle <- detectGenes(aml_monocle, min_expr = 0.1)#至少10%的表达，这一步完成之后，会出现num_cells_expressed这一列
print(head(fData(aml_monocle)))
print(head(pData(aml_monocle)))
expressed_genes <- row.names(subset(fData(aml_monocle), num_cells_expressed >= 10))
pData(aml_monocle)$Total_mRNAs <- Matrix::colSums(exprs(aml_monocle))#将UMI加入cds
aml_monocle <- aml_monocle[,pData(aml_monocle)$Total_mRNAs < 1e6]#阈值按照自己实际情况

#使用seurat计算差异基因
cds_seurat <- aml_monocle
seurat_variable_genes <- VariableFeatures(FindVariableFeatures(scrna_c3, assay = "RNA"), assay = "RNA")
cds_seurat <- setOrderingFilter(cds_seurat, seurat_variable_genes)
plot_ordering_genes(cds_DGT)
save(cds_seurat, file = 'cds_seurat.RData')

cds_seruat <- reduceDimension(cds_seurat, max_components = 2,reduction_method = 'DDRTree')
cds_seruat <- orderCells(cds_seruat)
plot_cell_trajectory(cds_seruat, color_by = "Pseudotime")
plot_cell_trajectory(cds_seruat, color_by = "orig.ident",cell_link_size = 1.0)
plot_cell_trajectory(cds_seruat, color_by = "celltype", cell_link_size = 1.0)
plot_cell_trajectory(cds_seruat, color_by = "State", cell_link_size = 1.0)

##可视化拟时序图
library(tibble)
#提取数据=======================================================================
data_df <- t(reducedDimS(cds_seruat)) %>% as.data.frame() %>% #提取坐标
  select_(Component_1 = 1, Component_2 = 2) %>% #重命名
  rownames_to_column("cells") %>% #rownames命名
  mutate(pData(cds_seruat)$State) %>% #添加State
  mutate(pData(cds_seruat)$Pseudotime, 
         pData(cds_seruat)$orig.ident, 
         pData(cds_seruat)$celltype,
         pData(cds_seruat)$cel)#将这些需要作图的有用信息都添加上

colnames(data_df) <- c("cells","Component_1","Component_2","State",
                       "Pseudotime","orig.ident","celltype","cel")

#==============================================================================
#轨迹数据提取---完全摘录于monocle包原函数
dp_mst <- minSpanningTree(cds_seruat)
reduced_dim_coords <- reducedDimK(cds_seruat)
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>% 
  mutate(sample_name = rownames(.), sample_state = rownames(.))


#构建一个做轨迹线图的数据
edge_df <- dp_mst %>% igraph::as_data_frame() %>% 
  select_(source = "from", target = "to") %>% 
  left_join(ica_space_df %>% select_(source = "sample_name", 
                                     source_prin_graph_dim_1 = "prin_graph_dim_1", 
                                     source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>% 
  left_join(ica_space_df %>% select_(target = "sample_name", 
                                     target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                     target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")




cols= c( '#11A579','#7F3C8D','#3969AC',
        '#E73F74', '#80BA5A', '#E68310'
)
library(ggplot2)
library(tidydr)
library(ggforce)
library(ggrastr)
library(viridis)
install.packages("tidydr")

g <- ggplot() + 
  geom_point_rast(data = data_df, aes(x = Component_1, 
                                      y = Component_2,
                                      color =State)) + #散点图
  scale_color_manual(values = alpha(cols,0.3))+ 
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               linewidth = 1, 
               linetype = "solid", na.rm = TRUE, data = edge_df)+#添加轨迹线
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


g <- ggplot() + 
  geom_point_rast(data = data_df, aes(x = Component_1, 
                                      y = Component_2,
                                      color =Pseudotime)) + #散点图
  scale_color_viridis()+ 
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               linewidth = 1, 
               linetype = "solid", na.rm = TRUE, data = edge_df)+#添加轨迹线
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

####饼图
library(plotrix)
cellratio <- as.data.frame(table(data_df$State,data_df$celltype))
clust1<- subset(cellratio, Var1=='1')
clust2<- subset(cellratio, Var1=='3')
clust3<- subset(cellratio, Var1=='4')
A = pie3D(x=clust1$Freq,
          radius=1,
          height=0.1,
          theta=pi/6,
          explode=0,
          main="BM",
          col=c( '#C5DCCD','#CDBAD4', '#C0C6DF',
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
      col=c( '#C5DCCD','#CDBAD4', '#C0C6DF',
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


A = pie3D(x=clust2$Freq,
          radius=1,
          height=0.1,
          theta=pi/6,
          explode=0,
          main="BM",
          col=c( '#C5DCCD','#CDBAD4', '#C0C6DF',
                 '#F3C7CD', '#D6E4C6', '#F2D3AF'),
          border = "black",
          shade = 0.5,
          labels=paste0(c(clust1$Var2),
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
      col=c( '#C5DCCD','#CDBAD4', '#C0C6DF',
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


A = pie3D(x=clust3$Freq,
          radius=1,
          height=0.1,
          theta=pi/6,
          explode=0,
          main="BM",
          col=c( '#C5DCCD','#CDBAD4', '#C0C6DF',
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

pie3D(x=clust3$Freq,
      labelpos=A,
      radius=1,
      height=0.1,
      theta=pi/6,
      explode=0.1,
      main="BM",
      col=c( '#C5DCCD','#CDBAD4', '#C0C6DF',
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


#############
library(plotrix)
cellratio <- as.data.frame(table(data_df$State,data_df$cel))
clust1<- subset(cellratio, Var1=='1')
clust2<- subset(cellratio, Var1=='3')
clust3<- subset(cellratio, Var1=='4')



A = pie3D(x=clust1$Freq,
          radius=1,
          height=0.1,
          theta=pi/6,
          explode=0,
          main="BM",
          col=c( '#67A9CC','#DA8A87',"grey"),
          border = "black",
          shade = 0.5,
          labels=paste0(c(clust1$Var2),
                        "\n",
                        round(clust1$Freq/sum(clust1$Freq) * 100,2), "%"),
          mar=c(2,2,2,3),
          labelcol = "black",
          labelcex = 0.8
)
A = pie3D(x=clust2$Freq,
          radius=1,
          height=0.1,
          theta=pi/6,
          explode=0,
          main="BM",
          col=c( '#67A9CC','#DA8A87',"grey"),
          border = "black",
          shade = 0.5,
          labels=paste0(c(clust2$Var2),
                        "\n",
                        round(clust2$Freq/sum(clust2$Freq) * 100,2), "%"),
          mar=c(2,2,2,3),
          labelcol = "black",
          labelcex = 0.8
)

A = pie3D(x=clust3$Freq,
          radius=1,
          height=0.1,
          theta=pi/6,
          explode=0,
          main="BM",
          col=c( '#67A9CC','#DA8A87',"grey"),
          border = "black",
          shade = 0.5,
          labels=paste0(c(clust3$Var2),
                        "\n",
                        round(clust3$Freq/sum(clust3$Freq) * 100,2), "%"),
          mar=c(2,2,2,3),
          labelcol = "black",
          labelcex = 0.8
)

cds_DGT_pseudotimegenes <- differentialGeneTest(cds_seruat,fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(cds_DGT_pseudotimegenes,file="cds_DGT_pseudotimegenes.csv")
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, qval < 1e-50)
Time_genes <- cds_DGT_pseudotimegenes_sig %>% pull(gene_short_name) %>% as.character()

p=plot_pseudotime_heatmap(cds_seruat[Time_genes,], num_clusters=4, cores=1, show_rownames=T)

p <-plot_pseudotime_heatmap(cds_seruat[Time_genes,],
                            num_clusters = 4,
                            cores = 2,
                            show_rownames = T,return_heatmap =T,
                            hmcols = viridis(256),
                            use_gene_short_name = T)
cds_subset <- cds_seruat
newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
                                       max(pData(cds_subset)$Pseudotime),
                                       length.out = 100)) 



m <- genSmoothCurves(cds_seruat[Time_genes,], 
                     trend_formula = '~sm.ns(Pseudotime, df=3)',  
                     relative_expr = T, new_data = newdata)
m=m[!apply(m,1,sum)==0,]
m <- log10(m+1) #log处理数据
#数据缩放
m=m[!apply(m,1,sd)==0,]
m=Matrix::t(scale(Matrix::t(m),center=TRUE))
m[m > 3] = 3#热图最大值
m[m <- 3] = -3#热图最小值，可自行调整



callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
library(pheatmap)
#先做一个普通热图
p1 <- pheatmap(m, 
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               clustering_method = "ward.D2",
               cutree_rows=4,
               filename=NA,
               border_color = NA,
               fontsize_row = 8,
               color=colorRampPalette(c("navy","white","firebrick3"))(100),
               clustering_callback = callback)


#列注释
annotation_col = data.frame(
  pseudotime = rescale(newdata$Pseudotime, to = c(-1, 1)))
row.names(annotation_col) <- colnames(m)


#行注释
annotation_row <- data.frame(Cluster=factor(cutree(p1$tree_row, 4)))
row.names(annotation_row) <- rownames(m)

rowcolor <- c( '#E55F55','#9FB5D7', '#80BB6B',
                   '#AF81B5', '#D6E4C6', '#F2D3AF') 
names(rowcolor) <- c("1","2","3","4") #类型颜色

#注释颜色修改
ann_colors <- list(pseudotime=viridis(100),
                   Cluster=rowcolor) #颜色设置


p2 <- pheatmap(m, 
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               clustering_method = "ward.D2",
               cutree_rows=4,
               filename=NA,
               border_color = NA,
               fontsize_row = 8,
               color=colorRampPalette(c("navy","white","firebrick3"))(100),
               annotation_col = annotation_col,
               annotation_colors=ann_colors,
               annotation_row = annotation_row,
               clustering_callback = callback,
               annotation_names_col = F,
               annotation_names_row = F,
               main="Pseudotime")




heat_gg <- m
heat_gg <- as.data.frame(heat_gg)
heat_gg <- heat_gg%>% mutate(gene=row.names(.)) %>% melt()#转化为ggplot画图需要的长列表

p3 <- ggplot(heat_gg,aes(x=variable,y=gene,fill=value))+
  geom_raster()+
  scale_fill_gradient2(low="#003366", high="#990033", mid="white",
                       name='Pseudotime')+
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = 'black',
                                   size = 8))+
  scale_y_discrete(position = "right")

#层次聚类
library(ggtree)
library(tidytree)
d <- m%>% as.data.frame()
hc <- hclust(dist(d), method = "ward.D2")
clus <-cutree(hc, 4)
d1 = data.frame(label=names(clus),member=factor(clus))
ph <- ggtree(as.phylo(hc))

#行注释
cluster <- d1 %>%
  ggplot(aes(x=1, y=label,fill=member))+
  geom_tile() + 
  scale_y_discrete(position="right") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank())+
  labs(fill = "cluster")+
  scale_fill_manual(values = c("#708090",'#68A180','#F3B1A0', '#D6E7A3'))


#列注释
group <- colnames(m) %>% as.data.frame() %>% 
  mutate(group=newdata$Pseudotime) %>%
  mutate(p="group") %>%
  ggplot(aes(.,y=1,fill=group))+
  geom_tile() + 
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text = element_blank(),
        panel.grid = element_blank())+
  scale_fill_gradientn(colours = c("#85B22E","#E29827",'#57C3F3',"#922927"))+
  scale_x_discrete(expand = c(0,0))+
  geom_segment(aes(x = 5, y = 1, xend = 95, yend = 1),
               arrow = arrow(length = unit(0.1, "inches"),
                             type = 'closed'))+
  theme(plot.margin = margin(0,0,0,0))
### gsea 富集分析 ####
library(clusterProfiler)
library(ggplot2)
source('./Monocle2_gene_enrichment.R')

module_gene <- as.data.frame(cutree(p1$tree_row,k=4))

module_gene$genes <- rownames(module_gene)
colnames(module_gene) <- c('module','gene')

module_gene4=module_gene[which(module_gene$module=="4"),]
geneList<-rownames(module_gene4) 
names(geneList) <- as.character(module_gene4[,2])
geneList<-rownames(d) 
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AnnotationHub)
library(dbplyr)
library(BiocFileCache)
geneList <- sort(geneList, decreasing = TRUE)

bitr(geneList, fromType="SYMBOL", toType=c("ENTREZID","UNIPROT"), OrgDb="org.Hs.eg.db")
Gene_ID <- bitr(geneList, fromType="SYMBOL", toType=c("ENTREZID","UNIPROT"), OrgDb="org.Hs.eg.db")
gene=Gene_ID$ENTREZID
ego <- enrichGO(gene = gene,
                OrgDb = org.Hs.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)


df <- data.frame(ego) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 200) %>%
  arrange(desc(pvalue))

ratio <- lapply(df$GeneRatio,function(x){as.numeric(eval(parse(text = x)))}) %>% unlist()
df$ratio <- ratio

df$Description <- factor(df$Description,levels = df$Description)

write.csv(df,file="df_4.csv")

module_gene3=module_gene[which(module_gene$module=="3"),]
geneList<-rownames(module_gene3) 
names(geneList) <- as.character(module_gene3[,2])
geneList<-rownames(d) 
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AnnotationHub)
library(dbplyr)
library(BiocFileCache)
geneList <- sort(geneList, decreasing = TRUE)

bitr(geneList, fromType="SYMBOL", toType=c("ENTREZID","UNIPROT"), OrgDb="org.Hs.eg.db")
Gene_ID <- bitr(geneList, fromType="SYMBOL", toType=c("ENTREZID","UNIPROT"), OrgDb="org.Hs.eg.db")
gene=Gene_ID$ENTREZID
ego <- enrichGO(gene = gene,
                OrgDb = org.Hs.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)

kk <- enrichKEGG(gene = gene,
                 organism = "hsa",
                 pvalueCutoff =0.05,
                 qvalueCutoff =1)


df <- data.frame(ego) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 100) %>%
  arrange(desc(pvalue))

ratio <- lapply(df$GeneRatio,function(x){as.numeric(eval(parse(text = x)))}) %>% unlist()
df$ratio <- ratio

df$Description <- factor(df$Description,levels = df$Description)

write.csv(df,file="df_3.csv")








module_gene2=module_gene[which(module_gene$module=="2"),]
geneList<-rownames(module_gene2) 
names(geneList) <- as.character(module_gene2[,2])
geneList<-rownames(d) 
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AnnotationHub)
library(dbplyr)
library(BiocFileCache)
geneList <- sort(geneList, decreasing = TRUE)

bitr(geneList, fromType="SYMBOL", toType=c("ENTREZID","UNIPROT"), OrgDb="org.Hs.eg.db")
Gene_ID <- bitr(geneList, fromType="SYMBOL", toType=c("ENTREZID","UNIPROT"), OrgDb="org.Hs.eg.db")
gene=Gene_ID$ENTREZID
ego <- enrichGO(gene = gene,
                OrgDb = org.Hs.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)

kk <- enrichKEGG(gene = gene,
                 organism = "hsa",
                 pvalueCutoff =0.05,
                 qvalueCutoff =1)


df <- data.frame(ego) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 100) %>%
  arrange(desc(pvalue))

ratio <- lapply(df$GeneRatio,function(x){as.numeric(eval(parse(text = x)))}) %>% unlist()
df$ratio <- ratio

df$Description <- factor(df$Description,levels = df$Description)

write.csv(df,file="df_2.csv")







module_gene1=module_gene[which(module_gene$module=="1"),]
geneList<-rownames(module_gene1) 
names(geneList) <- as.character(module_gene1[,2])
geneList<-rownames(d) 
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AnnotationHub)
library(dbplyr)
library(BiocFileCache)
geneList <- sort(geneList, decreasing = TRUE)

bitr(geneList, fromType="SYMBOL", toType=c("ENTREZID","UNIPROT"), OrgDb="org.Hs.eg.db")
Gene_ID <- bitr(geneList, fromType="SYMBOL", toType=c("ENTREZID","UNIPROT"), OrgDb="org.Hs.eg.db")
gene=Gene_ID$ENTREZID
ego <- enrichGO(gene = gene,
                OrgDb = org.Hs.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)




df <- data.frame(ego) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 200) %>%
  arrange(desc(pvalue))

ratio <- lapply(df$GeneRatio,function(x){as.numeric(eval(parse(text = x)))}) %>% unlist()
df$ratio <- ratio

df$Description <- factor(df$Description,levels = df$Description)
write.csv(df,file="df_1.csv")

df=kkkk
ggplot(df) +
  ggforce::geom_link(aes(x = 0,y = Description,
                         xend = -log10(pvalue),yend = Description,
                         alpha = stat(index),
                         color = ONTOLOGY,
                         size = after_stat(index)),
                     n = 500,
                     # color = "#FF0033",
                     show.legend = F) +
  geom_point(aes(x = -log10(pvalue),y = Description),
             color = "black",
             fill = "white",size = 6,shape = 21) +
  geom_line(aes(x = ratio*100,y = Description,group = 1),
            orientation = "y",linewidth = 1,color = "#FFCC00") +
  scale_x_continuous(sec.axis = sec_axis(~./100,
                                         labels = scales::label_percent(),
                                         name = "Percent of geneRatio")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        axis.text = element_text(color = "black")) +
  ylab("") + xlab("-log10 Pvalue") +
  facet_wrap(~ONTOLOGY,scales = "free",ncol = 1) +
  scale_color_brewer(palette = "Set1")
pdf("plot.pdf", width = 5, height = 10)
dev.off()

library(monocle)
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
library(dplyr)
library(tidytree)
library(viridis)
genes <- c("USP15")
genes_exp <- list()
for(i in 1:length(genes)){
  A <- log2(exprs(cds_seruat)[genes[i],]+1)
  A <- as.data.frame(A)
  genes_exp[[i]] <- A
}
gene_exp <- do.call(cbind, genes_exp)
colnames(gene_exp) <- genes
#将上述几个基因的拟时表达添加到monocle
pData(cds_seruat) = cbind(pData(cds_seruat), gene_exp)
#提取作图数据，只需要游基因表达和拟时即可
data <- pData(cds_seruat)
colnames(data)
#选择需要的列即可，我这里的origin.ident就是分组
data<-data %>% select("celltype","Pseudotime", "USP15" )
features <- c("USP15")
plist <- list()
#ggplot作图
for (i in 1:length(features)){
  df <- data[, colnames(data)%in% c("celltype","Pseudotime",features[i])]
  colnames(df) <- c("celltype","Pseudotime",'gene')
  p <- ggplot(df, aes(x=Pseudotime,
                      y=gene)) +
    geom_point(aes(color=celltype), size=0.8)+ scale_color_manual(values = alpha(cols,0.3))+
    geom_smooth(method = "loess",level = 0.95,
                formula = y~x, color='black',se=F)+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black",size = 1),
          axis.title = element_blank(),
          axis.text = element_text(size = 15, colour = 'black'),
          legend.position = 'none')+
    ggtitle(features[i])
  plist[[i]] <- p
}

####数据保存
save(scRNA_aml, file = 'scrna_aml_all.RData')
save(scRNA_aml1, file = 'scrna_aml_ma.RData')
save(scrna_c3,file="scrna_aml_c3.RData")
