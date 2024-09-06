######基于bulk的思想判定3个亚群的分化
##作者：yry
##时间：2024.7.4

scaml_beat = CreateSeuratObject(expr_beat,project="beat")
scaml_beat@meta.data$clust=beataml.ntp.pred$clust.res$clust

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




##monocle2 分析
aml_monocle_fun <- seurat_to_monocle(scaml_beat, assay = "RNA", slot = "counts")
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
seurat_variable_genes <- VariableFeatures(FindVariableFeatures(scaml_beat, assay = "RNA"), assay = "RNA")
cds_seurat <- setOrderingFilter(cds_seurat, seurat_variable_genes)
plot_ordering_genes(cds_DGT)
save(cds_seurat, file = 'cds_seurat_bulk_aml.RData')

cds_seruat <- reduceDimension(cds_seurat, max_components = 2,reduction_method = 'DDRTree',norm_method= "none")
cds_seruat <- orderCells(cds_seruat,root_state = 7)
plot_cell_trajectory(cds_seruat, color_by = "Pseudotime")
plot_cell_trajectory(cds_seruat, color_by = "clust",cell_link_size = 1.0)
plot_cell_trajectory(cds_seruat, color_by = "celltype", cell_link_size = 1.0)
plot_cell_trajectory(cds_seruat, color_by = "State", cell_link_size = 1.0)

####绘图
data_df <- t(reducedDimS(cds_seruat)) %>% as.data.frame() %>% #提取坐标
  select_(Component_1 = 1, Component_2 = 2) %>% #重命名
  rownames_to_column("cells") %>% #rownames命名
  mutate(pData(cds_seruat)$State) %>% #添加State
  mutate(pData(cds_seruat)$Pseudotime, 
         pData(cds_seruat)$orig.ident, 
         pData(cds_seruat)$clust
  )#将这些需要作图的有用信息都添加上

colnames(data_df) <- c("cells","Component_1","Component_2","State",
                       "Pseudotime","orig.ident","clust")

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




cols= c( "#2EC4B6","#E71D36","#FF9F1C"
)
library(ggplot2)
library(tidydr)
library(ggforce)
library(ggrastr)
library(viridis)


g <- ggplot() + 
  geom_point_rast(data = data_df, aes(x = Component_1, 
                                      y = Component_2,
                                      color =clust)) + #散点图
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


#####组成饼图

library(plotrix)
cellratio <- as.data.frame(table(data_df$State,data_df$clust))
clust3<- subset(cellratio, Var1=='1')

A = pie3D(x=clust3$Freq,
          radius=1,
          height=0.1,
          theta=pi/6,
          explode=0,
          main="BM",
          col=c( "#2EC4B6","#E71D36","#FF9F1C"),
          border = "black",
          shade = 0.5,
          labels=paste0(c(clust3$Var2),
                        "\n",
                        round(clust3$Freq/sum(clust3$Freq) * 100,2), "%"),
          mar=c(2,2,2,3),
          labelcol = "black",
          labelcex = 0.8
)

clust1_3<- subset(cellratio, Var1=='3')

A = pie3D(x=clust1_3$Freq,
          radius=1,
          height=0.1,
          theta=pi/6,
          explode=0,
          main="BM",
          col=c( "#2EC4B6","#E71D36","#FF9F1C"),
          border = "black",
          shade = 0.5,
          labels=paste0(c(clust1_3$Var2),
                        "\n",
                        round(clust1_3$Freq/sum(clust1_3$Freq) * 100,2), "%"),
          mar=c(2,2,2,3),
          labelcol = "black",
          labelcex = 0.8
)

clust1_6<- subset(cellratio, Var1=='6')
A = pie3D(x=clust1_6$Freq,
          radius=1,
          height=0.1,
          theta=pi/6,
          explode=0,
          main="BM",
          col=c( "#2EC4B6","#E71D36","#FF9F1C"),
          border = "black",
          shade = 0.5,
          labels=paste0(c(clust1_6$Var2),
                        "\n",
                        round(clust1_6$Freq/sum(clust1_6$Freq) * 100,2), "%"),
          mar=c(2,2,2,3),
          labelcol = "black",
          labelcex = 0.8
)




clust2_8<- subset(cellratio, Var1=='8')

clust2_7$Freq=c(50,107,7)
A = pie3D(x=clust2_7$Freq,
          radius=1,
          height=0.1,
          theta=pi/6,
          explode=0,
          main="BM",
          col=c( "#2EC4B6","#E71D36","#FF9F1C"),
          border = "black",
          shade = 0.5,
          labels=paste0(c(clust2_7$Var2),
                        "\n",
                        round(clust2_7$Freq/sum(clust2_7$Freq) * 100,2), "%"),
          mar=c(2,2,2,3),
          labelcol = "black",
          labelcex = 0.8
)

