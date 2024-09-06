#######figure 4A

####GSE71014
expr=exprSet[consam,]
expr1 <- apply(expr, 2, as.numeric)
expr1=as.data.frame(expr1)
rownames(expr1)=rownames(expr)
expr=expr1
expr=as.data.frame(t(expr))
subtype=target[,c(1,2)]
aas=merge(expr,subtype,by.x=0,by.y=1)
rownames(aas)=aas$Row.names
aas=aas[,-1]
aas$clust=paste0("clust",aas$clust)

expMat <- as.data.frame(t(apply(aas[,setdiff(colnames(aas), "clust")], 2, 
                                function(x) 
                                  tapply(x, 
                                         INDEX = factor(aas$clust), 
                                         FUN = median, 
                                         na.rm = TRUE)))) # 对同一亚型内的样本取中位数


library(ComplexHeatmap)
library(circlize)
library(ChAMPdata)
library(data.table)
??ChAMPdata
###
# 设置热图颜色
heatmap.BlWtRd <- c("#61AACF","#98CADD","#EAEFF6","#F9EFEF","#E9C6C6","#DA9599")

subt=c("clust1","clust2","clust3")
annCol <- data.frame(subtype = subt,
                     row.names = subt)
annCol <- annCol[order(annCol$subtype),,drop = F] # 按照亚型排序
annColors <- list()
annColors[["subtype"]] <- c("clust1" = "#2EC4B6",
                            "clust2" = "#E71D36",
                            "clust3" = "#FF9F1C")
top_anno <- HeatmapAnnotation(df                   = annCol,
                              col                  = annColors,
                              gp                   = gpar(col = "grey80"), # 每个单元格边框为灰色
                              simple_anno_size     = unit(3.5, "mm"), # 注释高3.5毫米
                              show_legend          = F, # 不显示亚型的图例，因为一目了然
                              show_annotation_name = F, # 不显示该注释的名称
                              border               = FALSE) # 不显示注释的外边框
rownames(immunomodulator)=immunomodulator$Gene
immunomodulator=immunomodulator[consam,]
# 创建行注释
annRow <- immunomodulator
annRow[which(annRow$Category == "Co-stimulator"),"Category"] <- "Co-stm" # 这里字符少一些，不会挤在一起，可以后期AI
annRow[which(annRow$Category == "Co-inhibitor"),"Category"] <- "Co-ihb"
annRow[which(annRow$Category == "Cell adhesion"),"Category"] <- "Cell\nadhesion" # 这里换行，不会挤在一起，可以后期AI
annRow[which(annRow$Category == "Antigen presentation"),"Category"] <- "Antigen\npresentation"
annRow$Category <- factor(annRow$Category, levels = c("Co-stm","Co-ihb","Ligand","Receptor","Cell\nadhesion","Antigen\npresentation","Other")) # 由于行需要按照类分割，所以需要定义因子顺序，否则按照字母表
annRow$ICI <- factor(annRow$ICI, levels = c("Inhibitory","N/A","Stimulatory"))
annRowColors <- list("ICI" = c("Inhibitory" = "#3D5C6F","N/A" = "#888888","Stimulatory" = "#E47159"))
left_anno <- HeatmapAnnotation(df                   = annRow[,"ICI",drop = F],
                               which                = "row", # 这里是行注释（默认为列）
                               gp                   = gpar(col = "grey80"), # 每个单元格边框为灰色
                               col                  = annRowColors,
                               simple_anno_size     = unit(3.5, "mm"), # 注释宽3.5毫米
                               show_annotation_name = F,
                               border               = F)


## 绘制表达谱热图（参数下同）

row_normalize <- function(df) {
  normalized_df <- t(apply(df, 1, function(x) {
    min_val <- min(x)
    max_val <- max(x)
    normalized_row <- (x - min_val) / (max_val - min_val) * 2 - 1
    return(normalized_row)
  }))
  colnames(normalized_df) <- colnames(df)
  return(normalized_df)
}

# 调用按行标准化函数
normalized_data <- row_normalize(expMat)

# 打印标准化后的数据
print(normalized_data)




col_expr <- colorRamp2(seq(-1, 1, length = 6), heatmap.BlWtRd) # 创建热图颜色（将热图输入矩阵的最大最小值取5个点，分配颜色红蓝色板；注意矩阵中可能存在的NA值）
hm.expr <- Heatmap(matrix             = as.matrix(normalized_data),
                   col                = col_expr,
                   border             = NA, # 无热图外边框
                   rect_gp = gpar(col = "grey80"), # 热图单元格边框为灰色
                   cluster_rows       = F, # 行不聚类
                   cluster_columns    = F, # 列不聚类
                   show_row_names     = T, # 显示行名
                   row_names_side     = "left", # 行名显示在左侧
                   row_names_gp       = gpar(fontsize = 10), # 行名字号为10
                   show_column_names  = F, # 不显示列名（可后期在颜色内AI使得亚型一目了然）
                   column_names_side  = "top", # 列名显示在顶部
                   row_split          = annRow$Category, # 行按照Category进行分割（因子顺序）
                   top_annotation     = top_anno, # 热图顶部注释
                   left_annotation    = left_anno, # 热图左侧注释
                   name               = "mRNA\nExpression", # 热图颜色图例的名称
                   width              = ncol(expMat) * unit(3.5, "mm"), # 热图单元格宽度（稍大于高度，因为所有注释都放在底部，平衡图形纵横比）
                   height             = nrow(expMat) * unit(2.4, "mm")) # 热图单元格高度

pdf(file = "complexheatmap of immunomodulator.pdf", width = 8,height = 12)
invisible(dev.off())
########GSE10358
expr=exprSet[consam,]
expr1 <- apply(expr, 2, as.numeric)
expr1=as.data.frame(expr1)
rownames(expr1)=rownames(expr)
expr=expr1
expr=na.omit(expr)
expr=as.data.frame(t(expr))
subtype=GSE10358.ntp.pred$clust.res
aas=merge(expr,subtype,by.x=0,by.y=1)
rownames(aas)=aas$Row.names
aas=aas[,-1]
aas$clust=paste0("clust",aas$clust)

expMat <- as.data.frame(t(apply(aas[,setdiff(colnames(aas), "clust")], 2, 
                                function(x) 
                                  tapply(x, 
                                         INDEX = factor(aas$clust), 
                                         FUN = median, 
                                         na.rm = TRUE)))) # 对同一亚型内的样本取中位数

library(ComplexHeatmap)
library(circlize)
library(ChAMPdata)
library(data.table)
??ChAMPdata
###
# 设置热图颜色
heatmap.BlWtRd <- c("#61AACF","#98CADD","#EAEFF6","#F9EFEF","#E9C6C6","#DA9599")

subt=c("clust1","clust2","clust3")
annCol <- data.frame(subtype = subt,
                     row.names = subt)
annCol <- annCol[order(annCol$subtype),,drop = F] # 按照亚型排序
annColors <- list()
annColors[["subtype"]] <- c("clust1" = "#2EC4B6",
                            "clust2" = "#E71D36",
                            "clust3" = "#FF9F1C")
top_anno <- HeatmapAnnotation(df                   = annCol,
                              col                  = annColors,
                              gp                   = gpar(col = "grey80"), # 每个单元格边框为灰色
                              simple_anno_size     = unit(3.5, "mm"), # 注释高3.5毫米
                              show_legend          = F, # 不显示亚型的图例，因为一目了然
                              show_annotation_name = F, # 不显示该注释的名称
                              border               = FALSE) # 不显示注释的外边框
rownames(immunomodulator)=immunomodulator$Gene
immunomodulator=immunomodulator[consam,]
# 创建行注释
annRow <- immunomodulator
annRow[which(annRow$Category == "Co-stimulator"),"Category"] <- "Co-stm" # 这里字符少一些，不会挤在一起，可以后期AI
annRow[which(annRow$Category == "Co-inhibitor"),"Category"] <- "Co-ihb"
annRow[which(annRow$Category == "Cell adhesion"),"Category"] <- "Cell\nadhesion" # 这里换行，不会挤在一起，可以后期AI
annRow[which(annRow$Category == "Antigen presentation"),"Category"] <- "Antigen\npresentation"
annRow$Category <- factor(annRow$Category, levels = c("Co-stm","Co-ihb","Ligand","Receptor","Cell\nadhesion","Antigen\npresentation","Other")) # 由于行需要按照类分割，所以需要定义因子顺序，否则按照字母表
annRow$ICI <- factor(annRow$ICI, levels = c("Inhibitory","N/A","Stimulatory"))
annRowColors <- list("ICI" = c("Inhibitory" = "#3D5C6F","N/A" = "#888888","Stimulatory" = "#E47159"))
left_anno <- HeatmapAnnotation(df                   = annRow[,"ICI",drop = F],
                               which                = "row", # 这里是行注释（默认为列）
                               gp                   = gpar(col = "grey80"), # 每个单元格边框为灰色
                               col                  = annRowColors,
                               simple_anno_size     = unit(3.5, "mm"), # 注释宽3.5毫米
                               show_annotation_name = F,
                               border               = F)


## 绘制表达谱热图（参数下同）

row_normalize <- function(df) {
  normalized_df <- t(apply(df, 1, function(x) {
    min_val <- min(x)
    max_val <- max(x)
    normalized_row <- (x - min_val) / (max_val - min_val) * 2 - 1
    return(normalized_row)
  }))
  colnames(normalized_df) <- colnames(df)
  return(normalized_df)
}

# 调用按行标准化函数
normalized_data <- row_normalize(expMat)

# 打印标准化后的数据
print(normalized_data)




col_expr <- colorRamp2(seq(-1, 1, length = 6), heatmap.BlWtRd) # 创建热图颜色（将热图输入矩阵的最大最小值取5个点，分配颜色红蓝色板；注意矩阵中可能存在的NA值）
hm.expr <- Heatmap(matrix             = as.matrix(normalized_data),
                   col                = col_expr,
                   border             = NA, # 无热图外边框
                   rect_gp = gpar(col = "grey80"), # 热图单元格边框为灰色
                   cluster_rows       = F, # 行不聚类
                   cluster_columns    = F, # 列不聚类
                   show_row_names     = T, # 显示行名
                   row_names_side     = "left", # 行名显示在左侧
                   row_names_gp       = gpar(fontsize = 10), # 行名字号为10
                   show_column_names  = F, # 不显示列名（可后期在颜色内AI使得亚型一目了然）
                   column_names_side  = "top", # 列名显示在顶部
                   row_split          = annRow$Category, # 行按照Category进行分割（因子顺序）
                   top_annotation     = top_anno, # 热图顶部注释
                   left_annotation    = left_anno, # 热图左侧注释
                   name               = "mRNA\nExpression", # 热图颜色图例的名称
                   width              = ncol(expMat) * unit(3.5, "mm"), # 热图单元格宽度（稍大于高度，因为所有注释都放在底部，平衡图形纵横比）
                   height             = nrow(expMat) * unit(2.4, "mm")) # 热图单元格高度

pdf(file = "complexheatmap of immunomodulator.pdf", width = 8,height = 12)
invisible(dev.off())

####BEAT_AML
expr=expr_beat[consam,]
expr1 <- apply(expr, 2, as.numeric)
expr1=as.data.frame(expr1)
rownames(expr1)=rownames(expr)
expr=expr1
expr=na.omit(expr)
expr=as.data.frame(t(expr))
subtype=beataml.ntp.pred$clust.res
aas=merge(expr,subtype,by.x=0,by.y=1)
rownames(aas)=aas$Row.names
aas=aas[,-1]
aas$clust=paste0("clust",aas$clust)

expMat <- as.data.frame(t(apply(aas[,setdiff(colnames(aas), "clust")], 2, 
                                function(x) 
                                  tapply(x, 
                                         INDEX = factor(aas$clust), 
                                         FUN = median, 
                                         na.rm = TRUE)))) # 对同一亚型内的样本取中位数

library(ComplexHeatmap)
library(circlize)
library(ChAMPdata)
library(data.table)
??ChAMPdata
###
# 设置热图颜色
heatmap.BlWtRd <- c("#61AACF","#98CADD","#EAEFF6","#F9EFEF","#E9C6C6","#DA9599")

subt=c("clust1","clust2","clust3")
annCol <- data.frame(subtype = subt,
                     row.names = subt)
annCol <- annCol[order(annCol$subtype),,drop = F] # 按照亚型排序
annColors <- list()
annColors[["subtype"]] <- c("clust1" = "#2EC4B6",
                            "clust2" = "#E71D36",
                            "clust3" = "#FF9F1C")
top_anno <- HeatmapAnnotation(df                   = annCol,
                              col                  = annColors,
                              gp                   = gpar(col = "grey80"), # 每个单元格边框为灰色
                              simple_anno_size     = unit(3.5, "mm"), # 注释高3.5毫米
                              show_legend          = F, # 不显示亚型的图例，因为一目了然
                              show_annotation_name = F, # 不显示该注释的名称
                              border               = FALSE) # 不显示注释的外边框
rownames(immunomodulator)=immunomodulator$Gene
immunomodulator=immunomodulator[consam,]
# 创建行注释
annRow <- immunomodulator
annRow[which(annRow$Category == "Co-stimulator"),"Category"] <- "Co-stm" # 这里字符少一些，不会挤在一起，可以后期AI
annRow[which(annRow$Category == "Co-inhibitor"),"Category"] <- "Co-ihb"
annRow[which(annRow$Category == "Cell adhesion"),"Category"] <- "Cell\nadhesion" # 这里换行，不会挤在一起，可以后期AI
annRow[which(annRow$Category == "Antigen presentation"),"Category"] <- "Antigen\npresentation"
annRow$Category <- factor(annRow$Category, levels = c("Co-stm","Co-ihb","Ligand","Receptor","Cell\nadhesion","Antigen\npresentation","Other")) # 由于行需要按照类分割，所以需要定义因子顺序，否则按照字母表
annRow$ICI <- factor(annRow$ICI, levels = c("Inhibitory","N/A","Stimulatory"))
annRowColors <- list("ICI" = c("Inhibitory" = "#3D5C6F","N/A" = "#888888","Stimulatory" = "#E47159"))
left_anno <- HeatmapAnnotation(df                   = annRow[,"ICI",drop = F],
                               which                = "row", # 这里是行注释（默认为列）
                               gp                   = gpar(col = "grey80"), # 每个单元格边框为灰色
                               col                  = annRowColors,
                               simple_anno_size     = unit(3.5, "mm"), # 注释宽3.5毫米
                               show_annotation_name = F,
                               border               = F)


## 绘制表达谱热图（参数下同）

row_normalize <- function(df) {
  normalized_df <- t(apply(df, 1, function(x) {
    min_val <- min(x)
    max_val <- max(x)
    normalized_row <- (x - min_val) / (max_val - min_val) * 2 - 1
    return(normalized_row)
  }))
  colnames(normalized_df) <- colnames(df)
  return(normalized_df)
}

# 调用按行标准化函数
normalized_data <- row_normalize(expMat)

# 打印标准化后的数据
print(normalized_data)




col_expr <- colorRamp2(seq(-1, 1, length = 6), heatmap.BlWtRd) # 创建热图颜色（将热图输入矩阵的最大最小值取5个点，分配颜色红蓝色板；注意矩阵中可能存在的NA值）
hm.expr <- Heatmap(matrix             = as.matrix(normalized_data),
                   col                = col_expr,
                   border             = NA, # 无热图外边框
                   rect_gp = gpar(col = "grey80"), # 热图单元格边框为灰色
                   cluster_rows       = F, # 行不聚类
                   cluster_columns    = F, # 列不聚类
                   show_row_names     = T, # 显示行名
                   row_names_side     = "left", # 行名显示在左侧
                   row_names_gp       = gpar(fontsize = 10), # 行名字号为10
                   show_column_names  = F, # 不显示列名（可后期在颜色内AI使得亚型一目了然）
                   column_names_side  = "top", # 列名显示在顶部
                   row_split          = annRow$Category, # 行按照Category进行分割（因子顺序）
                   top_annotation     = top_anno, # 热图顶部注释
                   left_annotation    = left_anno, # 热图左侧注释
                   name               = "mRNA\nExpression", # 热图颜色图例的名称
                   width              = ncol(expMat) * unit(3.5, "mm"), # 热图单元格宽度（稍大于高度，因为所有注释都放在底部，平衡图形纵横比）
                   height             = nrow(expMat) * unit(2.4, "mm")) # 热图单元格高度

pdf(file = "complexheatmap of immunomodulator.pdf", width = 8,height = 12)
invisible(dev.off())
####GSE1446
expr=expr_GSE14468[consam,]
expr1 <- apply(expr, 2, as.numeric)
expr1=as.data.frame(expr1)
rownames(expr1)=rownames(expr)
expr=expr1
expr=na.omit(expr)
expr=as.data.frame(t(expr))
subtype=GSE14468.ntp.pred$clust.res
aas=merge(expr,subtype,by.x=0,by.y=1)
rownames(aas)=aas$Row.names
aas=aas[,-1]
aas$clust=paste0("clust",aas$clust)

expMat <- as.data.frame(t(apply(aas[,setdiff(colnames(aas), "clust")], 2, 
                                function(x) 
                                  tapply(x, 
                                         INDEX = factor(aas$clust), 
                                         FUN = median, 
                                         na.rm = TRUE)))) # 对同一亚型内的样本取中位数

library(ComplexHeatmap)
library(circlize)
library(ChAMPdata)
library(data.table)
??ChAMPdata
###
# 设置热图颜色
heatmap.BlWtRd <- c("#61AACF","#98CADD","#EAEFF6","#F9EFEF","#E9C6C6","#DA9599")

subt=c("clust1","clust2","clust3")
annCol <- data.frame(subtype = subt,
                     row.names = subt)
annCol <- annCol[order(annCol$subtype),,drop = F] # 按照亚型排序
annColors <- list()
annColors[["subtype"]] <- c("clust1" = "#2EC4B6",
                            "clust2" = "#E71D36",
                            "clust3" = "#FF9F1C")
top_anno <- HeatmapAnnotation(df                   = annCol,
                              col                  = annColors,
                              gp                   = gpar(col = "grey80"), # 每个单元格边框为灰色
                              simple_anno_size     = unit(3.5, "mm"), # 注释高3.5毫米
                              show_legend          = F, # 不显示亚型的图例，因为一目了然
                              show_annotation_name = F, # 不显示该注释的名称
                              border               = FALSE) # 不显示注释的外边框
rownames(immunomodulator)=immunomodulator$Gene
immunomodulator=immunomodulator[consam,]
# 创建行注释
annRow <- immunomodulator
annRow[which(annRow$Category == "Co-stimulator"),"Category"] <- "Co-stm" # 这里字符少一些，不会挤在一起，可以后期AI
annRow[which(annRow$Category == "Co-inhibitor"),"Category"] <- "Co-ihb"
annRow[which(annRow$Category == "Cell adhesion"),"Category"] <- "Cell\nadhesion" # 这里换行，不会挤在一起，可以后期AI
annRow[which(annRow$Category == "Antigen presentation"),"Category"] <- "Antigen\npresentation"
annRow$Category <- factor(annRow$Category, levels = c("Co-stm","Co-ihb","Ligand","Receptor","Cell\nadhesion","Antigen\npresentation","Other")) # 由于行需要按照类分割，所以需要定义因子顺序，否则按照字母表
annRow$ICI <- factor(annRow$ICI, levels = c("Inhibitory","N/A","Stimulatory"))
annRowColors <- list("ICI" = c("Inhibitory" = "#3D5C6F","N/A" = "#888888","Stimulatory" = "#E47159"))
left_anno <- HeatmapAnnotation(df                   = annRow[,"ICI",drop = F],
                               which                = "row", # 这里是行注释（默认为列）
                               gp                   = gpar(col = "grey80"), # 每个单元格边框为灰色
                               col                  = annRowColors,
                               simple_anno_size     = unit(3.5, "mm"), # 注释宽3.5毫米
                               show_annotation_name = F,
                               border               = F)


## 绘制表达谱热图（参数下同）

row_normalize <- function(df) {
  normalized_df <- t(apply(df, 1, function(x) {
    min_val <- min(x)
    max_val <- max(x)
    normalized_row <- (x - min_val) / (max_val - min_val) * 2 - 1
    return(normalized_row)
  }))
  colnames(normalized_df) <- colnames(df)
  return(normalized_df)
}

# 调用按行标准化函数
normalized_data <- row_normalize(expMat)

# 打印标准化后的数据
print(normalized_data)




col_expr <- colorRamp2(seq(-1, 1, length = 6), heatmap.BlWtRd) # 创建热图颜色（将热图输入矩阵的最大最小值取5个点，分配颜色红蓝色板；注意矩阵中可能存在的NA值）
hm.expr <- Heatmap(matrix             = as.matrix(normalized_data),
                   col                = col_expr,
                   border             = NA, # 无热图外边框
                   rect_gp = gpar(col = "grey80"), # 热图单元格边框为灰色
                   cluster_rows       = F, # 行不聚类
                   cluster_columns    = F, # 列不聚类
                   show_row_names     = T, # 显示行名
                   row_names_side     = "left", # 行名显示在左侧
                   row_names_gp       = gpar(fontsize = 10), # 行名字号为10
                   show_column_names  = F, # 不显示列名（可后期在颜色内AI使得亚型一目了然）
                   column_names_side  = "top", # 列名显示在顶部
                   row_split          = annRow$Category, # 行按照Category进行分割（因子顺序）
                   top_annotation     = top_anno, # 热图顶部注释
                   left_annotation    = left_anno, # 热图左侧注释
                   name               = "mRNA\nExpression", # 热图颜色图例的名称
                   width              = ncol(expMat) * unit(3.5, "mm"), # 热图单元格宽度（稍大于高度，因为所有注释都放在底部，平衡图形纵横比）
                   height             = nrow(expMat) * unit(2.4, "mm")) # 热图单元格高度

pdf(file = "complexheatmap of immunomodulator.pdf", width = 8,height = 12)
invisible(dev.off())
