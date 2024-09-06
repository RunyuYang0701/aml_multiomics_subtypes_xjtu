###Figure 1 E
######
####c3subtypes

library(dplyr)
library(network)
library(ggraph)
library(igraph)

data <- read.csv("string_interactions_short.csv")
data=string_interactions_short
edges <- data.frame(from=data[,1],to=data[,2])
colnames(edges)=c("from","to")
rownames(gene_info) <- gene_info$`Gene name`
unique_gene <- c(edges$from,edges$to)
unique_gene <- unique(unique_gene)
##
gene_info <- gene_info[unique_gene, ]
View(gene_info)
gene_info <- na.omit(gene_info)
unique_gene %in% gene_info$`Gene name`

unique_gene[52]
DIPK1C <- data.frame(Gene_name='H3C12',
                     Regulated_Type = "CS3")

gene_vertex <- rbind(gene_vertex , DIPK1C)
colnames(gene_vertex)[2]="Regulated_Type"
g <- graph.data.frame(edges,gene_vertex,directed = FALSE)
plot(g)
rw <- walktrap.community(g, weights = E(g)$weight)
set.seed(123)
plot(rw, g)



plot(rw, g,
     layout=layout.auto,
     vertex.label.cex=0.6,
     margin=-0.03)
rw
rw$membership
rw$names
module_gene <- data.frame(group = rw$membership,
                          gene = rw$names)

module_gene <- module_gene[order(module_gene$group),]
table(module_gene$group)

library(clusterProfiler)
group <- data.frame(gene=module_gene$gene,
                    group=module_gene$group)

Gene_ID <- bitr(module_gene$gene, 
                fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")
library("clusterProfiler")
library("org.Hs.eg.db")
#构建文件并分析
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
#KEGG分析
module_KEGG <- compareCluster(ENTREZID~group, 
                              data=data, 
                              fun="enrichPathway", 
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05)
??compareCluster
dotplot(module_KEGG, showCategory=5,font.size = 8)
module_KEGG_results <- module_KEGG@compareClusterResult
write.csv(module_KEGG_results, file = "module_KEGG_results.csv")



####c1 subtypes
library(dplyr)
library(network)
library(ggraph)
library(igraph)

data <- read.csv("string_interactions_short.csv")
data=string_interactions_short
edges <- data.frame(from=data[,1],to=data[,2])
colnames(edges)=c("from","to")
rownames(gene_info) <- gene_info$`Gene name`
unique_gene <- c(edges$from,edges$to)
unique_gene <- unique(unique_gene)
##
gene_info <- gene_info[unique_gene, ]
View(gene_info)
gene_info <- na.omit(gene_info)
unique_gene %in% gene_info$`Gene name`

unique_gene[16]
DIPK1C <- data.frame(Gene_name='CUX1',
                     Regulated_Type = "CS1")
gene_vertex <- rbind(gene_vertex , DIPK1C)
colnames(gene_vertex)[2]="Regulated_Type"
g <- graph.data.frame(edges,gene_vertex,directed = FALSE)
plot(g)
rw <- walktrap.community(g, weights = E(g)$weight)
set.seed(123)
plot(rw, g)



plot(rw, g,
     layout=layout.auto,
     vertex.label.cex=0.6,
     margin=-0.03)
rw
rw$membership#模块
rw$names#基因
module_gene <- data.frame(group = rw$membership,
                          gene = rw$names)

module_gene <- module_gene[order(module_gene$group),]
table(module_gene$group)

#这里我们使用clusterprofile进行富集分析

library(clusterProfiler)
group <- data.frame(gene=module_gene$gene,
                    group=module_gene$group)

Gene_ID <- bitr(module_gene$gene, 
                fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")
library("clusterProfiler")
library("org.Hs.eg.db")
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
module_KEGG <- compareCluster(ENTREZID~group, 
                              data=data, 
                              fun="enrichPathway", 
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05)
??compareCluster
dotplot(module_KEGG, showCategory=5,font.size = 8)
module_KEGG_results <- module_KEGG@compareClusterResult
write.csv(module_KEGG_results, file = "module_c1_results.csv")




####c2subtypes
library(dplyr)
library(network)
library(ggraph)
library(igraph)

data <- read.csv("string_interactions_short.csv")
data=string_interactions_short
edges <- data.frame(from=data[,1],to=data[,2])
colnames(edges)=c("from","to")
rownames(gene_info) <- gene_info$`Gene name`
unique_gene <- c(edges$from,edges$to)
unique_gene <- unique(unique_gene)
##
gene_info <- gene_info[unique_gene, ]
View(gene_info)
gene_info <- na.omit(gene_info)
unique_gene %in% gene_info$`Gene name`
gene_vertex=gene_info
unique_gene[10]
DIPK1C <- data.frame(Gene_name='C16orf96',
                     Regulated_Type = "CS2")
gene_vertex <- rbind(gene_vertex , DIPK1C)
colnames(gene_vertex)[1]="Gene_name"
g <- graph.data.frame(edges,gene_vertex,directed = FALSE)
plot(g)
rw <- walktrap.community(g, weights = E(g)$weight)
set.seed(123)
plot(rw, g)



plot(rw, g,
     layout=layout.auto,
     vertex.label.cex=0.6,
     margin=-0.03)
rw
rw$membership
rw$names
module_gene <- data.frame(group = rw$membership,
                          gene = rw$names)

module_gene <- module_gene[order(module_gene$group),]
table(module_gene$group)


library(clusterProfiler)
group <- data.frame(gene=module_gene$gene,
                    group=module_gene$group)

Gene_ID <- bitr(module_gene$gene, 
                fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")
library("clusterProfiler")
library("org.Hs.eg.db")

data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')

module_KEGG <- compareCluster(ENTREZID~group, 
                              data=data, 
                              fun="enrichPathway", 
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05)
??compareCluster
dotplot(module_KEGG, showCategory=5,font.size = 8)
module_KEGG_results <- module_KEGG@compareClusterResult
write.csv(module_KEGG_results, file = "module_c2_results.csv")
