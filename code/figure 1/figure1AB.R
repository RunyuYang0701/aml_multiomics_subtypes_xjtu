###AML Multi-Omics Subtyping Project
##Author: yry
##Date: March 28, 2024
##Loading Multi-Omics Analysis R Packages
library(MOVICS)
aml.tcga=integrated_data
names(aml.tcga)
surv.info <- aml.tcga$clin.info
tmp       <- aml.tcga$fpkm_mrna# get expression data 
surv.info$futime=surv.info$OS.time
surv.info$fustat=surv.info$OS
elite.tmp <- getElites(dat       = tmp,
                       method    = "cox",
                       surv.info = surv.info, 
                       p.cutoff  = 0.001,
                       elite.num = 100) 
dim(elite.tmp$elite.dat) 
mRNA=elite.tmp$elite.dat
mRNA_cox=elite.tmp$unicox.res
tmp       <- aml.tcga$fpkm_lncrna# get expression data 
surv.info$futime=surv.info$OS.time
surv.info$fustat=surv.info$OS
elite.tmp <- getElites(dat       = tmp,
                       method    = "cox",
                       surv.info = surv.info, 
                       p.cutoff  = 0.001,
                       elite.num = 100) 
dim(elite.tmp$elite.dat) 
lncrna=elite.tmp$elite.dat
lncrna_cox=elite.tmp$unicox.res
##筛选cnv
tmp       <- aml.tcga$cna# get expression data 
surv.info$futime=surv.info$OS.time
surv.info$fustat=surv.info$OS
elite.tmp <- getElites(dat       = tmp,
                       method    = "cox",
                       surv.info = surv.info, 
                       p.cutoff  = 0.001,
                       elite.num = 100) 
dim(elite.tmp$elite.dat) 
cna=elite.tmp$elite.dat
cna_cox=elite.tmp$unicox.res
tmp       <- aml.tcga$meth.beta# get expression data 
surv.info$futime=surv.info$OS.time
surv.info$fustat=surv.info$OS
elite.tmp <- getElites(dat       = tmp,
                       method    = "cox",
                       surv.info = surv.info, # 生存信息，列名必须有'futime'和'fustat'
                       p.cutoff  = 0.001,
                       elite.num = 100) # 此时这个参数也是不起作用的

dim(elite.tmp$elite.dat) 
meth=elite.tmp$elite.dat
meth_cox=elite.tmp$unicox.res

tmp       <- aml.tcga$mut.status
rowSums(tmp)
elite.tmp <- getElites(dat       = tmp,
                       method    = "freq", 
                       elite.num = 5, 
                       elite.pct = 0.1) 
rowSums(elite.tmp$elite.dat)
mut=elite.tmp$elite.dat
library(dplyr)
mRNA=na.omit(mRNA)
mRNA <- mutate_if(mRNA, is.numeric, function(x) log2(x + 1))
mRNA=df_log
mo.data=list(mRNA,lncrna,meth,cna,mut)


optk.aml <- getClustNum(data        = mo.data, 
                        is.binary   = c(F,F,F,F,T), 
                        try.N.clust = 2:8, 
                        fig.name    = "CLUSTER NUMBER OF TCGA-LAML")
# perform multi-omics integrative clustering with the rest of 9 algorithms
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster","iClusterBayes"), # 9种算法
                         N.clust     = 3,
                         type        = c("gaussian", "gaussian", "gaussian", "gaussian","binomial"))
save(moic.res.list, file = "moic.res.list.rda")
save(moic.res.list, file = "cmoic.laml.rda")
cmoic.brca <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")
getSilhouette(sil      = cmoic.brca$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)
indata <- mo.data
indata$meth <- log2(indata$meth / (1 - indata$meth))
new_names <- c("mRNA", "lncRNA", "cna","meth","mut")
names(indata) <- new_names
# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,T,F)) # no scale for mutation
iClusterBayes.res=cmoic.brca$clust.res
feat   <- iClusterBayes.res
rownames(feat)=feat$feature
feat1  <- feat[which(feat$dataset == "dat1"),][1:20,"feature"] 
feat2  <- feat[which(feat$dataset == "dat2"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "dat4"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "dat3"),][1:5,"feature"]
feat5  <- feat[which(feat$dataset == "dat5"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4,feat5)

#set the colors for mult-omics
mRNA.col   <- c("#31B29E", "black","#EB6C5A")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
cna.col  <-  c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
meth.col   <- c("#0099CC", "white","#CC0033")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, cna.col,meth.col, mut.col)
rownames(anncol)=anncol$...1
annCol=anncol
write.csv(annCol,file="anncol.csv")
annCol    <- surv.info[,c(51,74,1), drop = FALSE]
annCol=annCol[,-c(1:2)]
# generate corresponding colors for sample annotation
table(annCol$cytogenetics)
annColors <- list(Age    = circlize::colorRamp2(breaks = c(min(annCol$Age),
                                                           max(annCol$Age)), 
                                                colors = c("#7FC97F",  "#E00115")),
                  
                  cytogenetics = c("Favorable"    = "#C3AEA5",
                                   "Intermediate/Normal"    = "#A6877A",
                                   "Poor"    = "#7B4B38",
                                   "unknow"    = "#C2BEBC"))

getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","CNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.expr","lncRNA.expr","CNA","M value","Mutated"),
             clust.res     = iClusterBayes.res$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")


clust=cmoic.brca$clust.res

clinical=surv.info[,c(3,4)]

aas=merge(clust,clinical,by.x=0,by.y=0)
aas$OS.time=aas$OS.time/30
fit <- survfit(Surv(aas$OS.time,aas$OS)~aas$clust,data=aas)
p2 <-ggsurvplot(fit, pval = TRUE,linetype =1,xlab="Months", ylab="Overall survival proability",
                palette =c("#2EC4B6","#E71D36","#FF9F1C"),
                font.x = c(14, "italic", "darkblue"),
                font.y = c(14, "italic", "darkblue"),
                font.tickslab= c(13, "italic", "black"),
                font.legend= c(10, "italic", "black") ,
                legend.title="",legend=c(0.7,0.7),censor.shape=124,censor.size=2,
                conf.int =F) 
p2
aas$OS.time=aas$OS.time/30
