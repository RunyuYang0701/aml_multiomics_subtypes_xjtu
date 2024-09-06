#####Figure 1D

clust=cmoic.brca$clust.res
clinical=integrated_data$clin.info
clinical=clinical[,c(1,74)]
aas=merge(clust,clinical,by.x=0,by.y=0)
aas=na.omit(aas)
table(aas$clust,aas$acute_myeloid_leukemia_calgb_cytogenetics_risk_category)
write.csv(aas,"tcga_eln.csv")
library(ggplot2)
library(gtools)
library(ggalluvial)
df=ELN
df$Celltype = factor(df$Celltype, levels = unique(df$Celltype))
df$Sample = factor(df$Sample, levels = mixedsort(unique(df$Sample)))

p <- ggplot(as.data.frame(df), aes(x = Sample, y = Proportion, fill = Celltype, stratum = Celltype, alluvium = Celltype)) +
  geom_flow(width = 0.5, alpha = 0.3, knot.pos=0, color = 'white') +
  geom_col(width = 0.5, color = 'white') +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values = c("#2EC4B6","#E71D36","#FF9F1C")) +
  xlab("") + ylab("Cell proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 20, face='bold'),
        axis.text.y = element_text(size = 20, face='bold'),
        axis.title.y = element_text(size = 24, face='bold'),
        legend.text=element_text(size=20),
        legend.title=element_text(size=24)
  )
ggsave("tcga_mutation.pdf", p, width = 8, height = 6)


####原始数据的投射---fab分型
clust=cmoic.brca$clust.res
clinical=integrated_data$clin.info
clinical=clinical[,c(1,51)]
aas=merge(clust,clinical,by.x=0,by.y=0)
aas=na.omit(aas)
table(aas$clust,aas$leukemia_french_american_british_morphology_code)
write.csv(aas,"tcga_eln.csv")

mut=integrated_data$mut.status
mut_flt=mut["FLT3",]
mult_npm1=mut["NPM1",]
mut_cebpa=mut["CEBPA",]
mut_TP53=mut["TP53",]
mut=rbind(mut_flt,mut_cebpa)
mut=rbind(mut,mut_TP53)
mut=rbind(mut,mult_npm1)
mut=as.data.frame(t(mut))
mut=merge(clust,mut,by.x=0,by.y=0)
table(mut$NPM1,mut$clust)
write.csv(aas,"tcga_eln.csv")
library(ggplot2)
library(gtools)
library(ggalluvial)
df=ELN
df$Celltype = factor(df$Celltype, levels = unique(df$Celltype))
df$Sample = factor(df$Sample, levels = mixedsort(unique(df$Sample)))

p <- ggplot(as.data.frame(df), aes(x = Sample, y = Proportion, fill = Celltype, stratum = Celltype, alluvium = Celltype)) +
  geom_flow(width = 0.5, alpha = 0.3, knot.pos=0, color = 'white') +
  geom_col(width = 0.5, color = 'white') +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values = c("#EE8172","#83D0E2","#4DBDAB","#7788AC")) +
  xlab("") + ylab("Cell proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 20, face='bold'),
        axis.text.y = element_text(size = 20, face='bold'),
        axis.title.y = element_text(size = 24, face='bold'),
        legend.text=element_text(size=20),
        legend.title=element_text(size=24)
  )
ggsave("tcga_ELN.pdf", p, width = 8, height = 6)
