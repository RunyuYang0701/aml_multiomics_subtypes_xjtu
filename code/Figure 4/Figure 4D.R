
library(ggplot2)
library(tidyverse)


expr=expr_beat[,aas$samID]
mtr=expr
data=aas
SE <- SummarizedExperiment::SummarizedExperiment(assays = list(fpkm = mtr),
                                                 colData = data)
aa_beat=pred_response(SE=SE,Signature = ipt,
                      method = "Weighted_mean",threshold = 0.8,
                      PT_drop = FALSE,sort_by = "Customed.Signature",
                      group_by = "Customed.Signature",show.Observed = TRUE,
                      rankscore = FALSE)
xx=aa_beat[[1]][["data"]]
data_beat=merge(data,xx,by.x=1,by.y=1)
table(data_beat$clust,data_beat$Value)




# 构建数据：
data <- data.frame(cluster1= c(167,60), cluster2 = c(111,74), cluster3 = c(167,70),
                   group = c('NR',"R")) %>% 
  pivot_longer(cols = !group, names_to = "X", values_to = "count")

# 绘图：
p1 <- ggplot(data)+
  geom_bar(aes(rev(X), count, fill = group), color = "#f3f4f4",
           position = "fill", stat = "identity", size = 1)+
  # 修改填充颜色：
  scale_fill_manual(values = c("#6AA8C5", "#EFDDDD", "#bbbdc0"),
                    # 图例标签：
                    labels=rev(c("TC 0(<1%)","TC 1(<5%)","TC 2+(>=5%)")))+
  # 添加星号注释：
  annotate("text", x = 2, y = 0.85, label="*", size = 5)+
  annotate("text", x = 1.5, y = 1.05, label=expression("*"~italic("P=0.003")), size = 4)+
  # 标题和副标题：
  ggtitle("Immune cell\nPD-L1", subtitle = "(by SP142 IHC)")+
  # 难点：学会使用expression函数：
  scale_x_discrete(labels = c(expression(atop(bold("Q4"), "(n=74)")),
                              expression(atop(bold("Q1"), "(n=74)"))))+
  # x轴和y轴标签
  xlab("")+
  ylab("")+
  # 设置主题：
  theme_bw()+
  theme(panel.grid = element_blank(),
        # 标题和副标题居中
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, face = "italic"),
        # 修改背景色：
        panel.background = element_rect(fill = "#f3f4f4")
  )+
  # 图例调整：
  # 图例顺序：
  guides(fill=guide_legend(reverse=TRUE))+
  # 图例标题：
  labs(fill="IC level")

p1
data <- matrix(c(167,70,111,74,167,60),
               nrow = 3, byrow = TRUE)
test_result <- chisq.test(data)

# 查看结果
print(test_result)





mtr=exprSet
data=clinical_GSE71014
colnames(data)[1]="sample id"
aa=GSE71014.ntp.pred$clust.res
data=merge(data,aa,by.x=1,by.y=1)
SE_GSE71014 <- SummarizedExperiment::SummarizedExperiment(assays = list(fpkm = mtr),
                                                          colData = data)
aa=pred_response(SE=SE_GSE71014,Signature = ipt,
                 method = "Weighted_mean",threshold = 0.8,
                 PT_drop = FALSE,sort_by = "Customed.Signature",
                 group_by = "Customed.Signature",show.Observed = TRUE,
                 rankscore = FALSE)
xx=aa[[1]][["data"]]
data_71014=merge(data,xx,by.x=1,by.y=1)
table(data_71014$clust,data_71014$Value)




# 构建数据：
data <- data.frame(cluster1= c(18,11), cluster2 = c(17,16), cluster3 = c(34,8),
                   group = c('NR',"R")) %>% 
  pivot_longer(cols = !group, names_to = "X", values_to = "count")

# 绘图：
p1 <- ggplot(data)+
  geom_bar(aes(rev(X), count, fill = group), color = "#f3f4f4",
           position = "fill", stat = "identity", size = 1)+
  # 修改填充颜色：
  scale_fill_manual(values = c("#6AA8C5", "#EFDDDD", "#bbbdc0"),
                    # 图例标签：
                    labels=rev(c("TC 0(<1%)","TC 1(<5%)","TC 2+(>=5%)")))+
  # 添加星号注释：
  annotate("text", x = 2, y = 0.85, label="*", size = 5)+
  annotate("text", x = 1.5, y = 1.05, label=expression("*"~italic("P=0.003")), size = 4)+
  # 标题和副标题：
  ggtitle("Immune cell\nPD-L1", subtitle = "(by SP142 IHC)")+
  # 难点：学会使用expression函数：
  scale_x_discrete(labels = c(expression(atop(bold("Q4"), "(n=74)")),
                              expression(atop(bold("Q1"), "(n=74)"))))+
  # x轴和y轴标签
  xlab("")+
  ylab("")+
  # 设置主题：
  theme_bw()+
  theme(panel.grid = element_blank(),
        # 标题和副标题居中
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, face = "italic"),
        # 修改背景色：
        panel.background = element_rect(fill = "#f3f4f4")
  )+
  # 图例调整：
  # 图例顺序：
  guides(fill=guide_legend(reverse=TRUE))+
  # 图例标题：
  labs(fill="IC level")

p1
data <- matrix(c(34,8,17,16,18,11),
               nrow = 3, byrow = TRUE)
test_result <- chisq.test(data)

# 查看结果
print(test_result)

