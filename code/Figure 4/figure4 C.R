data=tcell_state[,c(3,11)]
colnames(data) <- c("Group","Values")
data$Values <- as.numeric(data$Values)
data$Group=paste0("cluster",data$Group)
data$Group <- factor(data$Group,levels = c("cluster1","cluster2","cluster3"))
x_value <- c(0.9, 1.1)
for (i in 1:3) {
  x_value <- append(x_value, x_value[((i-1)*2+1):(i*2)]+1)
}

med_values <- rep(c(median(data[which(data$Group == "cluster1"),]$Values),
                    median(data[which(data$Group == "cluster2"),]$Values),
                    median(data[which(data$Group == "cluster3"),]$Values)),
                  rep(2,3))

tmp_data <- as.data.frame(cbind(x_value,med_values,
                                rep(c("cluster1","cluster2","cluster3"),
                                    rep(2,3))))
tmp_data$x_value <- as.numeric(tmp_data$x_value)
tmp_data$med_values <- as.numeric(tmp_data$med_values)
tmp_data$V3 <- factor(tmp_data$V3, levels = c("cluster1","cluster2","cluster3"))

# 绘图:
p1 <- ggplot(data,aes(Group,Values))+
  # 提琴图：
  geom_violin(aes(fill=Group),color="white")+
  # 内部箱线图：
  geom_boxplot(fill="#a6a7ac",color="#a6a7ac",
               width = 0.1,outlier.shape = NA)+
  # 填充颜色模式：
  scale_color_manual(values = c("#66C1B3", "#E47B85", "#F3C183", "#417bb9")) +
  scale_fill_manual(values = c("#66C1B3", "#E47B85", "#F3C183", "#417bb9")) +
  # 设置主题：
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +# 坐标轴label：
xlab("Status")+
  ylab("Gene expression")+
  # y轴刻度标签：
  scale_y_continuous(breaks = seq(0.5,1,0.1))+
  # 添加注释：
  annotate("text", label = "bolditalic(MTAP)", parse = TRUE, 
           x = 3, y = 13, size = 4, colour = "black")+
  annotate("text", label = "5346",
           x = 1.3, y = 7, size = 4, colour = "#a6a7ac")+
  annotate("text", label = "877",  
           x = 2.3, y = 8, size = 4, colour = "#a6a7ac")+
  annotate("text", label = "2619",
           x = 3.3, y = 7, size = 4, colour = "#a6a7ac")+
  annotate("text", label = "774",
           x = 4.3, y = 2.8, size = 4, colour = "#a6a7ac")+
  # 差异显著性检验：
  geom_signif(comparisons = list(c("LOH", "HD"), # 哪些组进行比较
                                 c("WT", "HD")),
              map_signif_level=T, # T显示显著性，F显示p值
              tip_length=c(0,0), # 修改显著性线两端的长短
              y_position = c(0,-0.5), # 设置显著性线的位置高度
              size = 0.5, # 修改线的粗细
              textsize = 4, # 修改显著性标记的大小
              extend_line = -0.1, # 线的长短：默认为0；
              color = "#a6a7ac", # 线的颜色
              test = "t.test") # 检验的类型
