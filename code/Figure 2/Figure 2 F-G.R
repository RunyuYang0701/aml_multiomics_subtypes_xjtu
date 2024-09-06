beataml.ntp.pred <- runNTP(expr       = expr_GSE10358,
                           templates  = marker.up$templates, # the template has been already prepared in runMarker()
                           scaleFlag  = TRUE, # scale input data (by default)
                           centerFlag = TRUE, # center input data (by default)
                           doPlot     = TRUE, # to generate heatmap
                           fig.name   = "NTP HEATMAP FOR YAU") 
save(beataml.ntp.pred,file = "GSE10358_ntp_pred.rda")
library("survival")
library("survminer")
fit <- survfit(Surv(aas$OS.time,aas$OS)~aas$clust,data=aas)
p2 <-ggsurvplot(fit, pval = TRUE,linetype =1,xlab="Months", ylab="Overall survival proability",
                palette =c("#2EC4B6","#E71D36","#FF9F1C"),
                font.x = c(14, "italic", "darkblue"),
                font.y = c(14, "italic", "darkblue"),
                font.tickslab= c(13, "italic", "black"),
                font.legend= c(10, "italic", "black") ,
                legend.title="",legend=c(0.7,0.7),censor.shape=124,censor.size=2,
                conf.int =F) #不显示置信区间
p2