####beat_aml
expr_beat=beataml_waves1to4_norm_exp_dbgap
expr_beat=expr_beat[,-c(1,3,4)]
rownames(expr_beat)=expr_beat$display_label
expr_beat=expr_beat[,-1]
rownames(expr_beat)=beataml_waves1to4_norm_exp_dbgap$display_label
expr_beat=as.data.frame(expr_beat)
save(expr_beat,file = "expr_beat.rda")
beataml.ntp.pred <- runNTP(expr       = expr_beat,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR YAU") 
save(beataml.ntp.pred,file = "beat_ntp_pred.rda")
####预后
clinical_beat=beataml_wv1to4_clinical[,c(3,62,63)]
clust=yau.ntp.pred$clust.res
aas=merge(clust,clinical_beat,by.x=1,by.y=1)
table(aas$vitalStatus)
aas=aas[which(aas$vitalStatus!="Unknown"),]
aas$OS=ifelse(aas$vitalStatus=="Dead",1,0)
aas$OS.time=as.numeric(aas$overallSurvival)
save(aas,file = "clust_os_beat.rda")
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
aas$OS.time=aas$OS.time/30
####根据投射明确临床指标
clinical_beat=beataml_wv1to4_clinical[,c(3,35,55,62,89,91,94,78)]
aas=merge(clust,clinical_beat,by.x=1,by.y=1)
###eln分级
aas1=aas[,c(1,2,9)]
aas1=na.omit(aas1)
table(aas1$responseToInductionTx)
aas_adver=aas1[which(aas1$ELN2017=="Adverse"),]
aas_Favor=aas1[which(aas1$ELN2017=="Favorable"),]
aas_inter=aas1[which(aas1$ELN2017=="Intermediate"),]
aas_elc=rbind(aas_adver,aas_Favor)
aas_elc=rbind(aas_elc,aas_inter)
table(aas_elc$clust,aas_elc$ELN2017)
###对诱导治疗的效果
aas_cr=aas1[which(aas1$responseToInductionTx=="Complete Response"),]
aas_rr=aas1[which(aas1$responseToInductionTx=="Refractory"),]
aas_induction=rbind(aas_cr,aas_rr)
table(aas_induction$clust,aas_induction$responseToInductionTx)
###fab分型
table(aas1$fabBlastMorphology)
aas_fab=aas1[which(aas1$fabBlastMorphology!="M5b"),]
table(aas_fab$clust,aas_fab$fabBlastMorphology)
###flit-3
aas1=aas[,c(1,2,6)]
aas1=na.omit(aas1)
table(aas1$`FLT3-ITD`)
table(aas1$clust,aas1$`FLT3-ITD`)
###npm1
aas1=aas[,c(1,2,7)]
aas1=na.omit(aas1)
table(aas1$NPM1)
table(aas1$clust,aas1$NPM1)
##tp53
aas1=aas[,c(1,2,8)]
table(aas1$TP53)
aas_p53_p=na.omit(aas1)
aas_p53_p$status="positive"
aas1[is.na(aas1)] <- "negative/unkonwn"
aas1_p53_n=aas1[which(aas1$TP53=="negative/unkonwn"),]
aas1_p53_n$status="negative/unkonwn"
aas_p53=rbind(aas1_p53_n,aas_p53_p)
table(aas_p53$clust,aas_p53$status)

m