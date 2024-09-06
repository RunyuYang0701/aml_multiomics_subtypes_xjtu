###亚型预测
library(MOVICS)
GSE132511.ntp.pred <- runNTP(expr       = expr_gse132511,
                             templates  = marker.up$templates, # the template has been already prepared in runMarker()
                             scaleFlag  = TRUE, # scale input data (by default)
                             centerFlag = TRUE, # center input data (by default)
                             doPlot     = TRUE, # to generate heatmap
                             fig.name   = "NTP HEATMAP FOR 132511") 