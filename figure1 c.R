##figure1 c
# mutational frequency comparison
mut.aml <- compMut(moic.res     = cmoic.aml
                    mut.matrix   = aml.tcga$mut.status, 
                    doWord       = TRUE, 
                    doPlot       = TRUE, 
                    freq.cutoff  = 0.01, 
                    p.adj.cutoff = 1, 
                    innerclust   = TRUE, 
                    annCol       = NULL,
                    annColors    = NULL, 
                    width        = 6, 
                    height       = 6,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")
