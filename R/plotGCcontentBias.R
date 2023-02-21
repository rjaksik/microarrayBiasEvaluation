#' Plot GC-content bias statistics
#'
#' @param expData - result of the procCustArray function
#' @param plotName - name of the output plot
#' @export
plotGCcontentBias = function(expData,plotName=NULL) {
  exptabAnnot_GC = expData$exptab %>% filter(ProbeType=='GC')

  exptab_mean = ddply(exptabAnnot_GC,'Feature',function(x) {colMeans(log2(x[,expData$sampleIDs]))})
  exptabAnnot_GC_meanMD = reshape2::melt(exptabAnnot_GC_mean)
  exptabAnnot_GC_meanMD$cells = limma::strsplit2(exptabAnnot_GC_meanMD$variable,'_')[,1]
  exptabAnnot_GC_meanMD$Feature = as.numeric(exptabAnnot_GC_meanMD$Feature )
  exptabAnnot_GC_meanMD = exptabAnnot_GC_meanMD[exptabAnnot_GC_meanMD$Feature>=10 & exptabAnnot_GC_meanMD$Feature<=50,]
  p = ggplot(exptabAnnot_GC_meanMD) + geom_line(aes(Feature/60,value,color=cells,group=variable)) +
    xlab('GC content [%]') + ylab('log2 expression')  +scale_color_manual(values=polslcolors) +theme_bw()
  if (!is.null(plotName)) {
    ggsave(p,file=plotName, scale=2, width=3, height=2)
  }
  return(p)
}
