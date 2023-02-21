#' Plot length bias statistics
#'
#' @param expData - result of the procCustArray function
#' @param plotName - name of the output plot
#' @export
plotLengthBias = function(expData,plotName=NULL) {

  #prepare data
  exptabAnnot_Len = expData$exptab %>% tibble() %>%
    filter(ProbeType=='LN') %>%
    mutate(Feature = as.numeric(Feature))

  #scale for each gene
  exptabAnnot_Len_norm = ddply(exptabAnnot_Len,'ProbeID',function(x) {
    scaleVal = x[x$Feature==sort(x$Feature)[1],]
    x[,'Feature'] = x[,'Feature']-scaleVal$Feature
    for (SampleID in SampleIDs) {
      x[,SampleID] = log2(x[,SampleID] / scaleVal[[SampleID]])
    }
    return(x)
  })
  #create histogram
  hst = hist(exptabAnnot_Len_norm$Feature,100,plot=F)
  winLen = hst$breaks[2] - hst$breaks[1]
  PlotTable = data.frame()
  for (int in hst$breaks[-1]) {
    texptabAnnot_Len_norm = exptabAnnot_Len_norm[exptabAnnot_Len_norm$Feature<int & exptabAnnot_Len_norm$Feature>(int-winLen),]
    tPlotTable = data.frame(int,N=nrow(texptabAnnot_Len_norm),SampleID=SampleIDs,Exprs=colMeans(texptabAnnot_Len_norm[,SampleIDs]))
    PlotTable = rbind(PlotTable,tPlotTable)
  }
  PlotTable$cells = strsplit2(PlotTable$SampleID,'_')[,1]

  #limit plots to 2500bp
  PlotTable_sub = PlotTable[PlotTable$int<2500,]

  #create both plots
  p1 = ggplot(PlotTable_sub) + geom_line(aes(int,Exprs,color=cells,group=SampleID)) + xlab('Distance from 3\'') + ylab('standardized signal')  +scale_color_manual(values=polslcolors) +theme_bw()+ scale_x_reverse()
  p2 = ggplot(PlotTable_sub) + geom_line(aes(int,N,color=cells))+theme_bw()+ xlab('Distance from the 3\'-end of transcript') + ylab('probe sets')+scale_color_manual(values=polslcolors)+ scale_x_reverse()
  library(ggpubr)
  figure <- ggarrange(p1, p2,labels = c("A", "B"),ncol = 1)

  if (!is.null(plotName)) {
    ggsave(figure,file=, scale=2, width=3, height=3)
  }

  return(figure)
}
