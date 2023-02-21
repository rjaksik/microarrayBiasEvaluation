#' Compare probes from two distinct groups
#'
#' @param ProbeData probe data obtained in the
#' @param SampleIDs ids of all samples in the dataset
#' @param label1 label of the first compared feature
#' @param label2 label of the second compared feature
#' @param label1_desc description of the first compared feature
#' @param label2_desc description of the second compared feature
#' @param FilePrefix prefix of the output files
#' @param FeatureValLabel description of the shown feature values used to label the plot axis
#' @param SampleID ID of the sample which will be used to plot detailes sample-level data
#' @param ExpCut expression level cutoff (only probes with signals above this value will be used)
#' @param Version version of the result files
#' @export
probediff = function(ProbeData, SampleIDs, label1, label2, label1_desc, label2_desc, FilePrefix, FeatureValLabel, SampleID, ExpCut=500, Version=1) {

  library('reshape2')

  ProbeData_nodup = ddply(ProbeData,c('ProbeID','Feature','FeatureVal'),function(xx) colMeans(xx[,SampleIDs]))

  G1 = ProbeData_nodup[ProbeData_nodup$Feature==label1,]
  G2 = ProbeData_nodup[ProbeData_nodup$Feature==label2,]
  rownames(G1) = G1$ProbeID
  rownames(G2) = G2$ProbeID
  G1$Feature = G2$Feature = G1$ProbeID = G2$ProbeID = G1$ProbeName  = G2$ProbeName  = G1$ProbeType = G2$ProbeType = NULL
  G1$SystematicName = G2$SystematicName = G1$FeatureVal = G2$FeatureVal = NULL
  ProbeData_diff = log2(G1[rownames(G2),]) - log2(G2)
  ProbeData_mean = (log2(G1[rownames(G2),]) + log2(G2))/2
  comp_label = paste(label1_desc, label2_desc, sep=' vs ')

  #high expression
  sel = rowMeans(G1)>ExpCut
  idx = names(sel[sel])
  ProbeData_diff_high = ProbeData_diff[idx,]
  ProbeData_mean_high = ProbeData_mean[idx,]

  #boxplots of LFC for all samples
  mdata = melt(ProbeData_diff_high)
  p = ggplot(mdata,aes(variable,value)) + geom_boxplot() +theme_bw() +coord_flip() + xlab('Sample') + ylab(paste0('logFC - ',comp_label))
  ggsave(p,file=paste0(FilePrefix,'_CombinedLFC_ExpHigher',ExpCut,'_v',Version,'.png'), scale=2, width=4, height=3, dpi=600)

  mdata = melt(ProbeData_diff)
  p = ggplot(mdata,aes(variable,value)) + geom_boxplot() +theme_bw() +coord_flip() + xlab('Sample') + ylab(paste0('logFC - ',comp_label))
  ggsave(p,file=paste0(FilePrefix,'_CombinedLFC_v',Version,'.png'), scale=2, width=4, height=3, dpi=600)


  #boxplots of norm signal of all samples
  mdata1 = melt(G1)
  mdata1$label = label1_desc
  mdata2 = melt(G2)
  mdata2$label = label2_desc
  mdata = rbind(mdata1,mdata2)
  p = ggplot(mdata,aes(variable,log2(value),color=label)) + geom_boxplot() +theme_bw() +coord_flip() + xlab('Sample') + ylab('log2 Expression')  +
    scale_color_manual(values=polslcolors)
  ggsave(p,file=paste0(FilePrefix,'_CombinedExprs_v',Version,'.png'), scale=2, width=4, height=3, dpi=600)

  #scatter plot - single sample
  mdata2 = data.frame(log2(G1[,SampleID]),log2(G2[,SampleID]))
  mdata2$color='red'
  mdata2$color[mdata2$G1>mdata2$G2]='blue'
  colnames(mdata2) = c('G1','G2')
  p = ggplot(mdata2,aes(G1,G2)) + geom_point(alpha=0.2) +theme_bw() +coord_flip() + xlab(label1_desc) + ylab(label2_desc) +
    geom_abline(intercept = 0, slope=1,linetype=2,color='red') + scale_color_manual(values=polslcolors)
  ggsave(p,file=paste0(FilePrefix,'_Scatter_',SampleID,'_v',Version,'.png'), scale=2, width=3, height=3, dpi=600)


  #MA plot - single sample
  library(affy)
  png(paste0(FilePrefix,'_MAplot_',SampleID,'_v',Version,'.png'),width = 7, height = 7, units = "in", res=600)
  ma.plot(ProbeData_mean[,SampleID], ProbeData_diff[,SampleID], cex=1 )
  dev.off()


  #LOESS regression - all samples
  loessDF = data.frame()
  for (samp in colnames(ProbeData_diff)) {
    df = data.frame(diff=ProbeData_diff[,samp],mean=ProbeData_mean[,samp])
    lo2 <- loess(diff ~ mean, df, control = loess.control(surface = "direct"))
    loessDF = rbind(loessDF,data.frame(x=lo2$x[,1],fitted=lo2$fitted,samp=samp))
  }
  DSplit = strsplit2(loessDF$samp,'_')
  loessDF$cells = DSplit[,1]
  loessDF$dye = DSplit[,3]
  p1=ggplot(loessDF,aes(x,fitted,group=samp,color=cells,linetype=dye))+geom_line()+theme_bw()+scale_color_manual(values=polslcolors) +
    xlab(paste0('log2 Expression - ',label1_desc)) +ylab('logFC')
  p2=ggplot(loessDF,aes(x,group=samp,color=cells,linetype=dye))+geom_density()+theme_bw()+scale_color_manual(values=polslcolors) +
    xlab(paste0('log2 Expression - ',label1_desc)) + ylab('Fraction of probes')
  library(ggpubr)
  figure <- ggarrange(p1, p2,labels = c("A", "B"),ncol = 1)
  ggsave(figure,file=paste0(FilePrefix,'_MAplot-LOESS_v',Version,'.PDF'), scale=2, width=3, height=3) #, dpi=600


  #association between LFC and FeatureVal
  featureTab = unique(ProbeData_nodup[,c('ProbeID','FeatureVal')])
  ProbeData_diff_label = data.frame(ProbeID = rownames(ProbeData_diff),ProbeData_diff)
  ProbeData_diff_label_annot = merge(ProbeData_diff_label,featureTab)
  loessDF = data.frame()
  for (samp in colnames(ProbeData_diff)) {
    df = data.frame(diff=ProbeData_diff_label_annot[,samp],FeatureVal=ProbeData_diff_label_annot$FeatureVal)
    lo2 <- loess(diff~FeatureVal, df, control = loess.control(surface = "direct"))
    loessDF = rbind(loessDF,data.frame(x=lo2$x[,1],fitted=lo2$fitted,samp=samp))
  }
  DSplit = strsplit2(loessDF$samp,'_')
  loessDF$cells = DSplit[,1]
  loessDF$dye = DSplit[,3]
  p1=ggplot(loessDF,aes(x,fitted,group=samp,color=cells,linetype=dye))+geom_line()+theme_bw()+scale_color_manual(values=polslcolors) +
    xlab(FeatureValLabel) +ylab('logFC')
  p2=ggplot(loessDF,aes(x,group=samp,color=cells,linetype=dye))+geom_density()+theme_bw()+scale_color_manual(values=polslcolors) +
    xlab(FeatureValLabel) + ylab('Fraction of probes')
  figure <- ggarrange(p1, p2,labels = c("A", "B"),ncol = 1)
  ggsave(figure,file=paste0(FilePrefix,'_FeatureVal-LFC-LOESS_v',Version,'.pdf'), scale=2, width=3, height=3)

}
