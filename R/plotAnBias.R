#' Plot An bias statistics
#'
#' @param expData - result of the procCustArray function
#' @param probeData - data.frame with detailed probe information
#' @param plotNamePrefix - prefix used to name all output plots
#' @export
plotAnBias = function(expData,probeData,plotNamePrefix=NULL) {

  exptabAnnot_An = expData$exptab %>% filter(ProbeType=='An')

  exptabAnnot_An_annot = merge(exptabAnnot_An,probeData,by.x='ProbeName',by.y='Probe.nID',all.x=T,all.y=F)
  exptabAnnot_An_annot$FeatureVal = exptabAnnot_An_annot$An.Length
  exptabAnnot_An_annot = exptabAnnot_An_annot[exptabAnnot_An_annot$An.Length<30,]

  probediff(exptabAnnot_An_annot, expData$sampleIDs, label1='U', label2='D',
            label1_desc = 'Upstream from (A)n', label2_desc = 'Downstream from (A)n', FeatureValLabel = '(A)n length',
            SampleID=expData$sampleIDs[1], FilePrefix = plotNamePrefix, ExpCut=500, Version=2)
}
