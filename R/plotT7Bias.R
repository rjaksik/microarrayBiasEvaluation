#' Plot T7 bias statistics
#'
#' @param expData - result of the procCustArray function
#' @param probeData - data.frame with detailed probe information
#' @param plotNamePrefix - prefix used to name all output plots
#' @export
plotT7Bias = function(expData,probeData,plotNamePrefix=NULL) {

  exptabAnnot_T7 = expData$exptab %>% filter(ProbeType=='T7')
  exptabAnnot_T7_annot = merge(exptabAnnot_T7,probeData,by.x='ProbeName',by.y='Probe.nID',all.x=T,all.y=F)
  exptabAnnot_T7_annot$FeatureVal = exptabAnnot_T7_annot$distance

  probediff(exptabAnnot_T7_annot, SampleIDs, label1='Y', label2='N',
            label1_desc = 'With T7 motif', label2_desc = 'No T7 motif', FeatureValLabel = 'Distance from T7',
            SampleID=expData$sampleIDs[1], FilePrefix = plotNamePrefix, ExpCut=500, Version=1)
}
