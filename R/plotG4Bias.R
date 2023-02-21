#' Plot G4 bias statistics
#'
#' @param expData - result of the procCustArray function
#' @param probeData - data.frame with detailed probe information
#' @param plotNamePrefix - prefix used to name all output plots
#' @export
plotG4Bias = function(expData,probeData,plotNamePrefix=NULL) {

  exptabAnnot_G4 = exptabAnnotSub[exptabAnnotSub$ProbeType=='G4',]
  exptabAnnot_G4_annot = merge(exptabAnnot_G4,probeData,by.x='ProbeName',by.y='Probe.nID',all.x=T,all.y=F)
  exptabAnnot_G4_annot$FeatureVal = exptabAnnot_G4_annot$G4.Score

  probediff(exptabAnnot_G4_annot, SampleIDs, label1='Y', label2='N',
            label1_desc = 'With G4 motif', label2_desc = 'No G4 motif', FeatureValLabel = 'G4 Score',
            SampleID = expData$sampleIDs[1], FilePrefix = plotNamePrefix, ExpCut=500, Version=1)
}
