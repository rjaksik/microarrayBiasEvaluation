#' Preprocess data from custom Agilent microarray
#'
#' @param targets_file - table specifying association between samples and data files, for details see limma::readTargets function description
#' @param path - directory containing the targets file and all data files
#' @export
procCustArray = function(targets_file,path) {

  #read targets and data files
  targets <- readTargets(targets_file,path=path)
  RG <- read.maimages(targets,
                      columns = list(G = "gMedianSignal", Gb = "gBGMedianSignal", R = "rMedianSignal",
                                     Rb = "rBGMedianSignal"),
                      other.columns = c('gIsWellAboveBG',	'rIsWellAboveBG'),
                      annotation = c("Row", "Col","FeatureNum", "ControlType","ProbeName","SystematicName"),path=path)

  #data preprocessing
  RG.bg <- backgroundCorrect(RG, method="minimum", offset=1)
  MA <- normalizeWithinArrays(RG.bg, method="loess")
  MA <- normalizeBetweenArrays(MA, method="Aquantile")
  #do not average replcicates
  #MA.avg <- avereps(MA, ID=MA$genes$SystematicName)

  #extract data from separate channels
  RG.pq <- RG.MA(MA)
  exptab <- cbind(RG.pq$R, RG.pq$G)
  colnames(exptab)=c(paste(MA$targets$Cy3,'R',sep='_'), paste(MA$targets$Cy5,'G',sep='_'))
  SampleIDs = colnames(exptab)
  annotCols = c('ProbeName','SystematicName','ProbeID','ProbeType','Feature')

  #add above BG information
  exptab_ab = cbind(RG$other$rIsWellAboveBG, RG$other$gIsWellAboveBG)
  colnames(exptab_ab)=paste0(colnames(exptab),'.AboveBG')
  exptabAnnot = cbind(MA$genes[,c('ProbeName','SystematicName')],exptab,exptab_ab)

  #extract probe details
  ids = strsplit2(exptabAnnot$ProbeName,'_')
  exptabAnnot$Feature = ids[,4]
  exptabAnnot$ProbeID = paste(ids[,1],ids[,2],ids[,3],sep="_")
  exptabAnnot$ProbeType = ids[,2]

  #extract only the factor-specific probes
  exptabAnnotSub = exptabAnnot[ids[,2] %in% c('GC','LN','An','G4','T7'),]
  exptabAnnotSub = exptabAnnotSub[order(exptabAnnotSub$ProbeID),]

  return(list(exptab = exptabAnnotSub, sampleIDs = SampleIDs, annotCols = annotCols))
}
