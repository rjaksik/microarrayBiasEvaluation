
<!-- README.md is generated from README.Rmd. Please edit that file -->

# microarrayBiasEvaluation

<!-- badges: start -->
<!-- badges: end -->

Functions used to analyze custom Agilent microarrays to determine the
extent of specific bias sources

## Instalation:

The package requires R 3.5.0 or later

``` r
install.packages("devtools")  
devtools::install_github("rjaksik/microarrayBiasEvaluation")
```

## Example use:

``` r
library('microarrayBiasEvaluation')

path='/RawData'
targets_file = "sample_desc.txt"
expData = procCustArray(targets_file,path)

probeData = read.table(system.file("resources", "A-MTAB-703_probe_details.txt", package = "microarrayBiasEvaluation"),header=T,sep='\t')
plotGCcontentBias(expData$exptab, 'CustMicro_GC_DistExprs_v2.pdf')
plotLengthBias(expData,'CustMicro_Len_DistExprs_v1.pdf')
plotAnBias(expData,probeData,'CustMicro_An')
plotT7Bias(expData,probeData,'CustMicro_T7')
plotG4Bias(expData,probeData,'CustMicro_G4')
```

## Acknowledgment

The creation of this software was supported by the Polish National Science Centre, grant No. 2016/23/D/ST7/03665.
