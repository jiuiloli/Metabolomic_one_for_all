library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)

## Get the full path to the CDF files
cdfs <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
            recursive = TRUE)
## Create a phenodata data.frame
pd <- data.frame(sample_name = sub(basename(cdfs), pattern = ".CDF",
                                   replacement = "", fixed = TRUE),
                 sample_group = c(rep("KO", 6), rep("WT", 6)),
                 stringsAsFactors = FALSE) 
raw_data <- readMSData(files = cdfs, pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk")
cwp <- CentWaveParam(peakwidth = c(30, 80), noise = 1000)
xdata <- findChromPeaks(raw_data, param = cwp) 
plotChromPeakImage(xdata)
typeof(xdata)
