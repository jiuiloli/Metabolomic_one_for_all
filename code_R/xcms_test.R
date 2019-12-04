library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)

# -----------this line is buggy, don't use it for now
#Find out how many cores are available (if you don't already know)
# cores<-detectCores()
# #Create cluster with desired number of cores, leave one open for the machine         
# #core processes
# cl <- makeCluster(cores[1]-1)
# #Register cluster
# registerDoParallel(cl)



setwd("~/Documents/GitHub/Metabolomic_one_for_all/data_mzXML/")
######################################## data import
## Get the full path to the CDF files
mzXMLs <- dir(("~/Documents/GitHub/Metabolomic_one_for_all/data_mzXML/"),
              full.names = TRUE, recursive = TRUE)
## Create a phenodata data.frame
md <- data.frame(sample_name = sub(basename(mzXMLs), pattern = ".mzXML",
                                   replacement = "", fixed = TRUE),
                 sample_group = c(rep("GG", 16), rep("JS", 20), rep("LX", 19), rep("QC", 18),
                                  rep("SY", 20), rep("WC", 25), rep("YF", 33)),
                 stringsAsFactors = FALSE)

raw_data_originial <- readMSData(files = mzXMLs, pdata = new("NAnnotatedDataFrame", md),
                       mode = "onDisk")

raw_data <- raw_data_originial
raw_data <- filterRt(raw_data_originial, rt = c(0, 1500))


####################################### initial data inspection

## these are just checking how many data points out there
# head(rtime(raw_data))
# mzs <- mz(raw_data)
# mzs_by_file <- split(mzs, f = fromFile(raw_data))
# length(mzs_by_file)



## Get the base peak chromatograms. This reads data from the files.
#returning spectrum the maximal intensity, or use sum for total ion chromatogram
bpis <- chromatogram(raw_data, aggregationFun = "max")

# give different group some colors
group_colors <- paste0(brewer.pal(7, "Set1")[1:7], "60")
names(group_colors) <- c("GG", "JS", "LX", "QC", "SY", "WC", "YF")

# this give all chromatographys
plot(bpis, col = group_colors[raw_data$sample_group])
# this one just gives first 2 chromatography
# plot(bpis[1,2], col = group_colors[raw_data1$sample_group])
# Plot all chromatograms.
# plot(bpis, col = group_colors[raw_data1$sample_group])

# -------------this doesn't work out really well -> may be only a small fraction needed?
# work_data <- filterRt(raw_data1, rt = c(700,1200))
# tc <- split(tic(work_data), f = fromFile(work_data))
# boxplot(tc, col = group_colors[work_data$sample_group],
#         ylab = "intensity", main = "Total ion current")


## Bin the BPC
bpis_bin <- bin(bpis, binSize = 2)

## Calculate correlation on the log2 transformed base peak intensities
cormat <- cor((do.call(cbind, lapply(bpis_bin, intensity))))
colnames(cormat) <- rownames(cormat) <- raw_data$sample_name

## Define which phenodata columns should be highlighted in the plot
ann <- data.frame(group = raw_data$sample_group)
rownames(ann) <- raw_data$sample_name

## Perform the cluster analysis
pheatmap(cormat, annotation = ann,
         annotation_color = list(group = group_colors))



######################################## chromatographic peak detection
# peak width optimization
# rtr usually is the 4x normal width of the peak
rtr <- c(616, 624)
mzr <- c(146.9, 147.1)
# extract the cromatogram
chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)
par(mfrow = c(2,1))
plot(chr_raw, col = group_colors[chr_raw$sample_group])
# peak width is 5 second
# reference value: fwhm 3


# ppm alignment
png(file="drink.png", width = 4000, height = 4000, bg = "white")
raw_data %>%
  filterRt(rt = rtr) %>%
  filterMz(mz = mzr) %>%
  plot(type = "XIC")
dev.off()
# ppm of choice is about 100
# reference value is 0.1 from website
xchr <- findChromPeaks(chr_raw, param = CentWaveParam(snthresh = 2))
# reference value for snthresh is 2
sample_colors <- group_colors[xchr$sample_group]
plot(xchr, col = sample_colors,
     peakBg = sample_colors[chromPeaks(xchr)[, "column"]])
# implementing centwave algorithm for peak detection
mfp <- MatchedFilterParam(fwhm = 3, snthresh = 2, mzdiff = 0.01)
# cwp <- CentWaveParam(peakwidth = c(5,10), snthresh = 6, ppm = 100, mzdiff = 0.01)
# xdata <- findChromPeaks(raw_data, param = cwp)
xdata <- findChromPeaks(raw_data, param = mfp)
chr_ex <- chromatogram(xdata, mz = mzr, rt = rtr)
data <- as.data.frame(chromPeaks(chr_ex))
pander(chromPeaks(xdata, mz = mzr, rt = rtr),
       caption = paste("Identified chromatographic peaks in a selected ","m/z and retention time range."))

# check if there is any peak filled up
# which(peakdata$is_filled == "TRUE")

# get some part chopped

# look for identified peaks in selected samples
plotChromPeaks(xdata, file =1)

# summary

summary_fun <- function(z)
  c(peak_count = nrow(z), rt = quantile(z[, "rtmax"] - z[, "rtmin"]))

T <- lapply(split.data.frame(
  chromPeaks(xdata), f = chromPeaks(xdata)[, "sample"]),
  FUN = summary_fun)
T <- do.call(rbind, T)
rownames(T) <- basename(fileNames(xdata))
pandoc.table(T,
             caption = paste0("Summary statistics on identified chromatographic",
                              " peaks. Shown are number of identified peaks per",
                              " sample and widths/duration of chromatographic ",
                              "peaks."))

# general heatplot
plotChromPeakImage(xdata)

sample_colors <- group_colors[chr_ex$sample_group]
plot(chr_ex, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_ex)[, "sample"]],
     peakBg = NA)




ints <- split(log2(chromPeaks(xdata)[, "into"]),
              f = chromPeaks(xdata)[, "sample"])
boxplot(ints, varwidth = TRUE, col = group_colors[xdata$sample_group],
        ylab = expression(log[2]~intensity), main = "Peak intensities", names = basename(fileNames(xdata)), las = 2)
grid(nx = NA, ny = NULL)
par(mfrow = c(1,1))

######################################## alignment
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 1))
# reference value is 1
head(adjustedRtime(xdata))
bpis_adj <- chromatogram(xdata, aggregationFun = "max")
par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
plot(bpis_adj, col = group_colors[bpis_adj$sample_group], peakType = "none")
plotAdjustedRtime(xdata, col = group_colors[xdata$sample_group])


# visual check
par(mfrow = c(2,1))
plot(chr_raw, col = group_colors[chr_raw$sample_group])

## Extract the chromatogram from the adjusted object and compare it with the original one
chr_adj <- chromatogram(xdata, rt = rtr, mz = mzr)
plot(chr_adj, col = group_colors[chr_raw$sample_group], peakType = "none")

if(hasAdjustedRtime(xdata) > 0){
  xdata <- applyAdjustedRtime(xdata)
}
hasAdjustedRtime(xdata)
######################################## correspondence
mzr <- c(147.0, 147.1)
chr_mzr<- chromatogram(xdata, mz = mzr, rt = c(900,1100))
par(mfrow = c(3,1), mar = c(1,4,1,0.5))
cols <- group_colors[chr_mzr$sample_group]
plot(chr_mzr, col = cols, xaxt = "n", xlab = "",
     peakBg = sample_colors[chromPeaks(chr_mzr)[, "sample"]], xlim = c(900,1100))

pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                        minFraction = 0.1, bw = 2)
plotChromPeakDensity(xdata, mz = mzr, col = cols, param = pdp,
                     pch = 16, xlim = c(900, 1100))
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                        minFraction = 0.1, bw = )
plotChromPeakDensity(xdata, mz = mzr, col = cols, param = pdp,
                     pch = 16, xlim = c(900, 1100))
# lets try 1 first
## Perform the correspondence
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                        minFraction = 0.2, bw = 2)
xdata <- groupChromPeaks(xdata, param = pdp)
data <- featureValues(xdata, value = "into")
# fill up the missing values
xdata <- fillChromPeaks(xdata)
data <- featureValues(xdata, value = "into")
# check for how many missing values were there and are there
apply(featureValues(xdata, filled = FALSE), MARGIN = 2,
      FUN = function(z) sum(is.na(z)))
apply(featureValues(xdata), MARGIN = 2,
      FUN = function(z) sum(is.na(z)))
# fruther data processing and analysis







# additionals