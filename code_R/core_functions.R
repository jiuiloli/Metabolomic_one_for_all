# this is core functions only
library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)
library(devtools)
library(ggbiplot)
library(data.table)
library(factoextra)
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
mzs <- mz(raw_data)

## Split the list by file
mzs_by_file <- split(mzs, f = fromFile(raw_data))

length(mzs_by_file)

bpis <- chromatogram(raw_data, aggregationFun = "max")
group_colors <- paste0(brewer.pal(7, "Set1")[1:7], "60")
names(group_colors) <- c("GG", "JS", "LX", "QC", 
                         "SY", "WC", "YF")

## Plot all chromatograms.
plot(bpis, col = group_colors[raw_data$sample_group])
# 
rtr <- c(616, 624)
mzr <- c(146.9, 147.1)
# extract the cromatogram
chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)
plot(chr_raw, col = group_colors[chr_raw$sample_group])

# peak detection using matchedfilter algorithm
cwp <- CentWaveParam(ppm = 100, peakwidth = c(3,10), mzdiff = 0.01, snthresh = 6)
xdata <- findChromPeaks(raw_data, param = cwp)
chr_ex <- chromatogram(xdata, mz = mzr, rt = rtr)
sample_colors <- group_colors[chr_ex$sample_group]
plot(chr_ex, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_ex)[, "sample"]],
     peakBg = NA)
par(mfrow = c(1,1), mar = c(4.1, 8.1, 4.1, 2.1))
plotChromPeakImage(xdata)
if(hasAdjustedRtime(xdata)>0){
  xdata <- dropAdjustedRtime(xdata)
}
plotChromPeaks(xdata, file = 3)
## Extract a list of per-sample peak intensities (in log2 scale)
par(mfrow = c(1, 1), mar = c(8.1, 4.1, 4.1, 2.1))
# bottom, left, top, right, default c(5.1, 4.1, 4.1, 2.1)
ints <- split(log2(chromPeaks(xdata)[, "into"]),
              f = chromPeaks(xdata)[, "sample"])
boxplot(ints, varwidth = TRUE, col = group_colors[xdata$sample_group],
        ylab = expression(log[2]~intensity), main = "Peak intensities",
        names = basename(fileNames(xdata)), las = 2)
grid(nx = NA, ny = NULL)




xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 1))

## Get the base peak chromatograms.
bpis_adj <- chromatogram(xdata, aggregationFun = "max")
par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
plot(bpis_adj, col = group_colors[bpis_adj$sample_group], peakType = "none")
## Plot also the difference of adjusted to raw retention time.
plotAdjustedRtime(xdata, col = group_colors[xdata$sample_group])



# alignment
for(i in c(5,4,3,2,1)){
  pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                          minFraction = 0.5, bw = i, minSamples = 1, binSize = 0.25,
                          maxFeatures = 100)
  xdata <- groupChromPeaks(xdata, param = pdp)
}

# normalizaiton has not been done!!!!
# pgp <- PeakGroupsParam(minFraction = 0.85, smooth = "linear", family = "gaussian")
# xdata <- adjustRtime(xdata, param = pgp)

t1 <- as.data.frame(apply(featureValues(xdata, filled = FALSE), MARGIN = 2,
      FUN = function(z) sum(is.na(z))))
xdata <- fillChromPeaks(xdata)
t1<- as.data.frame(apply(featureValues(xdata, filled = TRUE), MARGIN = 2,
                         FUN = function(z) sum(is.na(z))))

summary <- as.data.frame(featureSummary(xdata, group = xdata$sample_group))
feature_chroms <- featureChromatograms(xdata, features = 1:4)

feature_chroms
plot(feature_chroms)
full_data <- as.data.frame(t(ft_ints))
full_data <- setDT(full_data, keep.rownames = TRUE)[]
ft_ints <- log2(featureValues(xdata, value = "into"))
data <- as.data.frame(t(na.omit(ft_ints)))
data <- setDT(data, keep.rownames = TRUE)[]
data$rn <- substr(data$rn, 0, 2)
pca <-prcomp(data[,-1])
fviz_pca_ind(pca, geom.ind = "point", 
             pointshape = 21,
             pointsize = 2,
             fill.ind = data$rn,
             col.ind = "black",
             palette = "jco", 
             addEllipses = TRUE,
             ellipse.level = 0.95,
             invisible ="quali",
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "",) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom",
        panel.border=element_rect(fill=NA), axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20, face = "bold"))


par(mfrow = c(2, 1), mar = c(5.1, 4.1, 4.1, 2.1))
## Plot the raw data
plot(chr_raw, col = group_colors[chr_raw$sample_group])

## Extract the chromatogram from the adjusted object
chr_adj <- chromatogram(xdata, rt = rtr, mz = mzr)
plot(chr_adj, col = group_colors[chr_raw$sample_group], peakType = "none")
