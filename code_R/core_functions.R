# this is core functions only
xdata <- findChromPeaks(raw_data, param = mfp)
pander(chromPeaks(xdata, mz = mzr, rt = rtr),
       caption = paste("Identified chromatographic peaks in a selected ","m/z and retention time range."))
if(hasAdjustedRtime(xdata)>0){
  xdata <- dropAdjustedRtime(xdata)
}
pdp<-PeakDensityParam(sampleGroups = xdata$sample_group,
                      minFraction = 0.8)
xdata <- groupChromPeaks(xdata, param = pdp)

pgp <- PeakGroupsParam(minFraction = 0.9, span = 0.4)
xdata <- adjustRtime(xdata, param = pgp)
par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 2))
plotAdjustedRtime(xdata, col = group_colors[xdata$sample_group],
                  peakGroupsCol = "grey", peakGroupsPch = 1)
par(mfrow = c(2, 1))
## Plot the raw data
plot(chr_raw, col = group_colors[chr_raw$sample_group])

## Extract the chromatogram from the adjusted object
chr_adj <- chromatogram(xdata, rt = rtr, mz = mzr)
plot(chr_adj, col = group_colors[chr_raw$sample_group], peakType = "none")
par(mfrow = c(1, 1))


# alignment
par(mfrow = c(1,1))

mzr <- c(147.0, 147.1)
chr_mzr<- chromatogram(xdata, mz = mzr, rt = c(900,1100))
par(mfrow = c(3,1), mar = c(1,4,1,0.5))
cols <- group_colors[chr_mzr$sample_group]
plot(chr_mzr, col = cols, xaxt = "n", xlab = "",
     peakBg = sample_colors[chromPeaks(chr_mzr)[, "sample"]], xlim = c(900,1100))

pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                        minFraction = 0.4, bw = 2)
plotChromPeakDensity(xdata, mz = mzr, col = cols, param = pdp,
                     pch = 16, xlim = c(900, 1100))
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                        minFraction = 0.4, bw = 3)
plotChromPeakDensity(xdata, mz = mzr, col = cols, param = pdp,
                     pch = 16, xlim = c(900, 1100))
# let's use bw = 2

pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                        minFraction = 0.4, bw = 3)
xdata <- groupChromPeaks(xdata, param = pdp)
apply(featureValues(xdata, filled = FALSE), MARGIN = 2,
      FUN = function(z) sum(is.na(z)))
xdata <- fillChromPeaks(xdata)
apply(featureValues(xdata), MARGIN = 2,
      FUN = function(z) sum(is.na(z)))
feature_chroms <- featureChromatograms(xdata, features = 1:4)
plot(feature_chroms, col = sample_colors,
     peakBg = sample_colors[chromPeaks(feature_chroms)[, "sample"]])
data <- log2(featureValues(xdata, value = "into"))
