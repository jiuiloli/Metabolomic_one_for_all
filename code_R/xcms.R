setwd("~/Desktop/GFSC/Projects/Rice Authentication Project/GC_data_msZXL/")
 
library(xcms)
cdfpath <- file.path('.')  
cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE, pattern='*.mzXML')
xset <- xcmsSet(cdffiles, method = "matchedFilter", fwhm=3, snthresh=10, max=100, step=0.25, steps=2, mzdiff=0.5)
xset <- group(xset)
# pdf(paste('retentiontime_correction','.pdf',sep=''))
xset2 <- retcor(xset, missing = 1, extra = 1, smooth = "linear", family = "gaussian", plottype = "mdevden")
png(file = "retentiontime_correction.png")
dev.off()
xset2 <- group(xset2, method="density", mzwid= 0.25, minfrac=0.5, bw = 3)
xset3 <- fillPeaks(xset2)
peak_out <- peakTable(xset3, filebase="peakList")
rt_minutes <- peak_out$rt/60
rtmin_minutes <- peak_out$rtmin/60
rtmax_minutes <- peak_out$rtmax/60
mz_rt <- paste("mz", signif(peak_out$mz,digits = 4), "_rt", signif(rt_minutes,digits = 4), sep = "")
peak_final_out <- cbind(mz_rt, peak_out[, 1:3], rt_minutes, rtmin_minutes, rtmax_minutes, peak_out[, 7:ncol(peak_out)])
peak_final_out <- peak_final_out[order(rt_minutes), ]
write.csv(peak_final_out, file='Metabolite Table.csv')
warnings()
warnings_out <- warnings()
 