# This script takes in the mzML files supplied in the mzMLs/ directory
# and performs peakpicking and other xcms functions on them before exporting
# the peak/feature boundaries in a useful format into made_data/

# The full script takes a minute or two on my local machine

# Setup ----
library(tidyverse)
library(xcms)
library(RaMS)
options(pillar.sigfig=7)

mzML_files <- list.files("mzMLs/", full.names=TRUE)

file_data <- data.frame(filename=mzML_files) %>%
  mutate(samp_type=str_extract(filename, "Blk|Smp|Std|Poo")) %>%
  mutate(depth=str_extract(filename, "25m|DCM")) %>%
  mutate(depth=ifelse(is.na(depth), "", depth)) %>%
  mutate(colid=factor(paste0(depth, samp_type), levels=c("Blk", "25mSmp", "DCMSmp", "Std", "Poo"))) %>%
  mutate(col=alpha(c("red", "blue", "green", "black", "#008080"), 0.8)[colid]) %>%
  mutate(lwd=c(2, 1, 1, 1, 2)[colid])

# XCMS things ----
register(BPPARAM = SnowParam(workers = 6, tasks = nrow(file_data), progressbar = TRUE))
register(BPPARAM = SerialParam(progressbar = TRUE))
msnexp <- readMSData(
  files = file_data$filename, 
  pdata = new("NAnnotatedDataFrame", file_data), 
  msLevel. = 1, 
  mode = "onDisk"
)
cwp <- CentWaveParam(
  ppm = 5, 
  peakwidth = c(20, 80), 
  prefilter = c(5, 1e7), 
  snthresh = 0, 
  verboseColumns = TRUE, 
  extendLengthMSW = TRUE, 
  integrate = 2
)
msnexp_withpeaks <- findChromPeaks(msnexp, cwp)

obp <- ObiwarpParam(
  binSize = 0.1, 
  centerSample = round(nrow(file_data)/2), 
  response = 1, 
  distFun = "cor_opt"
)
msnexp_rtcor <- adjustRtime(msnexp_withpeaks, obp)

pdp <- PeakDensityParam(
  sampleGroups = file_data$colid, 
  bw = 12, 
  minFraction = 0.1, 
  binSize = 0.001, 
  minSamples = 2
)
msnexp_grouped <- groupChromPeaks(msnexp_rtcor, pdp)

fpp <- FillChromPeaksParam(ppm = 2.5)
msnexp_filled <- fillChromPeaks(msnexp_grouped, fpp)



# Organizing into useful format ----
init_rts <- msnexp_filled %>%
  dropAdjustedRtime() %>%
  rtime() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  set_names(c("scanid", "init_rt"))
rt_corrections <- msnexp_filled %>%
  rtime() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  set_names(c("scanid", "new_rt")) %>%
  mutate(filename=mzML_files[as.numeric(str_extract(scanid, "\\d+"))]) %>%
  left_join(init_rts, by="scanid") %>%
  mutate(new_rt=new_rt/60, rt=init_rt/60) %>%
  select(-init_rt, -scanid)

peak_data_long <- msnexp_filled %>%
  chromPeaks() %>%
  as.data.frame() %>%
  rownames_to_column()
feature_bounds <- msnexp_filled %>%
  featureDefinitions() %>%
  as.data.frame() %>%
  select(mzmed, rtmed, npeaks, peakidx) %>%
  rownames_to_column("id") %>%
  unnest_longer(peakidx) %>%
  rename_with(~paste0("feat_", .x)) %>%
  mutate(peak_data=peak_data_long[feat_peakidx,]) %>%
  unnest_wider(peak_data) %>%
  mutate(filename=mzML_files[sample]) %>%
  select(feature=feat_id, filename, mzmin, mzmax, rtmin, rtmax)

# Exporting ----
write.csv(file_data, "made_data/file_data.csv", row.names = FALSE)
saveRDS(msnexp_filled, "made_data/msnexp_filled.rds")
write.csv(rt_corrections, "made_data/rt_corrections.csv", row.names = FALSE)
write.csv(feature_bounds, "made_data/feature_bounds.csv", row.names = FALSE)
