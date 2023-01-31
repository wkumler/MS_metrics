# This script takes in the mzML files supplied in the mzMLs/ directory
# and performs peakpicking and other xcms functions on them before exporting
# the peak/feature boundaries in a useful format into made_data/

# The full script takes a minute or two on my local machine
# Longer for the FT2040 - the msdata extraction alone takes ~20 mins each

# This is the first script in the pipeline

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

# dataset_version <- "FT350"
dataset_version <- "FT2040"
if(dataset_version=="FT350"){
  prefilter_versioned <- c(5, 1e7)
} else if(dataset_version=="FT2040") {
  prefilter_versioned <- c(3, 1e6)
} else {
  stop(paste("Version", dataset_version, "not yet supported!"))
}
output_folder <- paste0("made_data_", dataset_version, "/")
if(!dir.exists(output_folder))dir.create(output_folder)

# XCMS things ----
register(BPPARAM = SerialParam(progressbar = TRUE))
msnexp <- readMSData(
  files = file_data$filename, 
  pdata = new("NAnnotatedDataFrame", file_data), 
  msLevel. = 1, 
  mode = "onDisk"
)

register(BPPARAM = SnowParam(workers = 6, tasks = nrow(file_data), progressbar = TRUE))
cwp <- CentWaveParam(
  ppm = 5, 
  peakwidth = c(20, 80), 
  prefilter = prefilter_versioned, 
  snthresh = 0, 
  verboseColumns = TRUE, 
  extendLengthMSW = TRUE, 
  integrate = 2
)
msnexp_withpeaks <- findChromPeaks(msnexp, cwp)

register(BPPARAM = SerialParam(progressbar = TRUE))
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



# Extracting corrected retention times for later use ----
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
  mutate(filename=basename(mzML_files[as.numeric(str_extract(scanid, "\\d+"))])) %>%
  left_join(init_rts, by="scanid") %>% 
  mutate(new_rt=round(new_rt/60, 10), rt=round(init_rt/60, 10)) %>%
  select(-init_rt, -scanid)
# rt_corrections %>%
#   distinct(rt, new_rt, filename) %>%
#   mutate(diff=rt-new_rt) %>%
#   left_join(file_data) %>%
#   ggplot() +
#   geom_line(aes(x=rt, y=diff, group=filename, color=col)) +
#   scale_color_identity()

# Organizing into useful format ----
peak_data_long <- msnexp_filled %>%
  chromPeaks() %>%
  as.data.frame() %>%
  rownames_to_column()
peak_data <- msnexp_filled %>%
  featureDefinitions() %>%
  as.data.frame() %>%
  select(mzmed, rtmed, npeaks, peakidx) %>%
  rownames_to_column("id") %>%
  unnest_longer(peakidx) %>%
  rename_with(~paste0("feat_", .x)) %>%
  mutate(peak_data=peak_data_long[feat_peakidx,]) %>%
  unnest_wider(peak_data) %>%
  mutate(filename=mzML_files[sample])
peak_bounds <- peak_data %>%
  select(feature=feat_id, filename, mzmin, mzmax, rtmin, rtmax)

# Some visualization ----
peak_data %>%
  group_by(feat_id) %>%
  summarize(min_mz=min(mzmin), max_mz=max(mzmin), 
            min_rt=min(rtmin), max_rt=max(rtmax),
            max_int=max(maxo), max_area=max(into)) %>%
  mutate(int_group=cut(max_int, breaks = 10^c(3:9))) %>%
  ggplot() +
  geom_rect(aes(xmin=min_rt/60, xmax=max_rt/60, ymin=min_mz, ymax=max_mz, color=int_group),
            linewidth=1, fill=NA) +
  scale_color_viridis_d(name="Max feature\nintensity") +
  scale_x_continuous(name="Retention time (min)", breaks = seq(0, 1500, 120), 
                     sec.axis = sec_axis(~.*60, name = "Retention time (s)"), ) +
  scale_y_continuous(breaks = seq(60, 720, 30), name = "m/z ratio") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.75), legend.background = element_rect(fill = NA)) +
  ggtitle("Overall distributions of features")
peak_data %>%
  filter(feat_mzmed%between%pmppm(118.0865)) %>%
  ggplot() +
  geom_rect(aes(xmin=rtmin/60, xmax=rtmax/60, ymin=mzmin, ymax=mzmax, color=feat_id),
            linewidth=1, fill=NA)

peak_data %>%
  filter(feat_mzmed%between%pmppm(118.0865)) %>%
  ggplot() +
  geom_bar(aes(x=filename)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
peak_data %>%
  filter(feat_mzmed%between%pmppm(118.0865)) %>%
  ggplot() +
  geom_bar(aes(x=filename)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


# Exporting ----
write.csv(file_data, paste0(output_folder, "file_data.csv"), row.names = FALSE)
saveRDS(msnexp_filled, paste0(output_folder, "msnexp_filled.rds"))
write.csv(rt_corrections, paste0(output_folder, "rt_corrections.csv"), row.names = FALSE)
write.csv(peak_bounds, paste0(output_folder, "peak_bounds.csv"), row.names = FALSE)

# Extract relevant msdata from each for export and use later ----
feature_df <- peak_bounds %>%
  group_by(feature) %>%
  summarize(min_mz=min(mzmin), max_mz=max(mzmin), 
            min_rt=min(rtmin), max_rt=max(rtmax)) %>%
  mutate(mean_mz=(min_mz+max_mz)/2)

msdata <- file_data$filename %>%
  grabMSdata(verbosity = 1, grab_what = "EIC",
             mz=feature_df$mean_mz, ppm = 50)
msdata$EIC2 <- msdata$EIC %>%
  mutate(rt=round(rt, 10)) %>%
  left_join(rt_corrections, by=c("filename", "rt")) %>%
  select(rt=new_rt, mz, int, filename)
saveRDS(msdata, file = paste0(output_folder, "msdata.rds"))

msdata_isoc <- file_data$filename %>%
  grabMSdata(verbosity = 1, grab_what = "EIC",
             mz=feature_df$mean_mz+1.003355, ppm = 50)
msdata_isoc$EIC2 <- msdata_isoc$EIC %>%
  mutate(rt=round(rt, 10)) %>%
  left_join(rt_corrections, by=c("filename", "rt")) %>%
  select(rt=new_rt, mz, int, filename)
saveRDS(msdata_isoc, file = paste0(output_folder, "msdata_isoc.rds"))
