# This script takes in the mzML files supplied in the mzMLs/ directory
# and performs peakpicking and other xcms functions on them before exporting
# the peak/feature boundaries in a useful format into made_data/

# The full script takes a minute or two on my local machine
# Longer for the FT2040 - the msdata extraction alone takes ~20 mins each
# Even longer for MS3000 - msdata takes ~45 mins each.

# This is the first script in the pipeline

# Setup ----
library(tidyverse)
library(xcms)
library(RaMS)
options(pillar.sigfig=7)

# dataset_version <- "FT2040"
# dataset_version <- "MS3000"
dataset_version <- "CultureData"
# dataset_version <- "Pttime"

output_folder <- paste0("made_data_", dataset_version, "/")
mzML_files <- list.files(paste0(output_folder, "mzMLs/"), full.names=TRUE)
prefilter_versioned <- c(3, 1e6)

if(dataset_version%in%c("FT2040", "MS3000")){
  file_data <- data.frame(filename=mzML_files) %>%
    mutate(samp_type=str_extract(filename, "Blk|Smp|Std|Poo")) %>%
    mutate(depth=str_extract(filename, "25m|DCM|175m|15m")) %>%
    mutate(depth=ifelse(is.na(depth), "", depth)) %>%
    mutate(samp_group=factor(paste0(depth, samp_type), levels=c("Blk", "25mSmp", "DCMSmp", "175mSmp", "15mSmp", "Std", "Poo"))) %>%
    mutate(col=alpha(c("red", "blue", "green", "purple", "blue", "black", "#008080"), 0.8)[samp_group]) %>%
    mutate(lwd=c(2, 1, 1, 1, 1, 2)[samp_group])
} else if(dataset_version=="CultureData"){
  organism_vec <- c("Cy", "Fc", "Np", "Pc55x", "Pt", "To", "Tp", "1314", "1545", 
                    "1771", "2021", "3430", "449", "1314P", "7803", "8102", "8501", 
                    "As9601", "MED4", "Nat", "116-1", "2090", "371", "Ca", "Cr", 
                    "Cs", "P3", "P5-5", "DSS3", "Och", "SA11", "SA16", "SA36", "SA44", "SA42", 
                    "SA48", "SA53", "SA55", "SA59", "SA60", "SA7", "SAR11", "SCM1"
  )
  file_data <- data.frame(filename=mzML_files) %>%
    mutate(samp_type=str_extract(filename, "Blk|Smp|Std|Poo")) %>%
    mutate(organism=str_extract(filename, paste(organism_vec, collapse = "|"))) %>%
    mutate(organism=ifelse(organism=="SA42", "SA44", organism)) %>%
    mutate(organism=ifelse(samp_type=="Smp", organism, NA)) %>%
    mutate(org_class=case_when(
      samp_type%in%c("Blk", "Std") ~ NA,
      samp_type=="Poo" ~ str_extract(filename, "Diatom|DinosGreens|Cyanos|Bacteria|Haptophytes"),
      organism%in%c("Cy", "Fc", "Np", "Pc55x", "Pt", "To", "Tp") ~ "Diatom",
      organism%in%c("1314", "1545", "1771", "2021", "449", "3430") ~ "DinosGreens",
      organism%in%c("DSS3", "Och", paste0("SA", c(7, 11, 16, 36, 44, 48, 53, 55, 59, 60, 67)),
                    "SAR11", "SCM1") ~ "Bacteria",
      organism%in%c("1314P", "7803", "8102", "8501", "As9601", "MED4", "Nat") ~ "Cyanos",
      organism%in%c("116-1", "2090", "371", "Ca", "Cr", "Cs", "P3", "P5-5") ~ "Haptos",
      TRUE~"UNKNOWN" 
    )) %>%
    mutate(samp_group=str_remove(basename(filename), "_?(St[1-4]_pos|Exp[1-3]_pos|MediaBlk.*|[A-I]{1,2}_pos|\\d)\\.mzML$"))
} else if(dataset_version=="Pttime"){
  file_data <- data.frame(filename=mzML_files) %>%
    mutate(samp_id=str_remove(basename(filename), "20190429_JJ_VB_BioSFA_Pttime_1_QE144_Ag68377-924_USHXG01160_POS_MSMS-v2_")) %>%
    mutate(samp_type=str_extract(samp_id, "exudate|pellet|extr")) %>%
    mutate(timepoint=str_extract(samp_id, "\\d+d")) %>%
    mutate(timepoint=ifelse(is.na(timepoint), "ctrl", timepoint)) %>%
    mutate(samp_group=paste(samp_type, timepoint, sep = "-")) %>%
    arrange(samp_group)
} else {
  stop("Dataset not recognized!")
}

# XCMS things ----
register(BPPARAM = SerialParam(progressbar = TRUE))
msnexp <- readMSData(
  files = file_data$filename, 
  pdata = new("NAnnotatedDataFrame", file_data), 
  msLevel. = 1, 
  mode = "onDisk"
)

register(BPPARAM = SnowParam(workers = parallel::detectCores()-1, 
                             tasks = nrow(file_data), progressbar = TRUE))
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
if(dataset_version%in%c("FT2040", "MS3000", "Pttime")){
  obp <- ObiwarpParam(
    binSize = 0.1, 
    centerSample = round(nrow(file_data)/2), 
    response = 1, 
    distFun = "cor_opt"
  )
  msnexp_rtcor <- adjustRtime(msnexp_withpeaks, obp)
} else {
  msnexp_rtcor <- msnexp_withpeaks
}

if(dataset_version%in%c("FT2040", "MS3000", "Pttime")){
  pdp <- PeakDensityParam(
    sampleGroups = file_data$samp_group, 
    bw = 12, 
    minFraction = 0.1, 
    binSize = 0.001, 
    minSamples = 2
  )
} else if(dataset_version=="CultureData"){
  pdp <- PeakDensityParam(
    sampleGroups = file_data$samp_group, 
    bw = 12, 
    minFraction = 0.4, 
    binSize = 0.001, 
    minSamples = 3
  )
}
msnexp_grouped <- groupChromPeaks(msnexp_rtcor, pdp)

register(BPPARAM = SnowParam(workers = parallel::detectCores()-1, 
                             tasks = nrow(file_data), progressbar = TRUE))
fpp <- FillChromPeaksParam(ppm = 2.5)
msnexp_filled <- fillChromPeaks(msnexp_grouped, fpp)
register(BPPARAM = SerialParam(progressbar = TRUE))


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
  mutate(filename=basename(fileNames(msnexp_filled)[as.numeric(str_extract(scanid, "\\d+"))])) %>%
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
  summarize(min_mz=min(mzmin), max_mz=max(mzmax), 
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
