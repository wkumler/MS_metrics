
# Setup ----
library(xcms)
library(tidyverse)
library(RaMS)
library(data.table)
options(pillar.sigfig=7)

file_data <- read_csv("made_data/file_data.csv") %>%
  mutate(filename=basename(filename))
peak_bounds <- read_csv("made_data/peak_bounds.csv")

msnexp_filled <- readRDS("made_data/msnexp_filled.rds")
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
  mutate(filename=file_data$filename[sample])
classified_feats <- read_csv("made_data/classified_feats.csv") %>%
  select(feature, feat_class)



# Simple feature extraction - those provided by XCMS directly ----
simple_feats <- peak_data %>%
  group_by(feat_id) %>%
  summarise(mean_mz=unique(feat_mzmed), sd_ppm=sd(mz)/feat_mzmed,
            mean_rt=unique(feat_rtmed), sd_rt=sd(rt),
            mean_pw=mean(rtmax-rtmin), sd_pw=sd(rtmax-rtmin),
            mean_area=mean(into), med_area=median(into), IQR(into),
            log_mean_height=log10(mean(maxo)), log_sd_height=log10(sd(maxo)),
            across(c(sn, f, scale, scpos, lmin, lmax), mean, na.rm=TRUE),
            n_found=41-sum(is.na(intb)), 
            samps_found=1-sum(is.na(intb) & str_detect(filename, "Smp"))/24,
            stans_found=1-sum(is.na(intb) & str_detect(filename, "Std"))/10,
            blank_found=any(!is.na(intb) & str_detect(filename, "Blk")))

# Calculate peak shape metrics from EICs ----
msdata <- readRDS("made_data/msdata.rds")
eic_dt <- peak_bounds %>%
  mutate(rtmin=rtmin/60, rtmax=rtmax/60, filename=basename(filename)) %>%
  pmap_dfr(function(...){
    row_data <- list(...)
    msdata$EIC2[filename==row_data$filename
                ][mz%between%c(row_data$mzmin, row_data$mzmax)
                  ][rt%between%c(row_data$rtmin, row_data$rtmax)][
                    , feature:=row_data$feature][]
  })
qscoreCalculator <- function(rt, int){
  #Check for bogus EICs
  if(length(rt)<5){
    return(list(0, 0))
  }
  #Calculate where each rt would fall on a beta dist (accounts for missed scans)
  scaled_rts <- (rt-min(rt))/(max(rt)-min(rt))
  
  # Create a couple different skews and test fit
  maybe_skews <- c(2.5,3,4,5) #Add 7 to catch more multipeaks and more noise
  #Add 2 to catch very slopey peaks and more noise
  best_skew <- maybe_skews[which.max(sapply(maybe_skews, function(x){
    cor(dbeta(scaled_rts, shape1 = x, shape2 = 5), int)
  }))]
  perf_peak <- dbeta(scaled_rts, shape1 = best_skew, shape2 = 5)
  peak_cor <- cor(perf_peak, int)
  
  
  #Calculate the normalized residuals
  residuals <- int/max(int)-perf_peak/max(perf_peak)
  #Calculate the minimum SD, after normalizing for any shape discrepancy
  old_res_sd <- sd(residuals)
  norm_residuals <- diff(residuals)
  new_res_sd <- sd(norm_residuals)
  while(new_res_sd<old_res_sd){
    old_res_sd <- new_res_sd
    norm_residuals <- diff(residuals)
    new_res_sd <- sd(residuals)
  }
  #Calculate SNR
  SNR <- (max(int)-min(int))/sd(norm_residuals*max(int))
  #Return the quality score
  return(list(SNR=SNR, peak_cor=peak_cor))
}
peakshape_mets <- eic_dt %>%
  group_by(feature, filename) %>%
  summarise(qscores=list(qscoreCalculator(rt, int))) %>%
  unnest_wider(qscores) %>%
  summarise(med_SNR=median(SNR, na.rm=TRUE), 
            med_cor=median(peak_cor, na.rm=TRUE))
peakshape_mets %>%
  left_join(classified_feats) %>%
  mutate(med_SNR=cut(med_SNR, breaks = seq(0, 30, 6))) %>%
  ggplot() +
  geom_bar(aes(x=med_SNR, fill=feat_class), position = "fill")
peakshape_mets %>%
  left_join(classified_feats) %>%
  mutate(med_cor=cut(med_cor, breaks = seq(0, 1, 0.1))) %>%
  ggplot() +
  geom_bar(aes(x=med_cor, fill=feat_class), position = "fill")


scan_time_diff <- diff(sort(unique(msdata$EIC2[filename==filename[1]]$rt)))
hist(scan_time_diff)
scan_time <- mean(scan_time_diff)

n_scans <- eic_dt[, .N, .(filename, feature)]
hist(n_scans$N, breaks = 100)
med_missed_scans <- peak_bounds %>%
  mutate(filename=basename(filename)) %>%
  mutate(rtmin=rtmin/60, rtmax=rtmax/60) %>%
  mutate(expected_scans=round((rtmax-rtmin)/scan_time)) %>%
  left_join(n_scans) %>%
  # with(hist(expected_scans-N, breaks = 100))
  group_by(feature) %>%
  summarise(med_missed_scans=median(expected_scans-N)) %>%
  # with(hist(med_missed_scans, breaks = 100))
  ungroup()

rt_dt <- rtime(msnexp_filled) %>%
  as.data.frame() %>%
  setNames("rt") %>%
  rownames_to_column("filenum") %>%
  mutate(filenum=as.numeric(str_extract(filenum, "(?<=F)\\d+"))) %>%
  mutate(filename=basename(fileNames(msnexp_filled))[filenum]) %>%
  mutate(rt=round(rt/60, digits = 7)) %>%
  select(filename, rt) %>%
  as.data.table()
file_feat_n_missed <- eic_dt %>%
  mutate(rt=round(rt, digits = 7)) %>%
  group_by(filename, feature) %>%
  group_split() %>%
  pbapply::pblapply(function(eic){
    fname_i <- eic$filename[1]
    min_rt=min(eic$rt)
    max_rt=max(eic$rt)
    psb_rts <- rt_dt[rt%between%c(min_rt, max_rt)][filename==fname_i]
    merged_dt <- merge(unique(eic), psb_rts, all.y=TRUE)
    summarise(merged_dt, 
              filename=unique(filename), 
              feature=unique(feature),
              n_missed=sum(is.na(merged_dt$int)),
              n_scans=n()) %>%
      na.omit()
  }) %>%
  bind_rows()
med_missed_scans_2 <- file_feat_n_missed %>%
  group_by(feature) %>%
  summarise(med_missed_scans_2=median(n_missed, na.rm=TRUE))
med_missed_scans_2 %>%
  left_join(classified_feats) %>%
  mutate(med_missed_scans_2=cut(med_missed_scans_2, breaks = c(0, 1, 3, 5, 10, 20, 80))) %>%
  ggplot() +
  geom_bar(aes(x=med_missed_scans_2, fill=feat_class), position = "fill")


# Calculate presence/absence of isotope info ----



# Calculate DOE metrics ----
depth_diffs <- peak_data %>%
  select(feat_id, filename, into) %>%
  left_join(file_data) %>%
  filter(samp_type=="Smp") %>%
  nest(data=-feat_id) %>%
  mutate(t_pval=map_dbl(data, function(x){
    if(any(table(x$depth)<2))return(1)
    broom::tidy(t.test(x$into~x$depth))$p.value
  })) %>%
  select(feat_id, t_pval)
depth_diffs %>%
  left_join(classified_feats, by=c(feat_id="feature")) %>%
  # mutate(t_pval=cut(log10(t_pval), breaks = seq(-5, 0, 1))) %>%
  mutate(t_pval=cut(t_pval, breaks = seq(0, 1, 0.1))) %>%
  ggplot() +
  geom_bar(aes(x=t_pval, fill=feat_class), position = "fill")

blank_diffs <- peak_data %>%
  select(feat_id, filename, into) %>%
  left_join(file_data) %>%
  filter(samp_type%in%c("Smp", "Blk")) %>%
  select(feat_id, samp_type, into) %>%
  group_by(feat_id, samp_type) %>%
  summarise(mean_area=mean(into)) %>%
  pivot_wider(names_from = samp_type, values_from = mean_area) %>%
  mutate(smp_to_blk=Smp/Blk) %>%
  select(feat_id, smp_to_blk)
blank_diffs %>%
  left_join(classified_feats, by=c(feat_id="feature")) %>%
  mutate(smp_to_blk=cut(log10(smp_to_blk), breaks = -1:4)) %>%
  ggplot() +
  geom_bar(aes(x=smp_to_blk, fill=feat_class), position = "fill")



# Join feature classes and write out ----
features_extracted <- simple_feats %>%
  left_join(peakshape_mets, by=c(feat_id="feature")) %>%
  left_join(med_missed_scans, by=c(feat_id="feature")) %>%
  left_join(depth_diffs) %>%
  left_join(blank_diffs) %>%
  left_join(classified_feats, by=c(feat_id="feature"))
write.csv(features_extracted, "made_data/features_extracted.csv", row.names = FALSE)














library(plotly)
plot_ly(features_extracted, x=~mean_pw, y=~log10(sn), z=~n_found, color=~feat_class,
        type = "scatter3d", mode="markers")
plot_ly(features_extracted, x=~sd_rt, y=~log10(mean_area), z=~sd_mz, color=~feat_class,
        type = "scatter3d", mode="markers")
plot_ly(features_extracted, x=~n_found, y=~samps_found, z=~blank_found, color=~feat_class,
        type = "scatter3d", mode="markers")
features_extracted %>%
  plot_ly(x=~med_SNR, y=~med_cor, z=~log10(mean_area), 
          color=~feat_class, text=~feat_id,
          type = "scatter3d", mode="markers")
features_extracted %>%
  plot_ly(x=~sqrt(med_SNR), y=~med_cor^4, 
          color=~feat_class, text=~feat_id,
          type = "scatter", mode="markers")
