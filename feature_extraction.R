# This script extracts ML features from the peakpicking output and the
# raw data itself. These features include peak-specific features such
# as signal-to-noise, peak shape, etc. as well as isotope-specific
# features such as the similarity in peak shape between the base peak and
# the C13 peak as well as feature-specific features such as the difference
# between sample types (via t-test) and the difference between sample area
# and blank area

# This script requires access to the raw files (in mzMLs/), the xcms output
# (in made_data/, made by peakpicking_and_prep.R), and the classified 
# features (in made_data/, made by training_tool.R (although only for
# the very last step to make the output data frame look pretty)).

# This is the third script in the pipeline

# Setup ----
library(xcms)
library(tidyverse)
library(RaMS)
library(data.table)
options(pillar.sigfig=7)

trapz <- function(x, y) {
  m <- length(x)
  xp <- c(x, x[m:1])
  yp <- c(numeric(m), y[m:1])
  n <- 2*m
  p1 <- sum(xp[1:(n-1)]*yp[2:n]) + xp[n]*yp[1]
  p2 <- sum(xp[2:n]*yp[1:(n-1)]) + xp[1]*yp[n]
  
  return(0.5*(p1-p2))
}

# dataset_version <- "FT350"
# dataset_version <- "FT2040"
dataset_version <- "MS3000"
output_folder <- paste0("made_data_", dataset_version, "/")


file_data <- read_csv(paste0(output_folder, "file_data.csv")) %>%
  mutate(filename=basename(filename))
peak_bounds <- read_csv(paste0(output_folder, "peak_bounds.csv"))

msnexp_filled <- readRDS(paste0(output_folder, "msnexp_filled.rds"))
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


# Simple feature extraction - those provided by XCMS directly ----
n_files <- nrow(file_data)
n_samps <- sum(file_data$samp_type=="Smp")
n_stans <- sum(file_data$samp_type=="Std")
simple_feats <- peak_data %>%
  group_by(feat_id) %>%
  mutate(sn=log10(sn)) %>%
  summarise(mean_mz=unique(feat_mzmed), sd_ppm=sd(mz/feat_mzmed),
            mean_rt=unique(feat_rtmed), sd_rt=sd(rt),
            mean_pw=mean(rtmax-rtmin), sd_pw=sd(rtmax-rtmin),
            log_mean_height=log10(mean(maxo)), log_sd_height=log10(sd(maxo)),
            across(c(sn, f, scale, lmin), function(x)mean(x, na.rm=TRUE)),
            feat_npeaks=unique(feat_npeaks)/n_files,
            n_found=(sum(is.na(intb))-n_files)/n_files,
            samps_found=1-sum(is.na(intb) & str_detect(filename, "Smp"))/n_samps,
            stans_found=1-sum(is.na(intb) & str_detect(filename, "Std"))/n_stans,
            blank_found=any(!is.na(intb) & str_detect(filename, "Blk"))
            )
write.csv(simple_feats, paste0(output_folder, "simple_feats.csv"), row.names = FALSE)

# Calculate peak shape metrics from EICs ----
msdata <- readRDS(paste0(output_folder, "msdata.rds"))
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
    return(list(SNR=0, peak_cor=0))
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
            med_cor=median(peak_cor, na.rm=TRUE)) %>%
  mutate(log_med_cor=log10(1-med_cor))
write.csv(peakshape_mets, paste0(output_folder, "peakshape_mets.csv"), row.names = FALSE)


scan_time_diff <- diff(sort(unique(msdata$EIC2[filename==filename[1]]$rt)))
hist(scan_time_diff)
scan_time <- mean(scan_time_diff)

n_scans <- eic_dt[, .N, .(filename, feature)]
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
write.csv(med_missed_scans, paste0(output_folder, "med_missed_scans.csv"), row.names = FALSE)


# Currently broken because NAs are included during RT correction
# rt_dt <- rtime(msnexp_filled) %>%
#   as.data.frame() %>%
#   setNames("rt") %>%
#   rownames_to_column("filenum") %>%
#   mutate(filenum=as.numeric(str_extract(filenum, "(?<=F)\\d+"))) %>%
#   mutate(filename=basename(fileNames(msnexp_filled))[filenum]) %>%
#   mutate(rt=round(rt/60, digits = 7)) %>%
#   select(filename, rt) %>%
#   as.data.table()
# file_feat_n_missed <- eic_dt %>%
#   mutate(rt=round(rt, digits = 7)) %>%
#   group_by(filename, feature) %>%
#   group_split() %>%
#   pbapply::pblapply(function(eic){
#     fname_i <- eic$filename[1]
#     min_rt=min(eic$rt)
#     max_rt=max(eic$rt)
#     psb_rts <- rt_dt[rt%between%c(min_rt, max_rt)][filename==fname_i]
#     merged_dt <- merge(unique(eic), psb_rts, all.y=TRUE)
#     summarise(merged_dt, 
#               filename=unique(filename), 
#               feature=unique(feature),
#               n_missed=sum(is.na(merged_dt$int)),
#               n_scans=n()) %>%
#       na.omit()
#   }) %>%
#   bind_rows()
# med_missed_scans_2 <- file_feat_n_missed %>%
#   group_by(feature) %>%
#   summarise(med_missed_scans_2=median(n_missed, na.rm=TRUE))


# Calculate presence/absence of isotope info ----
msdata <- readRDS(paste0(output_folder, "msdata.rds"))
msdata$EIC2 <- msdata$EIC2[order(filename, rt, int)]
msdata_isoc <- readRDS(paste0(output_folder, "msdata_isoc.rds"))
msdata_isoc$EIC2 <- msdata_isoc$EIC2[order(filename, rt, int)]
filesplit_msdata_EIC <- split(msdata$EIC2, msdata$EIC2$filename)
filesplit_msdata_isoc_EIC <- split(msdata_isoc$EIC2, msdata_isoc$EIC2$filename)
peak_isodata <- peak_bounds %>%
  mutate(filename=basename(filename)) %>%
  mutate(rtmin=rtmin/60, rtmax=rtmax/60) %>%
  # filter(feature=="FT0356") %>%
  # filter(filename=="190715_Std_4uMStdsMix2InH2O_2.mzML") %>%
  pbapply::pbmapply(FUN = function(df, feature_i, filename_i, mzmin_i, mzmax_i, 
                                    rtmin_i, rtmax_i){
    # print(filename_i)
    init_eic <- filesplit_msdata_EIC[[filename_i]][
      mz%between%pmppm(sum(mzmin_i, mzmax_i)/2)][
        rt%between%c(rtmin_i, rtmax_i), c("rt", "int")][
          , rt:=round(rt, 5)]
    addiso_eic <- filesplit_msdata_isoc_EIC[[filename_i]][
      mz%between%pmppm(sum(mzmin_i, mzmax_i)/2+1.003355)][
        rt%between%c(rtmin_i, rtmax_i), c("rt", "int")][
          , rt:=round(rt, 5)]
    if(nrow(addiso_eic)<5){
      iso_cor <- 0
      if(nrow(init_eic)>0){
        init_area <- trapz(init_eic$rt, init_eic$int)
      } else {
        init_area <- 0
      }
      iso_area <- 0
    } else {
      eic <- merge(init_eic, addiso_eic, by="rt", suffixes=c("_init", "_iso"))
      if(nrow(eic)==0){
        iso_cor <- 0
        init_area <- 0
        iso_area <- 0
      } else {
        # par(mfrow=c(2,1), mar=c(2.1, 2.1, 0.1, 0.1))
        # plot(eic$rt, eic$int_init, xlab="", ylab="")
        # plot(eic$rt, eic$int_iso, xlab="", ylab="")
        iso_cor <- cor(eic$int_init, eic$int_iso, use="pairwise")
        init_area <- trapz(eic$rt, eic$int_init)
        iso_area <- trapz(eic$rt, eic$int_iso)
      }
    }
    list(feature=feature_i, filename=filename_i, iso_cor=iso_cor, 
      init_area=init_area, iso_area=iso_area)
  }, .$feature, .$filename, .$mzmin, .$mzmax, .$rtmin, .$rtmax, SIMPLIFY = FALSE) %>%
  bind_rows()
write.csv(peak_isodata, file = paste0(output_folder, "peak_isodata.csv"), row.names = FALSE)
peak_isodata <- read_csv(paste0(output_folder, "peak_isodata.csv"))
feat_isodata <- peak_isodata %>%
  group_by(feature) %>%
  summarise(shape_cor=log10(1-median(iso_cor)), area_cor=log10(1-cor(init_area, iso_area))) %>%
  mutate(area_cor=ifelse(is.na(area_cor), 0, area_cor))
write.csv(feat_isodata, paste0(output_folder, "feat_isodata.csv"), row.names = FALSE)

# Calculate DOE metrics ----
depth_diffs <- peak_data %>%
  select(feat_id, filename, into) %>%
  left_join(file_data) %>%
  filter(samp_type=="Smp") %>%
  nest(data=-feat_id) %>%
  mutate(t_pval=map_dbl(data, function(x){
    depthtab <- table(x$depth)
    if(any(depthtab<2))return(1)
    if(length(depthtab)<2)return(0.001)
    broom::tidy(aov(x$into~x$depth))$p.value[1]
  })) %>%
  mutate(t_pval=log10(t_pval)) %>%
  select(feat_id, t_pval)
write.csv(depth_diffs, paste0(output_folder, "depth_diffs.csv"), row.names = FALSE)

blank_diffs <- peak_data %>%
  select(feat_id, filename, into) %>%
  left_join(file_data) %>%
  filter(samp_type%in%c("Smp", "Blk")) %>%
  select(feat_id, samp_type, into) %>%
  group_by(feat_id, samp_type) %>%
  summarise(mean_area=mean(into)) %>%
  pivot_wider(names_from = samp_type, values_from = mean_area) %>%
  mutate(smp_to_blk=Smp/Blk) %>%
  mutate(smp_to_blk=ifelse(is.na(smp_to_blk), 10000, smp_to_blk)) %>%
  select(feat_id, smp_to_blk) %>%
  mutate(smp_to_blk=log10(smp_to_blk))
write.csv(blank_diffs, paste0(output_folder, "blank_diffs.csv"), row.names = FALSE)

stan_diffs <- peak_data %>%
  select(feat_id, filename, into) %>%
  left_join(file_data) %>%
  filter(samp_type%in%c("Smp", "Std")) %>%
  select(feat_id, samp_type, into) %>%
  group_by(feat_id, samp_type) %>%
  summarise(mean_area=mean(into)) %>%
  pivot_wider(names_from = samp_type, values_from = mean_area) %>%
  mutate(smp_to_std=Smp/Std) %>%
  mutate(smp_to_std=ifelse(is.na(smp_to_std), 10, smp_to_std)) %>%
  select(feat_id, smp_to_std) %>%
  mutate(smp_to_std=log10(smp_to_std))
write.csv(stan_diffs, paste0(output_folder, "stan_diffs.csv"), row.names = FALSE)


# Join feature classes and write out ----
simple_feats <- read_csv(paste0(output_folder, "simple_feats.csv"))
peakshape_mets <- read_csv(paste0(output_folder, "peakshape_mets.csv"))
med_missed_scans <- read_csv(paste0(output_folder, "med_missed_scans.csv"))
depth_diffs <- read_csv(paste0(output_folder, "depth_diffs.csv"))
blank_diffs <- read_csv(paste0(output_folder, "blank_diffs.csv"))
stan_diffs <- read_csv(paste0(output_folder, "stan_diffs.csv"))
feat_isodata <- read_csv(paste0(output_folder, "feat_isodata.csv"))
classified_feats <- read_csv(paste0(output_folder, "classified_feats.csv")) %>%
  select(feature, feat_class)


features_extracted <- simple_feats %>%
  left_join(peakshape_mets, by=c(feat_id="feature")) %>%
  left_join(med_missed_scans, by=c(feat_id="feature")) %>%
  left_join(depth_diffs) %>%
  left_join(blank_diffs) %>%
  left_join(stan_diffs) %>%
  left_join(feat_isodata, by=c(feat_id="feature")) %>% 
  left_join(classified_feats, by=c(feat_id="feature")) %>%
  drop_na()
write.csv(features_extracted, paste0(output_folder, "features_extracted.csv"), 
          row.names = FALSE)


# Visualization ----

metric_plots <- features_extracted %>%
  names() %>%
  setdiff(c("feat_id", "blank_found", "feat_class")) %>%
  lapply(function(col_name){
    features_extracted %>%
      mutate(col_to_plot=cut(get(col_name), pretty(get(col_name), n = 10))) %>%
      ggplot() +
      geom_bar(aes(x=col_to_plot, fill=feat_class)) +
      labs(x=col_name) +
      theme(legend.position = "none")
  })

fullplot <- do.call(gridExtra::arrangeGrob, c(metric_plots, ncol=1))
ggsave("all_metrics.pdf", plot = fullplot, device = "pdf", path = "figures",
       width=5, height = 50, units = "in", limitsize = FALSE)




library(plotly)
plot_ly(features_extracted, x=~mean_pw, y=~log10(sn), z=~n_found, color=~feat_class,
        type = "scatter3d", mode="markers")
plot_ly(features_extracted, x=~n_found, y=~samps_found, z=~blank_found, color=~feat_class,
        type = "scatter3d", mode="markers")
features_extracted %>%
  plot_ly(x=~med_SNR, y=~med_cor, z=~log_mean_height, 
          color=~feat_class, text=~feat_id,
          type = "scatter3d", mode="markers")
features_extracted %>%
  plot_ly(x=~log10(1-shape_cor), y=~med_SNR, z=~log10(1-med_cor), color=~feat_class,
          type = "scatter3d", mode="markers")


library(GGally)
gp <- features_extracted %>%
  select(-feat_id) %>%
  ggpairs(aes(color=feat_class), lower=list(
    combo=wrap("facethist",  bins=30))
    )
ggsave(plot = gp, filename = "pairsplot.pdf", device = "pdf", width = 25, height = 25, 
       units = "in", limitsize = FALSE, path = "figures")

