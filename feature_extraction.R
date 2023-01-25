
# Setup ----
library(xcms)
library(tidyverse)
library(RaMS)
options(pillar.sigfig=7)

file_data <- read_csv("made_data/file_data.csv") %>%
  mutate(filename=basename(filename))

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

features_extracted <- peak_data %>%
  group_by(feat_id) %>%
  summarise(mean_mz=unique(feat_mzmed), sd_mz=sd(mz), range_mz=max(mz)-min(mz),
            mean_rt=unique(feat_rtmed), sd_rt=sd(rt), range_rt=max(rt)-min(rt),
            mean_pw=mean(rtmax-rtmin), sd_pw=sd(rtmax-rtmin),
            # mean_area=mean(into), med_area=median(into), IQR(into),
            log_mean_height=log10(mean(maxo)), log_sd_height=log10(sd(maxo)),
            across(c(sn, f, scale, scpos, lmin, lmax), mean, na.rm=TRUE),
            n_found=41-sum(is.na(intb)), 
            samps_found=1-sum(is.na(intb) & str_detect(filename, "Smp"))/24,
            stans_found=1-sum(is.na(intb) & str_detect(filename, "Std"))/10,
            blank_found=any(!is.na(intb) & str_detect(filename, "Blk"))) %>%
  left_join(classified_feats, by=c(feat_id="feature"))
write.csv(features_extracted, "made_data/features_extracted.csv", row.names = FALSE)

with(features_extracted, plot(mean_mz, sd_mz))

library(plotly)
plot_ly(features_extracted, x=~mean_pw, y=~log10(sn), z=~n_found, color=~feat_class,
        type = "scatter3d", mode="markers")
plot_ly(features_extracted, x=~sd_rt, y=~log10(mean_area), z=~sd_mz, color=~feat_class,
        type = "scatter3d", mode="markers")


