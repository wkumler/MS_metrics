

# Setup ----
library(tidyverse)
library(xcms)
library(data.table)
file_data <- read_csv("made_data_MS3000/file_data.csv") %>%
  mutate(filename=basename(filename))
msnexp_filled <- readRDS("made_data_MS3000/msnexp_filled.rds")

# Reorganize peak data from XCMS into clean format ----
peak_data_long <- msnexp_filled %>%
  chromPeaks() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(peakidx=row_number()) %>%
  select(sample, peakidx, into)
peak_data <- msnexp_filled %>%
  featureDefinitions() %>%
  as.data.frame() %>%
  filter(rtmed%between%c(30, 1200)) %>%
  # filter(mzmed%between%pmppm(116.071154, 5))
  select(peakidx) %>%
  rownames_to_column("feature") %>%
  unnest_longer(peakidx) %>%
  left_join(peak_data_long, by = join_by(peakidx)) %>%
  mutate(filename=basename(fileNames(msnexp_filled))[sample]) %>%
  complete(feature, filename, fill = list(into=0)) %>%
  left_join(file_data, by = join_by(filename)) %>%
  filter(str_detect(filename, "180821"))

# Perform B-MIS ----
source("manuscript/IS_integrations/is_combine.R")
IS_areas <- read_csv("manuscript/IS_integrations/all_IS.csv") %>%
  mutate(filename=paste0(filename, ".mzML")) %>%
  rename(area="IS_area") %>%
  filter(!str_detect(filename, "QC")) %>%
  pivot_wider(names_from=IS_name, values_from = "IS_area") %>%
  mutate("Inj_vol"=ifelse(str_detect(filename, "Poo.*Half"), 0.5, 1)) %>%
  pivot_longer(-filename, names_to = "IS_name", values_to = "IS_area")
pooled_IS_areas <- IS_areas %>%
  filter(str_detect(filename, "_Poo_"))

# Check using Proline
peak_data %>%
  select(feature, filename, into, depth) %>%
  filter(str_detect(filename, "_Poo_")) %>%
  filter(feature=="FT0532") %>%
  full_join(pooled_IS_areas, multiple = "all", by = join_by(filename)) %>%
  group_by(IS_name) %>%
  mutate(plot_area=into/IS_area*mean(IS_area)) %>%
  ggplot() +
  geom_col(aes(x=filename, y=plot_area)) +
  facet_wrap(~IS_name, ncol = 2)

# init_cvs <- peak_data %>%
#   filter(str_detect(filename, "_Poo_")) %>%
#   group_by(feature) %>%
#   summarize(init_cv=sd(into)/mean(into))
  
BMIS <- peak_data %>%
  filter(str_detect(filename, "_Poo_")) %>%
  full_join(pooled_IS_areas, multiple = "all", by = join_by(filename)) %>%
  group_by(feature, IS_name) %>%
  mutate(plot_area=into/IS_area*mean(IS_area)) %>%
  summarise(cv=sd(plot_area)/mean(plot_area)) %>%
  # left_join(init_cvs) %>%
  arrange(cv) %>%
  slice(1) %>%
  select(feature, IS_name)

BMISed_data <- peak_data %>%
  filter(samp_type=="Smp") %>%
  select(feature, filename, into) %>%
  left_join(BMIS, by = join_by(feature)) %>%
  left_join(IS_areas, by = join_by(filename, IS_name)) %>%
  group_by(feature) %>%
  mutate(norm_area=into/IS_area*mean(IS_area)) %>%
  select(feature, filename, norm_area)

norm_area_data <- BMISed_data %>%
  arrange(norm_area) %>%
  group_by(feature, filename) %>%
  slice(1) %>%
  ungroup()

# Write out ----
if(any(is.na(norm_area_data$norm_area)))stop("Found NAs in norm_area!")
if(!all(table(norm_area_data$feature)==99))stop("Found not 99 samps in plotready!")
write_csv(norm_area_data, "manuscript/MS3000_norm_areas.csv")
