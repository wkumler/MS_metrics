
library(tidyverse)
filled_file_metadata <- read_csv("MW_upload/filled_file_metadata.csv")
untarg_peaks <- read_csv("MW_upload/untarg_peaks.csv")
untarg_feats <- read_csv("MW_upload/untarg_feats.csv")
targ_peaks <- read_csv("MW_upload/targ_peaks.csv")

## Falkor ----
### Metadata
filled_file_metadata %>%
  filter(cruise=="FK") %>%
  filter(polarity=="pos") %>%
  mutate(Sample_ID=filename) %>%
  mutate(Raw_file_name=filename) %>%
  select(Raw_file_name, Sample_ID, station, cast, depth, abs_depth,
         time, lat, lon, temp, sal, sla_corr, PC_um) %>%
  write.table("MW_upload/falkor_metadata.txt", sep = "\t", row.names = FALSE)

### Untargeted data
untarg_peaks %>%
  filter(cruise=="FK") %>%
  select(-cruise, -polarity) %>%
  mutate(norm_area=ifelse(!is.finite(norm_area), NA, norm_area)) %>%
  pivot_wider(names_from = filename, values_from = norm_area) %>%
  left_join(all_feats %>% select(feature, mzmed, rtmed)) %>%
  mutate(mz_rt=paste0(round(mzmed, 5), "_", round(rtmed, 2))) %>%
  select(mz_rt, starts_with("190715")) %>%
  write.table("MW_upload/falkor_untargeted.txt", sep = "\t", row.names = FALSE)

### Targeted data
metabs <- unique(targ_peaks$compound_name)
metabs_enc <- URLencode(paste(metabs, collapse = "\n"))
metab_url <- "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet.php"
metab_body <- paste0("metabolite_name=", metabs_enc, "&submit=Submit")
heads <- c(`Content-Type`= 'application/x-www-form-urlencoded')
stan_names <- httr::POST(metab_url, body = metab_body, config = httr::add_headers(heads)) %>%
  httr::content() %>%
  rvest::html_table() %>%
  `[[`(2)

targ_peaks %>%
  filter(str_detect(filename, "190715")) %>%
  mutate(filename=paste0(filename, ".mzML")) %>%
  mutate(nM=ifelse(!is.finite(nM), NA, nM)) %>%
  mutate(nM=signif(nM, 4)) %>%
  select(-polarity) %>%
  pivot_wider(names_from = filename, values_from = nM) %>%
  left_join(stan_names %>% select(compound_name=`Input name`, stan_name=`Standardized name`)) %>%
  mutate(compound_name=ifelse(stan_name=="-", compound_name, stan_name)) %>%
  select(-stan_name) %>%
  write.table("MW_upload/falkor_targeted.txt", sep = "\t", row.names = FALSE)

### Zip mzMLs

## MESOSCOPE ----

## CultureData ----

## Pttime ----
file_order <- c("1_extr-ctrl--extr-ctrl-_1", 
                "2_extr-ctrl--extr-ctrl-_2", "3_extr-ctrl--extr-ctrl-_3", "4_extr-ctrl--extr-ctrl-_4", 
                "10_exudate-0d_2", "11_exudate-0d_3", "12_exudate-0d_4", "9_exudate-0d_1", 
                "25_exudate-12d_1", "26_exudate-12d_2", "27_exudate-12d_3", "28_exudate-12d_4", 
                "13_exudate-3d_1", "14_exudate-3d_2", "15_exudate-3d_3", "16_exudate-3d_4", 
                "17_exudate-6d_1", "18_exudate-6d_2", "19_exudate-6d_3", "20_exudate-6d_4", 
                "21_exudate-9d_1", "22_exudate-9d_2", "23_exudate-9d_3", "24_exudate-9d_4", 
                "5_exudate-med-ctrl-_1", "6_exudate-med-ctrl-_2", "7_exudate-med-ctrl-_3", 
                "8_exudate-med-ctrl-_4", "33_pellet-0d_1", "34_pellet-0d_2", 
                "35_pellet-0d_3", "36_pellet-0d_4", "49_pellet-12d_1", "50_pellet-12d_2", 
                "51_pellet-12d_3", "52_pellet-12d_4", "37_pellet-3d_1", "38_pellet-3d_2", 
                "39_pellet-3d_3", "40_pellet-3d_4", "41_pellet-6d_1", "42_pellet-6d_2", 
                "43_pellet-6d_3", "44_pellet-6d_4", "45_pellet-9d_1", "46_pellet-9d_2", 
                "47_pellet-9d_3", "48_pellet-9d_4", "29_pellet-med-ctrl-_1", 
                "30_pellet-med-ctrl-_2", "31_pellet-med-ctrl-_3", "32_pellet-med-ctrl-_4"
)
pttime_quals <- read_csv("made_data_Pttime/classified_feats.csv")
library(xcms)
pttime_msnexp <- readRDS("made_data_Pttime/msnexp_filled.rds")
peak_data_long <- pttime_msnexp %>%
  chromPeaks() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(peakidx=row_number())
peak_data <- pttime_msnexp %>%
  featureDefinitions() %>%
  as.data.frame() %>%
  select(mzmed, rtmed, npeaks, peakidx) %>%
  rownames_to_column("id") %>%
  unnest_longer(peakidx) %>%
  rename_with(~paste0("feat_", .x), .cols = -peakidx) %>%
  dplyr::rename(feature="feat_id") %>%
  left_join(peak_data_long, by = join_by(peakidx)) %>%
  mutate(filename=basename(fileNames(pttime_msnexp))[sample])
peak_data %>%
  mutate(feat_rtmed=feat_rtmed/60) %>%
  mutate(filename=str_remove(filename, "20190429_JJ_VB_BioSFA_Pttime_1_QE144_Ag68377-924_USHXG01160_POS_MSMS-v2_")) %>%
  mutate(filename=str_remove(filename, "(?<=_\\d).+")) %>%
  mutate(mz_rt=paste0(round(feat_mzmed, 5), "_", round(feat_rtmed, 2))) %>%
  left_join(pttime_quals %>% select(feature, feat_class)) %>%
  filter(feat_class=="Good") %>%
  select(mz_rt, filename, into) %>%
  distinct() %>%
  arrange(desc(into)) %>%
  group_by(mz_rt, filename) %>%
  slice(1) %>%
  ungroup() %>%
  pivot_wider(names_from=filename, values_from=into) %>%
  select(mz_rt, all_of(file_order)) %>%
  write.table("MW_upload/pttime_untargeted.txt", sep = "\t", row.names = FALSE)
