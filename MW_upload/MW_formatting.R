
library(tidyverse)
filled_file_metadata <- read_csv("MW_upload/filled_file_metadata.csv")
untarg_peaks <- read_csv("MW_upload/untarg_peaks.csv")
untarg_feats <- read_csv("MW_upload/untarg_feats.csv")
targ_peaks <- read_csv("MW_upload/targ_peaks.csv")
clean_stans <- read_csv("MW_upload/clean_stans.csv")

metabs <- unique(targ_peaks$compound_name)
metabs_enc <- URLencode(paste(metabs, collapse = "\n"))
metab_url <- "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet.php"
metab_body <- paste0("metabolite_name=", metabs_enc, "&submit=Submit")
heads <- c(`Content-Type`= 'application/x-www-form-urlencoded')
stan_names <- httr::POST(metab_url, body = metab_body, config = httr::add_headers(heads)) %>%
  httr::content() %>%
  rvest::html_table() %>%
  `[[`(2)


## Falkor ----
### Metadata
filled_file_metadata %>%
  filter(cruise=="FK") %>%
  filter(polarity=="pos") %>%
  mutate(Sample_ID=filename) %>%
  mutate(Raw_file_name=filename) %>%
  mutate(samp_type=case_when(
    samp_type=="Smp"~"Sample",
    samp_type=="Std"~"Standard mix",
    samp_type=="Poo"~"QC Pooled",
    samp_type=="Blk"~"Blank"
  )) %>%
  mutate(across(lat:PC_um, signif)) %>%
  select(Raw_file_name, Sample_ID, Sample_type=samp_type, 
         Station=station, Depth=depth, Absolute_depth=abs_depth,
         Time_sampled=time, Latitude=lat, Longitude=lon, 
         Temperature=temp, Salinity=sal, Sea_level_anomaly=sla_corr, 
         Particulate_carbon=PC_um) %>%
  # distinct(Station, Time_sampled)
  write.table("MW_upload/falkor_metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

### Untargeted data
untarg_peaks %>%
  filter(cruise=="FK") %>%
  filter(polarity=="pos") %>%
  select(-cruise, -polarity) %>%
  mutate(norm_area=ifelse(!is.finite(norm_area), NA, norm_area)) %>%
  pivot_wider(names_from = filename, values_from = norm_area) %>%
  left_join(untarg_feats %>% select(feature, mzmed, rtmed)) %>%
  mutate(mz_rt=paste0(round(mzmed, 5), "_", round(rtmed, 2))) %>%
  select(mz_rt, starts_with("190715")) %>%
  write.table("MW_upload/falkor_pos_untargeted.txt", sep = "\t", row.names = FALSE, quote = FALSE)
untarg_peaks %>%
  filter(cruise=="FK") %>%
  filter(polarity=="neg") %>%
  select(-cruise, -polarity) %>%
  mutate(norm_area=ifelse(!is.finite(norm_area), NA, norm_area)) %>%
  pivot_wider(names_from = filename, values_from = norm_area) %>%
  left_join(untarg_feats %>% select(feature, mzmed, rtmed)) %>%
  mutate(mz_rt=paste0(round(mzmed, 5), "_", round(rtmed, 2))) %>%
  select(mz_rt, starts_with("190715")) %>%
  write.table("MW_upload/falkor_neg_untargeted.txt", sep = "\t", row.names = FALSE, quote = FALSE)


### Targeted data
targ_peaks %>%
  filter(str_detect(filename, "190715")) %>%
  mutate(filename=paste0(filename, ".mzML")) %>%
  mutate(nM=ifelse(!is.finite(nM), NA, nM)) %>%
  mutate(nM=signif(nM, 4)) %>%
  filter(polarity=="pos") %>%
  select(-polarity) %>%
  pivot_wider(names_from = filename, values_from = nM) %>%
  left_join(stan_names %>% select(compound_name=`Input name`, stan_name=`Standardized name`)) %>%
  mutate(compound_name=ifelse(stan_name=="-", compound_name, stan_name)) %>%
  select(-stan_name) %>%
  write.table("MW_upload/falkor_pos_targeted.txt", sep = "\t", row.names = FALSE, quote = FALSE)
targ_peaks %>%
  filter(polarity=="pos") %>%
  distinct(compound_name) %>%
  left_join(clean_stans %>% filter(polarity=="pos")) %>%
  left_join(stan_names %>% select(compound_name=`Input name`, stan_name=`Standardized name`)) %>%
  mutate(compound_name=ifelse(stan_name=="-", compound_name, stan_name)) %>%
  select(compound_name, formula, stan_mz, stan_rt_expected=stan_rt, which_standard_mix=mix, 
         concentration_added=conc_um, kegg_id) %>%
  write.table("MW_upload/falkor_pos_targeted_metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

targ_peaks %>%
  filter(str_detect(filename, "190715")) %>%
  mutate(filename=paste0(filename, ".mzML")) %>%
  mutate(nM=ifelse(!is.finite(nM), NA, nM)) %>%
  mutate(nM=signif(nM, 4)) %>%
  filter(polarity=="neg") %>%
  select(-polarity) %>%
  pivot_wider(names_from = filename, values_from = nM) %>%
  left_join(stan_names %>% select(compound_name=`Input name`, stan_name=`Standardized name`)) %>%
  mutate(compound_name=ifelse(stan_name=="-", compound_name, stan_name)) %>%
  select(-stan_name) %>%
  write.table("MW_upload/falkor_neg_targeted.txt", sep = "\t", row.names = FALSE, quote = FALSE)
targ_peaks %>%
  filter(polarity=="neg") %>%
  distinct(compound_name) %>%
  left_join(clean_stans %>% filter(polarity=="neg")) %>%
  left_join(stan_names %>% select(compound_name=`Input name`, stan_name=`Standardized name`)) %>%
  mutate(compound_name=ifelse(stan_name=="-", compound_name, stan_name)) %>%
  select(compound_name, formula, stan_mz, stan_rt_expected=stan_rt, which_standard_mix=mix, 
         concentration_added=conc_um, kegg_id) %>%
  write.table("MW_upload/falkor_neg_targeted_metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)


### Zip mzMLs
# Done with just a "new compressed folder" and then copy/pasting the mzMLs to it

## MESOSCOPE ----
### Metadata
filled_file_metadata %>%
  filter(cruise=="MS") %>%
  filter(polarity=="pos") %>%
  mutate(Sample_ID=filename) %>%
  mutate(Raw_file_name=filename) %>%
  mutate(samp_type=case_when(
    samp_type=="Smp"~"Sample",
    samp_type=="Std"~"Standard mix",
    samp_type=="Poo"~"QC Pooled",
    samp_type=="Blk"~"Blank"
  )) %>%
  mutate(across(lat:PC_um, signif)) %>%
  select(Raw_file_name, Sample_ID, Sample_type=samp_type, 
         Station=station, Depth=depth, Absolute_depth=abs_depth,
         Time_sampled=time, Latitude=lat, Longitude=lon, 
         Temperature=temp, Salinity=sal, Sea_level_anomaly=sla_corr, 
         Particulate_carbon=PC_um) %>%
  distinct(Raw_file_name, Sea_level_anomaly) %>% print(n=Inf)
  write.table("MW_upload/mesoscope_metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

### Untargeted data
untarg_peaks %>%
  filter(cruise=="MS") %>%
  filter(polarity=="pos") %>%
  select(-cruise, -polarity) %>%
  mutate(norm_area=ifelse(!is.finite(norm_area), NA, norm_area)) %>%
  pivot_wider(names_from = filename, values_from = norm_area) %>%
  left_join(untarg_feats %>% select(feature, mzmed, rtmed)) %>%
  mutate(mz_rt=paste0(round(mzmed, 5), "_", round(rtmed, 2))) %>%
  select(mz_rt, everything(), -c(feature, mzmed, rtmed)) %>% 
  write.table("MW_upload/mesoscope_pos_untargeted.txt", sep = "\t", row.names = FALSE, quote = FALSE)
untarg_peaks %>%
  filter(cruise=="MS") %>%
  filter(polarity=="neg") %>%
  select(-cruise, -polarity) %>%
  mutate(norm_area=ifelse(!is.finite(norm_area), NA, norm_area)) %>%
  pivot_wider(names_from = filename, values_from = norm_area) %>%
  left_join(untarg_feats %>% select(feature, mzmed, rtmed)) %>%
  mutate(mz_rt=paste0(round(mzmed, 5), "_", round(rtmed, 2))) %>%
  select(mz_rt, everything(), -c(feature, mzmed, rtmed)) %>% 
  write.table("MW_upload/mesoscope_neg_untargeted.txt", sep = "\t", row.names = FALSE, quote = FALSE)

### Targeted data
ms_filenames <- untarg_peaks %>%
  filter(cruise=="MS") %>%
  filter(polarity=="pos") %>%
  distinct(filename) %>%
  pull(filename)
targ_peaks %>%
  mutate(filename=paste0(filename, ".mzML")) %>%
  filter(filename%in%ms_filenames) %>%
  mutate(nM=ifelse(!is.finite(nM), NA, nM)) %>%
  mutate(nM=signif(nM, 4)) %>%
  filter(polarity=="pos") %>%
  select(-polarity) %>%
  pivot_wider(names_from = filename, values_from = nM) %>%
  left_join(stan_names %>% select(compound_name=`Input name`, stan_name=`Standardized name`)) %>%
  mutate(compound_name=ifelse(stan_name=="-", compound_name, stan_name)) %>%
  select(-stan_name) %>% 
  write.table("MW_upload/mesoscope_pos_targeted.txt", sep = "\t", row.names = FALSE, quote = FALSE)
targ_peaks %>%
  filter(polarity=="pos") %>%
  distinct(compound_name) %>%
  left_join(clean_stans %>% filter(polarity=="pos")) %>%
  left_join(stan_names %>% select(compound_name=`Input name`, stan_name=`Standardized name`)) %>%
  mutate(compound_name=ifelse(stan_name=="-", compound_name, stan_name)) %>%
  select(compound_name, formula, stan_mz, stan_rt_expected=stan_rt, which_standard_mix=mix, 
         concentration_added=conc_um, kegg_id) %>%
  write.table("MW_upload/mesoscope_pos_targeted_metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

targ_peaks %>%
  mutate(filename=paste0(filename, ".mzML")) %>%
  filter(filename%in%ms_filenames) %>%
  mutate(nM=ifelse(!is.finite(nM), NA, nM)) %>%
  mutate(nM=signif(nM, 4)) %>%
  filter(polarity=="neg") %>%
  select(-polarity) %>%
  pivot_wider(names_from = filename, values_from = nM) %>%
  left_join(stan_names %>% select(compound_name=`Input name`, stan_name=`Standardized name`)) %>%
  mutate(compound_name=ifelse(stan_name=="-", compound_name, stan_name)) %>%
  select(-stan_name) %>%
  write.table("MW_upload/mesoscope_neg_targeted.txt", sep = "\t", row.names = FALSE, quote = FALSE)
targ_peaks %>%
  filter(polarity=="neg") %>%
  distinct(compound_name) %>%
  left_join(clean_stans %>% filter(polarity=="neg")) %>%
  left_join(stan_names %>% select(compound_name=`Input name`, stan_name=`Standardized name`)) %>%
  mutate(compound_name=ifelse(stan_name=="-", compound_name, stan_name)) %>%
  select(compound_name, formula, stan_mz, stan_rt_expected=stan_rt, which_standard_mix=mix, 
         concentration_added=conc_um, kegg_id) %>%
  write.table("MW_upload/mesoscope_neg_targeted_metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)


### Zip mzMLs
# Done with just a "new compressed folder" and then copy/pasting the mzMLs to it





## CultureData ----
### Metadata
culturedata_fullnames <- read_csv("made_data_CultureData/Culture Status - Sheet1.csv") %>%
  set_names("org_name", "organism", "org_class") %>%
  mutate(organism=case_when(
    organism=="K-0643"~"Ca",
    organism=="K-0644"~"Cr",
    organism=="K-0272"~"Cs",
    organism=="AS9601"~"As9601",
    organism=="NAT"~"Nat",
    organism=="1161"~"116-1",
    organism=="P5"~"P5-5",
    TRUE~organism
  )) %>%
  drop_na()
culturedata_files <- read_csv("made_data_CultureData/file_data.csv") %>%
  select(filename, samp_type, organism, run_group=org_class) %>%
  mutate(filename=basename(filename)) %>%
  mutate(Sample_ID=filename) %>%
  mutate(Raw_file_name=filename) %>%
  mutate(samp_type=case_when(
    samp_type=="Smp"~"Sample",
    samp_type=="Std"~"Standard mix",
    samp_type=="Poo"~"QC Pooled",
    samp_type=="Blk"~"Blank"
  )) %>%
  left_join(culturedata_fullnames) %>%
  select(Raw_file_name, Sample_ID, Sample_type=samp_type, Sample_type=samp_type,
         Species_strain=org_name, Broad_taxon=org_class, MS_run=run_group) %>%
  write.table("MW_upload/culturedata_metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

### Untargeted
culturedata_quals <- read_csv("made_data_CultureData/classified_feats.csv")
library(xcms)
culturedata_msnexp <- readRDS("made_data_CultureData/msnexp_filled.rds")
peak_data_long <- culturedata_msnexp %>%
  chromPeaks() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(peakidx=row_number())
peak_data <- culturedata_msnexp %>%
  featureDefinitions() %>%
  as.data.frame() %>%
  select(mzmed, rtmed, npeaks, peakidx) %>%
  rownames_to_column("id") %>%
  unnest_longer(peakidx) %>%
  rename_with(~paste0("feat_", .x), .cols = -peakidx) %>%
  dplyr::rename(feature="feat_id") %>%
  left_join(peak_data_long, by = join_by(peakidx)) %>%
  mutate(filename=basename(fileNames(culturedata_msnexp))[sample])
peak_data %>%
  mutate(feat_rtmed=feat_rtmed/60) %>%
  mutate(mz_rt=paste0(round(feat_mzmed, 5), "_", round(feat_rtmed, 2))) %>%
  left_join(culturedata_quals %>% select(feature, feat_class)) %>%
  filter(feat_class=="Good") %>%
  select(mz_rt, filename, into) %>%
  distinct() %>%
  arrange(desc(into)) %>%
  group_by(mz_rt, filename) %>%
  slice(1) %>%
  ungroup() %>%
  pivot_wider(names_from=filename, values_from=into) %>%
  write.table("MW_upload/culturedata_untargeted.txt", sep = "\t", row.names = FALSE, quote = FALSE)



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
