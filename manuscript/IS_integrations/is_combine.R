library(dplyr)

# Check that SkylineRunner.exe is installed
if(system("manuscript/IS_integrations/SkylineRunner.exe", show.output.on.console = FALSE)){
  # Website = https://skyline.ms/wiki/home/software/Skyline/page.view?name=SkylineInstall_64_21-1&submit=false
  stop("Please grab a copy of SkylineRunner from the website")
}

sky_neg <- paste0(
  "manuscript/IS_integrations/SkylineRunner.exe",
  " --in=manuscript/IS_integrations/HILIC_QE_NEGATIVE_Mesotransect.sky",
  ' --report-name="Ingalls_Lab_QE_Transition Results"', 
  ' --report-file="manuscript/IS_integrations/negative_mode.csv"'
)
sky_pos <- paste0(
  "manuscript/IS_integrations/SkylineRunner.exe",
  " --in=manuscript/IS_integrations/HILIC_QE_POSITIVE_Mesotransect.sky",
  ' --report-name="Ingalls_Lab_QE_Transition Results"', 
  ' --report-file="manuscript/IS_integrations/positive_mode.csv"'
)
system(sky_neg)
system(sky_pos)

IS_pos <- read.csv("manuscript/IS_integrations/positive_mode.csv")
IS_neg <- read.csv("manuscript/IS_integrations/negative_mode.csv")
IS_all <- bind_rows(IS_pos, IS_neg)

# Fix sample 180821_Smp_MS6C315m_C which didn't get extraction standards added
# Current fix is to replace with median Area of all other sample files
ext_stans_list <- c(
  "Isethionic Acid, 13C2",
  "L-Cysteic Acid, D3",
  "Sulfoacetic Acid, 13C2",
  "Sulfolactic Acid, 13C3",
  "Taurine, D4"
)
update_df <- IS_all %>%
  filter(grepl("Smp", Replicate.Name)) %>%
  filter(Precursor.Ion.Name%in%ext_stans_list) %>%
  mutate(Area=as.numeric(Area)) %>%
  group_by(Precursor.Ion.Name) %>%
  summarise(Area=median(Area)) %>%
  mutate(Replicate.Name="180821_Smp_MS6C315m_C")
IS_all <- IS_all %>%
  rows_update(update_df, by = c("Precursor.Ion.Name", "Replicate.Name")) %>%
  select(IS_name=Precursor.Ion.Name, filename=Replicate.Name, area=Area)
write.csv(IS_all, "manuscript/IS_integrations/all_IS.csv", row.names = FALSE)
