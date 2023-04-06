# This script takes in the mzML files supplied in the mzMLs/ directory
# and the output from the peakpicking_and_prep.R script (csv files with
# peak boundaries, retention time correction, and file metadata)
# and allows the user to interactively assess the quality of various peaks
# and use keyboard hotkeys (up/down/left/right) to classify them. Classified
# peaks are then exported into made_data/classified_feats.csv

# This is the second script in the pipeline

# Setup ----
library(tidyverse)
library(RaMS)
options(pillar.sigfig=7)

# dataset_version <- "FT2040"
# dataset_version <- "MS3000"
dataset_version <- "CultureData"
# dataset_version <- "Pttime"
output_folder <- paste0("made_data_", dataset_version, "/")

file_data <- read_csv(paste0(output_folder, "file_data.csv")) %>%
  mutate(filename=basename(filename))
peak_bounds <- read_csv(paste0(output_folder, "peak_bounds.csv"))
rt_corrections <- read_csv(paste0(output_folder, "rt_corrections.csv"))



# Extract feature data from raw mzML files ----
if(dataset_version%in%c("FT2040", "MS3000")){
  feature_df <- peak_bounds %>%
    group_by(feature) %>%
    summarize(min_mz=min(mzmin), max_mz=max(mzmax), 
              min_rt=min(rtmin), max_rt=max(rtmax)) %>%
    mutate(mean_mz=(min_mz+max_mz)/2)
} else {
  FT2040_features <- read_csv("made_data_FT2040/features_extracted.csv")
  MS3000_features <- read_csv("made_data_MS3000/features_extracted.csv")
  both_min_model <- rbind(FT2040_features, MS3000_features) %>%
    select(feat_class, med_cor, med_SNR) %>%
    filter(feat_class%in%c("Good", "Bad")) %>%
    mutate(feat_class=feat_class=="Good") %>%
    glm(formula=feat_class~., family = binomial)
  good_features <- read_csv(paste0(output_folder, "features_extracted.csv")) %>%
    mutate(pred_prob=predict(object=both_min_model,newdata = ., type = "response")) %>%
    filter(pred_prob>0.9)
  feature_df <- read_csv(paste0(output_folder, "peak_bounds.csv")) %>%
    group_by(feature) %>%
    summarize(min_mz=min(mzmin), max_mz=max(mzmax), 
              min_rt=min(rtmin), max_rt=max(rtmax)) %>%
    mutate(mean_mz=(min_mz+max_mz)/2) %>%
    filter(feature%in%good_features$feature)
  
  msdata <- file_data$filename %>%
    paste0(output_folder, "mzMLs/", .) %>% 
    grabMSdata(verbosity = 1, grab_what = "EIC",
               mz=feature_df$mean_mz, ppm = 10)
  msdata$EIC2 <- msdata$EIC %>%
    mutate(rt=round(rt, 10)) %>%
    left_join(rt_corrections, by=c("filename", "rt")) %>%
    select(rt=new_rt, mz, int, filename)
  saveRDS(msdata, file = paste0(output_folder, "msdata_only_good.rds"))
}
write.csv(feature_df, paste0(output_folder, "feature_df.csv"), row.names = FALSE)

feature_df <- read_csv(paste0(output_folder, "feature_df.csv"))


# Read in existing classification doc and find last classified feature ----
if(file.exists(paste0(output_folder, "classified_feats.csv"))){
  data_classified <- read_csv(paste0(output_folder, "classified_feats.csv"))
  feat_id <- tail(data_classified$feature, 1)
} else {
  data.frame(
    feature=character(),
    min_mz=numeric(),
    max_mz=numeric(),
    min_rt=numeric(),
    max_rt=numeric(),
    mean_mz=numeric(),
    feat_class=character()
  ) %>%
    write.csv(file = paste0(output_folder, "classified_feats.csv"), row.names = FALSE)
  feat_id <- feature_df$feature[1]
}

# Perform classification ----
if(dataset_version%in%c("FT2040", "MS3000")){
  msdata <- readRDS(paste0(output_folder, "msdata.rds"))
} else {
  msdata <- readRDS(paste0(output_folder, "msdata_only_good.rds"))
}

while(TRUE){
  row_data <- feature_df[feature_df$feature==feat_id,]
  mzbounds <- c(row_data$min_mz, row_data$max_mz)
  rtbounds <- c(row_data$min_rt/60, row_data$max_rt/60)
  peakwidth <- diff(rtbounds)
  plotbounds <- c(min(rtbounds, min(rtbounds)-peakwidth/2), 
                  max(rtbounds, max(rtbounds)+peakwidth/2))
  eic <- msdata$EIC2[mz%between%mzbounds] %>% 
    # right_join(rt_corrections, by=c(rt="new_rt", "filename")) %>% 
    # mutate(int=ifelse(is.na(int), 0, int)) %>%
    filter(rt%between%(plotbounds)) %>%
    arrange(rt)
  dev.new(width=6, height=4)
  plot.new()
  plot.window(xlim=plotbounds, ylim=c(0, max(eic$int)), mar=c(1.1, 0.1, 0.1, 0.1))
  for(j in unique(eic$filename)){
    lines(eic[eic$filename==j, c("rt", "int")],
          col=file_data$col[file_data$filename==j],
          lwd=file_data$lwd[file_data$filename==j])
  }
  axis(1)
  title(paste(feat_id, round(mzbounds[1], 7)))
  abline(v=rtbounds, col="red")
  keyinput <- getGraphicsEvent(prompt = "", onKeybd = function(x){
    return(x)
  })
  dev.off()
  if(keyinput=="ctrl-["){ # aka Esc button
    break
  }
  if(keyinput=="ctrl-H"){
    # Open the file, remove the last row, write it back out
    data_classified <- read_csv(paste0(output_folder, "classified_feats.csv"), show_col_types = FALSE)
    feat_id <- tail(data_classified$feature, 1)
    data_classified <- data_classified[1:(nrow(data_classified)-1),]
    write_csv(data_classified, file = paste0(output_folder, "classified_feats.csv"))
    next
  }
  feat_class <- switch(
    keyinput,
    "Right" = "Bad",
    "Left" = "Good",
    "Up" = "Meh",
    "Down" = "Stans only"
  )
  cbind(row_data, feat_class) %>%
    write_csv(paste0(output_folder, "classified_feats.csv"), append = TRUE)
  feat_id <- feature_df$feature[which(feature_df$feature==feat_id)+1]
  if(is.na(feat_id)){
    print("Yay!")
    break
  }
}

# Visualization ----
data_classified <- read_csv(paste0(output_folder, "classified_feats.csv"), show_col_types = FALSE)
table(data_classified$feat_class)

ggplot(data_classified) +
  geom_rect(aes(xmin=min_rt/60, xmax=max_rt/60, ymin=min_mz, ymax=max_mz, color=feat_class),
            fill=NA, linewidth=1) +
  theme_bw() +
  labs(x="Retention time (in minutes)", y="m/z ratio", color="Classification")
plotly::ggplotly()

