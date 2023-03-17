
# Setup ----
library(tidyverse)
dataset_version <- "FT2040"
# dataset_version <- "MS3000"
# dataset_version <- "CultureData"
# dataset_version <- "Pttime"
output_folder <- paste0("made_data_", dataset_version, "/")
features_extracted <- read_csv(paste0(output_folder, "features_extracted.csv"))
classified_feats <- read_csv(paste0(output_folder, "classified_feats.csv")) %>%
  select(feature, feat_class)
quickclass_feats <- read_csv(paste0(output_folder, "quickclass_feats.csv")) %>%
  select(feature=feat_id, quick_class=feat_class)

# features_extracted %>%
#   left_join(quickclass_feats, by=c(feat_id="feature")) %>%
#   ggplot() +
#   geom_point(aes(x=med_cor, y=mean_pw, color=quick_class))


# Check performance against classified dataset ----
left_join(classified_feats, quickclass_feats) %>%
  with(table(feat_class, quick_class))


# Model comparison
# Model keeps blowing up with all features, reducing to a subset of "best"
full_model <- features_extracted %>%
  select(-feat_class) %>%
  left_join(classified_feats, by=c(feat_id="feature")) %>%
  select(-feat_id) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=ifelse(feat_class=="Bad", 0, 1)) %>%
  glm(feat_class~., family=binomial, data=.)
pred_probs <- features_extracted %>%
  mutate(pred_prob=predict(full_model, newdata=., type="response"))
pred_probs %>%
  # filter(feat_class%in%c("Good", "Bad")) %>%
  ggplot() +
  geom_histogram(aes(x=pred_prob, fill=feat_class), bins=50)
pred_probs %>%
  mutate(pred_class=ifelse(pred_prob>0.5, "Good", "Bad")) %>%
  with(table(pred_class, feat_class))

quick_model <- features_extracted %>%
  select(-feat_class) %>%
  left_join(quickclass_feats, by=c(feat_id="feature")) %>%
  select(-feat_id) %>%
  filter(quick_class%in%c("Good", "Bad")) %>%
  mutate(quick_class=ifelse(quick_class=="Bad", 0, 1)) %>%
  glm(quick_class~., family=binomial, data=.)
features_extracted %>%
  left_join(quickclass_feats, by=c(feat_id="feature")) %>%
  mutate(pred_prob=predict(quick_model, newdata=., type="response")) %>%
  ggplot() +
  geom_histogram(aes(x=pred_prob, fill=feat_class), bins=50)
features_extracted %>%
  left_join(quickclass_feats, by=c(feat_id="feature")) %>%
  mutate(pred_prob=predict(quick_model, newdata=., type="response")) %>%
  mutate(pred_class=ifelse(pred_prob>0.5, "Good", "Bad")) %>%
  with(table(pred_class, feat_class))

# Best features only below - dodging model explosion this way
full_model <- features_extracted %>%
  select(feat_class, med_cor, med_SNR, sd_rt) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=ifelse(feat_class=="Bad", 0, 1)) %>%
  glm(feat_class~., family=binomial, data=.)
pred_probs <- features_extracted %>%
  mutate(pred_prob=predict(full_model, newdata=., type="response"))
pred_probs %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  ggplot() +
  geom_histogram(aes(x=pred_prob, fill=feat_class), bins=50)
pred_probs %>%
  mutate(pred_class=ifelse(pred_prob>0.5, "Good", "Bad")) %>%
  with(table(pred_class, feat_class))

quick_model <- features_extracted %>%
  select(-feat_class) %>%
  left_join(quickclass_feats, by=c(feat_id="feature")) %>%
  select(c("med_cor", "med_SNR", "sd_rt", 
           # "mean_mz", "sd_ppm", "mean_rt", "mean_pw", "sd_pw",
           # "log_mean_height", "log_sd_height","sn", "f", "scale",
           # "lmin", "feat_npeaks", "n_found", "samps_found", "stans_found", 
           # "med_missed_scans", "smp_to_blk", "smp_to_std",
           # "shape_cor", "area_cor", 
           "quick_class")) %>%
  filter(quick_class%in%c("Good", "Bad")) %>%
  mutate(quick_class=ifelse(quick_class=="Bad", 0, 1)) %>%
  glm(quick_class~., family=binomial, data=.)
pred_probs <- features_extracted %>%
  mutate(pred_prob=predict(quick_model, newdata=., type="response")) %>%
  left_join(quickclass_feats, by=c(feat_id="feature"))
pred_probs %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  ggplot() +
  geom_histogram(aes(x=pred_prob, fill=feat_class), bins=50)
pred_probs %>%
  mutate(pred_class=ifelse(pred_prob>0.5, "Good", "Bad")) %>%
  with(table(pred_class, feat_class))

list(manual=broom::tidy(full_model), quick=broom::tidy(quick_model)) %>%
  bind_rows(.id="label_type") %>%
  ggplot(aes(x=term, color=label_type)) +
  geom_point(aes(y=estimate), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=estimate-2*std.error, ymax=estimate+2*std.error), 
                position = position_dodge(width = 0.5)) +
  facet_wrap(~term, ncol=5, scales = "free") +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept = 0, color="black")

