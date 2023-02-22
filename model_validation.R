
# Setup ----
library(tidyverse)

FT2040_features <- read_csv("made_data_FT2040/features_extracted.csv") %>%
  mutate(sn=ifelse(is.infinite(sn), 0, sn)) %>%
  select(-blank_found)
MS3000_features <- read_csv("made_data_MS3000/features_extracted.csv") %>%
  mutate(sn=ifelse(is.infinite(sn), 0, sn)) %>%
  select(-blank_found)
# Visualize a single feature across both datasets
bind_rows(list(Falkor=FT2040_features, MESOSCOPE=MS3000_features), .id = "cruise") %>%
  ggplot() + geom_histogram(aes(x=t_pval, fill=feat_class)) +
  facet_wrap(~cruise, ncol = 1, scales = "free_y")

# Testing Falkor model fit on MESOSCOPE data ----
falkor_full_model <- FT2040_features %>% 
  select(-feat_id) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)
FT2040_preds <- FT2040_features %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(pred_prob=predict(object=falkor_full_model,newdata = ., type = "response")) %>%
  mutate(pred_class=pred_prob>0.5) %>%
  mutate(pred_class=ifelse(pred_class, "Good", "Bad"))
table(predicted=FT2040_preds$pred_class, actual=FT2040_preds$feat_class)
FT2040_preds %>% ggplot() + geom_histogram(aes(x=pred_prob, fill=feat_class))

MS3000_preds <- MS3000_features %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(pred_prob=predict(object=falkor_full_model,newdata = ., type = "response")) %>%
  mutate(pred_class=pred_prob>0.5) %>%
  mutate(pred_class=ifelse(pred_class, "Good", "Bad"))
table(predicted=MS3000_preds$pred_class, actual=MS3000_preds$feat_class)
MS3000_preds %>% ggplot() + geom_histogram(aes(x=pred_prob, fill=feat_class))


# Modeling MESOSCOPE alone and comparing parameter differences ----
falkor_full_model <- FT2040_features %>% 
  select(-feat_id) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)
meso_full_model <- MS3000_features %>% 
  select(-feat_id) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)

bind_rows(
  list(
    falkor=broom::tidy(falkor_full_model),
    meso=broom::tidy(meso_full_model)
  ), .id="cruise"
) %>%
  ggplot(aes(y=term, color=cruise)) +
  geom_point(aes(x=estimate)) +
  geom_errorbar(aes(xmin=estimate-2*std.error, xmax=estimate+2*std.error)) +
  facet_wrap(~term, ncol=5, scales = "free") +
  theme(axis.text.y = element_blank())
ggsave("model_comparison.pdf", device = "pdf", path = "figures", 
       width = 9, height = 6)


# Subsetting to figure out max performance ----
meso_sub_model <- MS3000_features %>% 
  select(-feat_id, -log_mean_height, -log_med_cor) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)
falkor_sub_model <- FT2040_features %>% 
  select(-feat_id, -log_mean_height, -log_med_cor) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)
bind_rows(
  list(
    falkor=broom::tidy(falkor_sub_model),
    meso=broom::tidy(meso_sub_model)
  ), .id="cruise"
) %>%
  ggplot(aes(y=term, color=cruise)) +
  geom_point(aes(x=estimate)) +
  geom_errorbar(aes(xmin=estimate-2*std.error, xmax=estimate+2*std.error)) +
  facet_wrap(~term, ncol=5, scales = "free") +
  theme(axis.text.y = element_blank())
