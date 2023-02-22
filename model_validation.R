
# Setup ----
library(tidyverse)

FT2040_features <- read_csv("made_data_FT2040/features_extracted.csv") %>%
  mutate(sn=ifelse(is.infinite(sn), 0, sn)) %>%
  select(-blank_found)
MS3000_features <- read_csv("made_data_MS3000/features_extracted.csv") %>%
  mutate(sn=ifelse(is.infinite(sn), 0, sn)) %>%
  select(-blank_found)


# Model fitting, both ----
falkor_full_model <- FT2040_features %>% 
  select(-feat_id) %>%
  # select(med_cor, med_SNR, feat_class) %>%
  # mutate(t_pval=t_pval/min(t_pval)) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)

MS3000_preds <- MS3000_features %>%
  # mutate(t_pval=t_pval/min(t_pval)) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(pred_class=predict(object=falkor_full_model,newdata = ., type = "response")>0.5) %>%
  mutate(pred_class=ifelse(pred_class, "Good", "Bad"))

table(MS3000_preds$pred_class, MS3000_preds$feat_class)

meso_full_model <- MS3000_features %>% 
  select(-feat_id, -log_mean_height, -log_med_cor) %>%
  # select(med_cor, med_SNR, feat_class) %>%
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


bind_rows(list(falkor=FT2040_features, meso=MS3000_features), .id="cruise") %>%
  mutate(t_pval=cut(t_pval, breaks = pretty(t_pval, n = 20))) %>%
  group_by(cruise, t_pval, feat_class) %>%
  summarise(n=log10(n())) %>%
  ggplot() +
  geom_col(aes(x=t_pval, y=n, fill=feat_class)) +
  facet_wrap(~cruise, ncol=1)


bind_rows(list(falkor=FT2040_features, meso=MS3000_features), .id="cruise") %>%
  group_by(cruise) %>%
  mutate(t_pval=t_pval/min(t_pval)) %>%
  # mutate(t_pval=10^t_pval) %>%
  ungroup() %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  mutate(feat_class=as.numeric(feat_class)) %>%
  ggplot(aes(x=t_pval, y=feat_class)) +
  geom_jitter(height = 0.05, width = 0) +
  geom_smooth(formula = y~x, method="glm", method.args = list(family=binomial)) +
  facet_wrap(~cruise, ncol=1)


# sapply(2:length(FT2040_features), function(i){
#   mod_subset <- FT2040_features %>%
#     filter(feat_class%in%c("Good", "Bad")) %>%
#     mutate(feat_class=feat_class=="Good") %>%
#     select(feat_class, any_of(names(FT2040_features)[2:i]))
#   falkor_subset_model <- mod_subset %>%
#     glm(formula=feat_class~., family = binomial)
#   MS3000_preds <- MS3000_features %>%
#     filter(feat_class%in%c("Good", "Bad")) %>%
#     mutate(pred_class=predict(object=falkor_subset_model, newdata = ., type = "response")>0.5) %>%
#     mutate(pred_class=ifelse(pred_class, "Good", "Bad"))
#   feat_added <- names(FT2040_features)[i]
#   n_false_pos <- table(c(MS3000_preds$pred_class, "Good", "Bad"), c(MS3000_preds$feat_class, "Bad", "Good"))[2,1]
#   data.frame(feat_added, n_false_pos)
# }) %>% t() %>% knitr::kable()
