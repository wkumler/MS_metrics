
# Setup ----
library(tidyverse)

# dataset_version <- "FT350"
dataset_version <- "FT2040"
output_folder <- paste0("made_data_", dataset_version, "/")

features_extracted <- read_csv(paste0(output_folder, "features_extracted.csv")) %>%
  mutate(sn=ifelse(is.infinite(sn), 0, sn))

traintestlist <- features_extracted %>%
  # select(-med_SNR, -med_cor, -med_missed_scans, -shape_cor, -area_cor) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=ifelse(feat_class=="Good", TRUE, FALSE)) %>%
  slice_sample(n = nrow(.)) %>%
  split(rbernoulli(nrow(.), 0.2)) %>%
  setNames(c("train", "test"))



# Single-(ML)-feature regression ----
pred_feat <- "f"
init_gp <- traintestlist$train %>%
  mutate(pred_cut=cut(get(pred_feat), pretty(get(pred_feat), n=15), include.lowest=TRUE)) %>%
  mutate(feat_class=factor(feat_class, labels = c("Bad", "Good"), levels=c(FALSE, TRUE))) %>%
  ggplot(aes(x=pred_cut, fill=feat_class)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  labs(fill="Classification")
bar_gp <- init_gp + 
  geom_bar() + 
  labs(y="Count") + 
  theme(axis.text.x = element_blank())
filled_gp <- init_gp + 
  geom_bar(position = "fill") + 
  labs(y="Proportion") +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5))
reg_gp <- traintestlist$train %>%
  ggplot(aes(x=get(pred_feat), y=as.numeric(feat_class))) +
  geom_jitter(height = 0.01, alpha=0.2, aes(color=feat_class)) +
  geom_smooth(method="glm", method.args = list(family=binomial), formula='y~x',
              color="black") +
  theme_bw() +
  labs(x="Metric values", y="Logistic likelihood", color="Classification")
cowplot::plot_grid(bar_gp, filled_gp, reg_gp, ncol = 1)
ggsave(paste0(pred_feat, ".png"), plot = last_plot(), path = "figures", 
       width = 5, height = 6, units = "in", device = "png")

single_model <- glm(feat_class~get(pred_feat), family = binomial, data=traintestlist$train)
summary(single_model)
traintestlist$test %>%
  cbind(pred_class=predict(single_model, ., type="response")>0.5) %>%
  with(table(feat_class, pred_class)) %>%
  # caret::confusionMatrix() %>%
  print()


# All-feature regression ----
full_model <- traintestlist$train %>% 
  select(-feat_id, -blank_found) %>%
  glm(formula=feat_class~., family = binomial)
summary(full_model)
traintestlist$test %>%
  cbind(pred_class=predict(full_model, ., type="response")>0.5) %>%
  with(table(feat_class, pred_class)) %>%
  # caret::confusionMatrix() %>%
  print()

# Just-best regression
custom_model <- traintestlist$train %>% 
  select(feat_class, med_cor, med_SNR) %>%
  glm(formula=feat_class~., family = binomial)
summary(custom_model)
traintestlist$test %>%
  mutate(pred_class=predict(custom_model, ., type="response")>=0.5) %>%
  with(table(feat_class, pred_class)) %>%
  # caret::confusionMatrix() %>%
  print()

# Just xcms regression
xcms_model <- traintestlist$train %>%
  select(mean_mz, sd_ppm, mean_rt, sd_rt, mean_pw, sd_pw,
         log_mean_height, log_sd_height, sn, f, scale, 
         lmin, feat_class) %>%
  glm(formula=feat_class~., family = binomial)
summary(xcms_model)
traintestlist$test %>%
  mutate(pred_class=predict(xcms_model, ., type="response")>=0.5) %>%
  with(table(feat_class, pred_class)) %>%
  # caret::confusionMatrix() %>%
  print()
