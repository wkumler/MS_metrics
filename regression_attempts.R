
# Setup ----
library(tidyverse)

# dataset_version <- "FT350"
dataset_version <- "FT2040"
output_folder <- paste0("made_data_", dataset_version, "/")

features_extracted <- read_csv(paste0(output_folder, "features_extracted.csv")) %>%
  mutate(sn=ifelse(is.infinite(sn), 0, sn))

set.seed(123)
traintestlist <- features_extracted %>%
  # select(-med_SNR, -med_cor, -med_missed_scans, -shape_cor, -area_cor) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=ifelse(feat_class=="Good", TRUE, FALSE)) %>%
  slice_sample(n = nrow(.)) %>%
  split(rbernoulli(nrow(.), 0.2)) %>%
  setNames(c("train", "test"))

# Single-(ML)-feature regression ----
single_model <- glm(feat_class~sd_rt, family = binomial, data=traintestlist$train)
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
xcms_model <- 

summary(full_model)
traintestlist$test %>%
  cbind(pred_class=predict(full_model, ., type="response")>0.9) %>%
  with(table(feat_class, pred_class)) %>%
  # caret::confusionMatrix() %>%
  print()



# Stepwise AIC optimization using MASS ----
all_model_feats <- select(traintestlist$train, -feat_id, -blank_found)
full_model <- glm(formula=feat_class~., family = binomial, data = all_model_feats)
step_model <- MASS::stepAIC(full_model, direction = "both")

setdiff(names(all_model_feats), broom::tidy(step_model)$term)



# Visualization (Single-metric, includes ggplot fit) ----
pred_feat <- "sd_rt"
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
  labs(x=pred_feat, y="Logistic likelihood", color="Classification")
cowplot::plot_grid(bar_gp, filled_gp, reg_gp, ncol = 1)
ggsave(paste0(pred_feat, ".png"), plot = last_plot(), path = "figures", 
       width = 5, height = 6, units = "in", device = "png")

# Visualization (Full model, excludes ggplot fit) ----
full_model <- traintestlist$train %>% 
  select(-feat_id, -blank_found) %>%
  glm(formula=feat_class~., family = binomial)
init_gp <- traintestlist$test %>%
  cbind(pred_prob=predict(full_model, ., type="response")) %>%
  mutate(pred_prob=cut(pred_prob, pretty(pred_prob, n=15), include.lowest=TRUE)) %>%
  mutate(feat_class=factor(feat_class, labels = c("Bad", "Good"), levels=c(FALSE, TRUE))) %>%
  ggplot(aes(x=pred_prob, fill=feat_class)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  labs(fill="Classification") +
  scale_x_discrete(drop=FALSE)
bar_gp <- init_gp + 
  geom_bar() + 
  labs(y="Count") + 
  theme(axis.text.x = element_blank()) +
  ggtitle("Full multiple regression model predictions on test set")
filled_gp <- init_gp + 
  geom_bar(position = "fill") + 
  labs(y="Proportion") +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5))
cowplot::plot_grid(bar_gp, filled_gp, ncol = 1)
ggsave("fullmodel_test.png", plot = last_plot(), path = "figures", 
       width = 8, height = 5, units = "in", device = "png")
