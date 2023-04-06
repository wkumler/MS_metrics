
# Setup ----
library(tidyverse)

# dataset_version <- "FT350"
dataset_version <- "FT2040"
# dataset_version <- "MS3000"
output_folder <- paste0("made_data_", dataset_version, "/")

features_extracted <- read_csv(paste0(output_folder, "features_extracted.csv"))

set.seed(123)
traintestlist <- features_extracted %>%
  # select(-med_SNR, -med_cor, -med_missed_scans, -shape_cor, -area_cor) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=ifelse(feat_class=="Good", TRUE, FALSE)) %>%
  slice_sample(n = nrow(.)) %>%
  split(runif(nrow(.))<0.2) %>%
  setNames(c("train", "test"))

# Single-(ML)-feature regression ----
single_model <- glm(feat_class~log_mean_area, family = binomial, data=traintestlist$train)
summary(single_model)
traintestlist$test %>%
  cbind(pred_class=predict(single_model, ., type="response")>0.5) %>%
  with(table(c(feat_class, FALSE), c(pred_class, TRUE))) %>%
  caret::confusionMatrix() %>%
  print()



# All-feature regression ----
full_model <- traintestlist$train %>% 
  select(-feature) %>%
  glm(formula=feat_class~., family = binomial)

summary(full_model)
traintestlist$test %>%
  cbind(pred_class=predict(full_model, ., type="response")>0.5) %>%
  with(table(feat_class, pred_class)) %>%
  caret::confusionMatrix() %>%
  print()



# Stepwise AIC optimization using MASS ----
all_model_feats <- select(traintestlist$train, -feature)
full_model <- glm(formula=feat_class~., family = binomial, data = all_model_feats)
step_model <- MASS::stepAIC(full_model, direction = "both")

setdiff(names(all_model_feats), broom::tidy(step_model)$term)



# Minimal model that still has good performance ----
min_model_feats <- traintestlist$train %>%
  select(feat_class, max_cor, max_SNR, sd_rt)
min_model <- glm(formula=feat_class~., family = binomial, data = min_model_feats)
traintestlist$test %>%
  mutate(pred_prob=predict(object=min_model, newdata = ., type = "response")) %>%
  mutate(pred_class=round(pred_prob)) %>%
  with(table(pred_class, feat_class))



# Visualization (Single-metric, includes ggplot fit) ----
pred_feat <- "med_cor"
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

# Visualization (Two-metric with plotly viz) ----
dual_model <- traintestlist$train %>% 
  select(feat_class, med_SNR, med_cor) %>%
  glm(formula=feat_class~., family = binomial)
broom::tidy(dual_model)

model <- glm(formula = feat_class~med_cor+med_SNR, 
             data = traintestlist$train, 
             family = binomial)
b0 <- model$coefficients[1]
b1 <- model$coefficients[2]
b2 <- model$coefficients[3]

SNR_range <- seq(-5, 30, 1)
cor_range <- seq(-0.5, 1, 0.1)
log_curve_data <- expand_grid(med_SNR=SNR_range, med_cor=cor_range) %>%
  mutate(feat_class=1/(1+exp(-(b0+b1*med_cor+b2*med_SNR)))) %>%
  pivot_wider(names_from = med_cor, values_from = feat_class) %>%
  column_to_rownames("med_SNR") %>%
  data.matrix()
raw_curve_data <- traintestlist$train %>% 
  mutate(feat_class=as.numeric(feat_class)) %>%
  mutate(med_SNR=((med_SNR-min(med_SNR))/(max(med_SNR)-min(med_SNR)))) %>%
  mutate(med_SNR=med_SNR*length(SNR_range)) %>%
  mutate(med_cor=((med_cor-min(med_cor))/(max(med_cor)-min(med_cor)))) %>%
  mutate(med_cor=med_cor*length(cor_range))
library(plotly)
plot_ly() %>%
  add_trace(x=~med_cor, y=~med_SNR, z=~feat_class, opacity=0.5,
            mode="markers", type="scatter3d", 
            data=raw_curve_data, marker=list(color="black")) %>%
  add_surface(z=log_curve_data)

traintestlist$test %>%
  cbind(pred_class=predict(dual_model, ., type="response")>0.5) %>%
  with(table(feat_class, pred_class)) %>%
  caret::confusionMatrix() %>%
  print()



# Visualization (Full model, excludes ggplot fit) ----
full_model <- traintestlist$train %>% 
  select(-feature) %>%
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

