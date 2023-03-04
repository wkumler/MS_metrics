
# Setup ----
library(tidyverse)

FT2040_features <- read_csv("made_data_FT2040/features_extracted.csv") %>%
  mutate(sn=ifelse(is.infinite(sn), 0, sn)) %>%
  select(-blank_found)
MS3000_features <- read_csv("made_data_MS3000/features_extracted.csv") %>%
  mutate(sn=ifelse(is.infinite(sn), 0, sn)) %>%
  select(-blank_found)
# Visualize a single metric across both datasets
bind_rows(list(Falkor=FT2040_features, MESOSCOPE=MS3000_features), .id = "cruise") %>%
  drop_na(-shape_cor, -area_cor) %>%
  ggplot() + 
  geom_histogram(aes(x=max_SNR, fill=feat_class), bins=40) +
  facet_wrap(~cruise, ncol = 1, scales = "free_y")
library(ggh4x)
bind_rows(list(Falkor=FT2040_features, MESOSCOPE=MS3000_features), .id = "cruise") %>%
  ggplot() + 
  geom_histogram(aes(x=t_pval, fill=feat_class), bins=40) +
  facet_grid2(cruise~feat_class, scales = "free_y", independent = "y")

# Testing Falkor full model fit on MESOSCOPE data ----
falkor_full_model <- FT2040_features %>% 
  select(-feat_id, -area_cor, -shape_cor) %>%
  drop_na() %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)
FT2040_preds <- FT2040_features %>%
  select(-feat_id, -area_cor, -shape_cor) %>%
  drop_na() %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(pred_prob=predict(object=falkor_full_model,newdata = ., type = "response")) %>%
  mutate(pred_class=pred_prob>0.5) %>%
  mutate(pred_class=ifelse(pred_class, "Good", "Bad"))
table(predicted=FT2040_preds$pred_class, actual=FT2040_preds$feat_class)
FT2040_preds %>% 
  ggplot() + 
  geom_histogram(aes(x=pred_prob, fill=feat_class))

MS3000_preds <- MS3000_features %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(pred_prob=predict(object=falkor_full_model,newdata = ., type = "response")) %>%
  mutate(pred_class=pred_prob>0.5) %>%
  mutate(pred_class=ifelse(pred_class, "Good", "Bad"))
table(predicted=MS3000_preds$pred_class, actual=MS3000_preds$feat_class)
MS3000_preds %>% 
  drop_na(pred_prob) %>%
  ggplot() + 
  geom_histogram(aes(x=pred_prob, fill=feat_class))


# Testing Falkor minimal model on MESOSCOPE ----
falkor_min_model <- FT2040_features %>% 
  select(feat_class, med_cor, med_SNR, sd_rt) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)
MS3000_features %>%
  mutate(pred_prob=predict(object=falkor_min_model, newdata = ., type = "response")) %>%
  mutate(pred_class=round(pred_prob)) %>%
  with(table(pred_class, feat_class))


# Modeling MESOSCOPE alone and comparing parameter differences ----
falkor_full_model <- FT2040_features %>% 
  select(-feat_id, -area_cor, -shape_cor) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)
meso_full_model <- MS3000_features %>% 
  select(-feat_id, -area_cor, -shape_cor) %>%
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


# Investigating specific peaks ----
library(RaMS)
msdata <- readRDS("made_data_MS3000/msdata.rds")
file_types <- c("Blk", "15m", "DCM", "175m", "Poo", "Std")

MS3000_preds %>%
  filter(feat_class=="Good") %>%
  arrange(pred_prob) %>% 
  select(feat_id, pred_prob, mean_mz, mean_rt)

given_peak <- "FT1955"
(row_data <- MS3000_features %>% 
  filter(feat_id==given_peak)) %>% 
  as.data.frame() %>% t()
msdata$EIC2[mz%between%pmppm(row_data$mean_mz, 5)][
  rt%between%(row_data$mean_rt/60+c(-1, 1))
  ] %>%
  mutate(line_color=str_extract(filename, paste(file_types, collapse="|"))) %>%
  mutate(line_color=factor(line_color, levels=file_types)) %>%
  # filter(line_color=="Poo") %>%
  filter(line_color!="Std") %>%
  ggplot() +
  geom_line(aes(x=rt, y=int, group=filename, color=line_color)) +
  # facet_wrap(~filename) +
  scale_y_continuous()
