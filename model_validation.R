
# Setup ----
library(tidyverse)

FT2040_features <- read_csv("made_data_FT2040/features_extracted.csv")
MS3000_features <- read_csv("made_data_MS3000/features_extracted.csv")
Pttime_features <- read_csv("made_data_Pttime/features_extracted.csv")
CultureData_features <- read_csv("made_data_CultureData/features_extracted.csv")

# Visualize a single metric across all datasets
bind_rows(list(Falkor=FT2040_features, MESOSCOPE=MS3000_features,
               CultureData=CultureData_features, Pttime=Pttime_features), 
          .id = "cruise") %>%
  ggplot() + 
  geom_histogram(aes(x=med_cor, fill=feat_class), bins=40) +
  facet_wrap(~cruise, ncol = 1, scales = "free_y") +
  scale_fill_manual(breaks = c("Good", "Bad", "Meh", "Stans only", "Unclassified"),
                    values = c("#f4bb23", "#a41118", "#0b505c", "#028e34", "grey50"))

library(ggh4x)
bind_rows(list(Falkor=FT2040_features, MESOSCOPE=MS3000_features), .id = "cruise") %>%
  ggplot() + 
  geom_histogram(aes(x=med_SNR, fill=feat_class), bins=40) +
  facet_grid2(cruise~feat_class, scales = "free_y", independent = "y")



# Testing Falkor full model fit on MESOSCOPE data ----
falkor_full_model <- FT2040_features %>% 
  select(-feature) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)
FT2040_preds <- FT2040_features %>%
  select(-feature) %>%
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
  ggplot() + 
  geom_histogram(aes(x=pred_prob, fill=feat_class))


# Testing Falkor minimal model on MESOSCOPE ----
falkor_min_model <- FT2040_features %>% 
  select(feat_class, med_cor, med_SNR) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)
MS3000_features %>%
  mutate(pred_prob=predict(object=falkor_min_model, newdata = ., type = "response")) %>%
  mutate(pred_class=ifelse(pred_prob>0.5, "Good", "Bad")) %>%
  with(table(pred_class, feat_class))


# Modeling MESOSCOPE alone and comparing parameter differences ----
falkor_full_model <- FT2040_features %>% 
  select(-feature) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)
meso_full_model <- MS3000_features %>% 
  select(-feature) %>%
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
  select(feature, pred_prob, mean_mz, mean_rt)

given_peak <- "FT1955"
(row_data <- MS3000_features %>% 
  filter(feature==given_peak)) %>% 
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

# Constructing a model from BOTH datasets for comparison to each ----
# using only the minimal model of med_cor, med_SNR
falkor_min_model <- FT2040_features %>% 
  select(feat_class, med_cor, med_SNR) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)
FT2040_preds <- FT2040_features %>%
  select(-feature) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(pred_prob=predict(object=falkor_min_model,newdata = ., type = "response"))
meso_min_model <- MS3000_features %>% 
  select(feat_class, med_cor, med_SNR) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)
MS3000_preds <- MS3000_features %>%
  select(-feature) %>%
  mutate(pred_prob=predict(object=meso_min_model,newdata = ., type = "response"))

FM_min_model <- bind_rows(MS3000_features, FT2040_features) %>%
  select(feat_class, med_cor, med_SNR, sd_rt) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)
FM_preds <- bind_rows(MS3000_features, FT2040_features) %>%
  select(-feature) %>%
  mutate(pred_prob=predict(object=FM_min_model,newdata = ., type = "response")) %>%
  mutate(pred_class=pred_prob>0.5) %>%
  mutate(pred_class=ifelse(pred_class, "Good", "Bad"))
with(FM_preds, table(pred_class, feat_class))

list(falkor=falkor_min_model, meso=meso_min_model, both=FM_min_model) %>%
  sapply(broom::tidy, simplify = FALSE) %>%
  bind_rows(.id="model") %>%
  ggplot() +
  geom_point(aes(x=estimate, y=model, color=model)) +
  geom_linerange(aes(xmin=estimate-2*std.error, xmax=estimate+2*std.error, 
                     y=model, color=model)) +
  facet_wrap(~term, ncol=1, scales = "free")



# Using combined model to visualize predictions in all 4 datasets ----
FT2040_features <- read_csv("made_data_FT2040/features_extracted.csv")
MS3000_features <- read_csv("made_data_MS3000/features_extracted.csv")
both_min_model <- rbind(FT2040_features, MS3000_features) %>%
  select(feat_class, med_cor, med_SNR) %>%
  filter(feat_class%in%c("Good", "Bad")) %>%
  mutate(feat_class=feat_class=="Good") %>%
  glm(formula=feat_class~., family = binomial)
bind_rows(list(Falkor=FT2040_features, MESOSCOPE=MS3000_features,
               CultureData=CultureData_features, Pttime=Pttime_features), 
          .id = "cruise")

pred_probs_all <- bind_rows(lst(Falkor=FT2040_features, MESOSCOPE=MS3000_features,
                                CultureData=CultureData_features, Pttime=Pttime_features), 
                            .id = "dataset") %>%
  mutate(pred_prob=predict(object=both_min_model,newdata = ., type = "response")) %>%
  mutate(pred_class=ifelse(pred_prob>0.9, "Good", "Bad")) %>%
  mutate(error_type=case_when(
    feat_class=="Good" & pred_class=="Good" ~ "%TP",
    feat_class=="Good" & pred_class=="Bad" ~ "%FN",
    feat_class=="Bad" & pred_class=="Good" ~ "%FP",
    feat_class=="Bad" & pred_class=="Bad" ~ "%TN",
    TRUE ~ "Incomparable"
  ))
ggplot(pred_probs_all) +
  geom_histogram(aes(x=pred_prob, fill=feat_class), bins=50) +
  facet_wrap(~dataset, ncol = 1, scales="free_y") +
  coord_cartesian(ylim = c(0, 200)) +
  geom_vline(xintercept = 0.9, color="red", linewidth=1) +
  scale_fill_manual(breaks = c("Good", "Bad", "Meh", "Stans only", "Unclassified"),
                    values = c("#f4bb23", "#a41118", "#0b505c", "#028e34", "grey50"))

pred_probs_all %>%
  filter(error_type!="Incomparable") %>%
  mutate(dataset=factor(dataset, levels = c("Falkor", "MESOSCOPE", "CultureData", "Pttime"))) %>%
  with(table(feat_class, pred_class, dataset))

pred_probs_all %>%
  filter(error_type!="Incomparable") %>%
  group_by(dataset, error_type) %>%
  summarise(n=n(), .groups = 'drop_last') %>%
  mutate(perc=n/sum(n)) %>%
  ungroup() %>%
  complete(dataset, error_type, fill=list(n=0, perc=0)) %>%
  ggplot(aes(x=error_type, y=perc)) +
  geom_col(color="#0b505c", fill="#0b505c") +
  geom_text(aes(label=n), vjust=-0.2, color="#0b505c") +
  facet_wrap(~dataset, nrow=1) +
  scale_y_continuous(limits=c(0, 1)) +
  theme_bw() +
  theme(axis.title = element_blank())

pred_probs_all %>%
  filter(error_type!="Incomparable") %>%
  group_by(dataset, error_type) %>%
  summarise(n=n(), .groups = 'drop_last') %>%
  mutate(perc=n/sum(n)) %>%
  ungroup() %>%
  complete(dataset, error_type, fill=list(n=0)) %>%
  select(-perc) %>%
  pivot_wider(names_from=error_type, values_from = n) %>%
  mutate(`%FDR`=`%FP`/(`%FP`+`%TP`)*100) %>%
  mutate(`%GPF`=`%TP`/(`%FN`+`%TP`)*100) %>%
  mutate(`%GPF`=ifelse(dataset%in%c("CultureData", "Pttime"), NA, `%GPF`)) %>%
  pivot_longer(c(`%FDR`, `%GPF`)) %>%
  mutate(dataset=factor(dataset, levels = c("Falkor", "MESOSCOPE", "CultureData", "Pttime"))) %>%
  ggplot() +
  geom_point(aes(x=dataset, y=value, color=name), size=3) +
  facet_wrap(~name, scales = "free_y", nrow = 1) +
  scale_y_continuous(limits = c(0, NA))



