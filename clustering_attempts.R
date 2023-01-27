

# Setup ----

library(tidyverse)
library(dbscan)

features_extracted <- read_csv("made_data/features_extracted.csv")
feature_mat <- features_extracted %>%
  select(where(is.numeric)) %>%
  select(-sn) %>%
  na.omit() %>%
  data.matrix() %>%
  scale()

res.dbscan <- dbscan::dbscan(feature_mat, eps = 3, minPts = 5)
factoextra::fviz_cluster(res.dbscan, feature_mat, geom = "point")

features_extracted %>%
  mutate(db_clust=res.dbscan$cluster) %>%
  na.omit() %>%
  # with(table(feat_class, db_clust))
  # filter(db_clust==0, feat_class=="Good")
  ggplot() +
  geom_jitter(aes(x=feat_class, y=db_clust), height = 0.1, width = 0.1)

feature_subset <- features_extracted %>%
  select(med_cor, med_SNR,n_found, med_missed_scans) %>%
  data.matrix() %>%
  na.omit()
