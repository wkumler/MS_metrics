

# Setup ----
library(tidyverse)
library(randomForest)

# dataset_version <- "FT350"
dataset_version <- "FT2040"
output_folder <- paste0("made_data_", dataset_version, "/")

features_extracted <- read_csv(paste0(output_folder, "features_extracted.csv")) %>%
  mutate(sn=ifelse(is.infinite(sn), 0, sn))

# Using ALL features up to this point ----
feature_subset <- features_extracted %>%
  select(-feat_id, -feat_class) %>%
  data.matrix()
rfmodel <- randomForest(feature_subset, y=factor(features_extracted$feat_class), 
                        importance=TRUE, proximity = TRUE, keep.forest = TRUE)
rfmodel
varImpPlot(rfmodel)
rf.mds <- MDSplot(rfmodel, factor(features_extracted$feat_class), k=3)
rf.mds$points %>%
  as.data.frame() %>%
  cbind(feat_class=features_extracted$feat_class) %>%
  plotly::plot_ly(x=~`Dim 1`, y=~`Dim 2`, z=~`Dim 3`, color=~feat_class)
# features_extracted %>%
#   select(-feat_id, -feat_class) %>%
#   data.matrix() %>%
#   partialPlot(x = rfmodel, x.var="med_SNR")
