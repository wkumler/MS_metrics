

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


# Using the xcms defaults ----
features_xcms <- features_extracted %>%
  select(mean_mz, sd_ppm, mean_rt, sd_rt, mean_pw, log_mean_height, 
         sn, f, scale, lmin, feat_npeaks, sd_pw) %>%
  data.matrix()
rfmodel <- randomForest(features_xcms, y=factor(features_extracted$feat_class), importance=TRUE)
rfmodel
varImpPlot(rfmodel)

# Using and plotting the best 3 ----
rfmodel <- features_extracted %>%
  select(med_cor, med_SNR, shape_cor) %>%
  randomForest(y=factor(features_extracted$feat_class), importance=TRUE)
rfmodel
varImpPlot(rfmodel)
features_extracted %>%
  plotly::plot_ly(x=~log10(1-med_cor), y=~med_SNR, z=~shape_cor, color=~feat_class,
                  mode="markers", type= "scatter3d")



# rpart?
library(rpart)
fit <- features_extracted %>% 
  select(-feat_id) %>%
  select(mean_mz, sd_ppm, mean_rt, sd_rt, mean_pw, log_mean_height, 
         sn, f, lmin, scale, feat_npeaks, sd_pw, feat_class) %>%
  rpart(formula = feat_class~.)
library(rpart.plot)
pdf("rpartplot.pdf", width=6, height = 6)
rpart.plot(fit)
dev.off()


# Normal dimensionality reductions
features_extracted %>%
  column_to_rownames("feat_id") %>%
  select(-feat_class) %>%
  data.matrix() %>%
  scale() %>%
  prcomp() %>%
  biplot()
features_extracted %>%
  select(-feat_id, -feat_class) %>%
  data.matrix() %>%
  scale() %>%
  prcomp() %>%
  `$`("x") %>%
  as.data.frame() %>%
  cbind(feat_class=features_extracted$feat_class) %>%
  ggplot() +
  geom_point(aes(x=PC1, y=PC2, color=feat_class))
# plotly::plot_ly(x=~PC1, y=~PC2, z=~PC3, color=~feat_class,
#                 mode="markers", type="scatter3d")
