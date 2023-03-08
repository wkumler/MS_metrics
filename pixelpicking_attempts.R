
# Setup ----
library(tidyverse)
library(RaMS)
classified_feats <- read_csv(paste0(output_folder, "classified_feats.csv")) %>%
  select(feature, feat_class)
msnexp_filled <- readRDS(paste0(output_folder, "msnexp_filled.rds"))
feature_centers <- featureDefinitions(msnexp_filled) %>%
  as.data.frame() %>%
  select(mzmed, rtmed) %>%
  rownames_to_column("feat_id") %>%
  mutate(rtmed=rtmed/60)

# Extract the raw data and coerce to pixel matrix ----
msdata <- readRDS(paste0(output_folder, "msdata.rds"))
interp_dt <- feature_centers %>%
  # slice(300:400) %>%
  pmap_dfr(function(...){
    row_data <- list(...)
    interp_range <- row_data$rtmed+c(-0.5, 0.5)
    interp_points <- seq(interp_range[1], interp_range[2], length.out=50)
    msdata$EIC2[mz%between%pmppm(row_data$mzmed)] %>%
      split(.$filename) %>%
      lapply(function(eic_file){
        if(nrow(eic_file)>2){
          setNames(approx(eic_file$rt, eic_file$int, xout=interp_points), c("rt", "int"))
        } else {
          data.frame(rt=numeric(), int=numeric())
        }
      }) %>%
      bind_rows(.id="filename") %>%
      mutate(feature=row_data$feat_id)
  }) %>%
  group_by(feature, filename) %>%
  mutate(int=int/max(int))

# Perform the PCA and check variance explained ----
pcaoutput <- interp_dt %>%
  group_by(feature, filename) %>%
  mutate(rt=rank(rt)) %>%
  # filter(feature%in%sprintf("FT%04d", 440:450)) %>%
  ungroup() %>%
  pivot_wider(names_from=feature, values_from = int) %>%
  select(which(colSums(is.na(.))==0)) %>%
  arrange(filename, rt) %>%
  select(-rt, -filename) %>%
  data.matrix() %>%
  prcomp()

plot(pcaoutput)
barplot(head(pcaoutput$sdev^2, 30))
round(sum(head(pcaoutput$sdev^2, 3))/sum(pcaoutput$sdev^2)*100)


# Visualize in two and three dimensions ----
pcaoutput$rotation %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>%
  left_join(classified_feats) %>%
  ggplot() +
  geom_text(aes(x=PC1, y=PC2, label=feature, color=feat_class))

pcaoutput$rotation %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>%
  left_join(classified_feats) %>%
  plot_ly(x=~PC1, y=~PC2, z=~PC3, color=~feat_class, mode="markers", type="scatter3d",
          opacity=0.5)


# Fact-check a single feature ----
row_data <- feature_centers %>% filter(feat_id=="FT1213")
msdata_gp <- msdata$EIC2[mz%between%pmppm(row_data$mzmed, 5)] %>%
  filter(rt%between%(row_data$rtmed+c(-1, 1))) %>%
  ggplot() +
  geom_line(aes(x=rt, y=int, group=filename)) +
  ggtitle(row_data$feat_id)
pixel_gp <- interp_dt %>%
  filter(feature==row_data$feat_id) %>%
  mutate(samp_type=str_extract(filename, "Smp|Std|Blk|Poo")) %>%
  ggplot() +
  geom_tile(aes(x=rt, y=filename, fill=int)) +
  facet_wrap(~samp_type, ncol=1, strip.position = "left", scales = "free_y") +
  theme(axis.text.y=element_blank())
plot(gridExtra::arrangeGrob(msdata_gp, pixel_gp, layout_matrix = matrix(c(
  1,2,2), ncol = 1)))

# I don't think the below is working?
hc <- hclust(dist(pcaoutput$rotation[,1:3]))
hc_colors <- classified_feats %>% 
  mutate(class_col=hcl.colors(6)[as.numeric(factor(feat_class))+1]) %>%
  pull(class_col, feature)
plot(hc, hang=-1)
points(hc$order, y=numeric(length(hc$order)), col=hc_colors[hc$labels[hc$order]], pch=19)
