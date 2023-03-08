
# Setup ----
library(tidyverse)
library(RaMS)
library(xcms)
dataset_version <- "FT2040"
# dataset_version <- "MS3000"
# dataset_version <- "CultureData"
output_folder <- paste0("made_data_", dataset_version, "/")
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
interp_df <- pbapply::pbmapply(function(mzmed_i, rtmed_i, feat_id_i){
    interp_range <- rtmed_i+c(-0.5, 0.5)
    interp_points <- seq(interp_range[1], interp_range[2], length.out=50)
    msdata$EIC2[mz%between%pmppm(mzmed_i)] %>%
      split(.$filename) %>%
      lapply(function(eic_file){
        if(nrow(eic_file)>2){
          setNames(approx(eic_file$rt, eic_file$int, xout=interp_points), c("rt", "int"))
        } else {
          data.frame(rt=numeric(), int=numeric())
        }
      }) %>%
      bind_rows(.id="filename") %>%
      mutate(feature=feat_id_i)
  }, feature_centers$mzmed, feature_centers$rtmed, feature_centers$feat_id, 
  SIMPLIFY = FALSE) %>%
  group_by(feature, filename) %>%
  mutate(int=int/max(int))
interp_dt <- as.data.table(interp_df)

# Perform the PCA and check variance explained ----
pcaoutput <- interp_df %>%
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
pixel_gp <- interp_df %>%
  filter(feature==row_data$feat_id) %>%
  mutate(samp_type=str_extract(filename, "Smp|Std|Blk|Poo")) %>%
  ggplot() +
  geom_tile(aes(x=rt, y=filename, fill=int)) +
  facet_wrap(~samp_type, ncol=1, strip.position = "left", scales = "free_y") +
  theme(axis.text.y=element_blank())
plot(gridExtra::arrangeGrob(msdata_gp, pixel_gp, layout_matrix = matrix(c(
  1,2,2), ncol = 1)))

# Write a shiny app to choose features ----
library(shiny)
library(data.table)
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      h3("Circle a group to see an average peak from the group"),
      numericInput("n_init_groups", label="Number of initial groups",
                   value = 2, min = 1, max = 10, step = 1),
      numericInput("n_kmeans_dims", label="How many dimensions to use for k-means?",
                   value = 3, min = 2, max = 10, step = 1),
      actionButton("chosen_good", label = "Click to choose selection as Good")
    ),
    mainPanel(
      plotlyOutput(outputId = "plotly_pca"),
      plotOutput(outputId = "avg_selected_peak")
    )
  )
)
server <- function(input, output, session){
  output$plotly_pca <- renderPlotly({
    gp <- pcaoutput$rotation[,1:input$n_kmeans_dims] %>%
      as.data.frame() %>%
      mutate(cluster=factor(kmeans(., centers=input$n_init_groups)$cluster)) %>%
      rownames_to_column("feature") %>%
      ggplot(aes(x=PC1, y=PC2, label=feature, color=cluster, key=feature)) +
      geom_text()
    ggplotly(gp, source = "plotlypca") %>% layout(dragmode="lasso")
  })
  output$avg_selected_peak <- renderPlot({
    ed <- event_data(source = "plotlypca", event = "plotly_selected")
    req(ed)
    plot_dt <- interp_dt[
      feature%chin%ed$key, .(rt=rank(rt, ties.method = "min"), int), by=feature][
        , .(int=mean(int), iqr_int=IQR(int)), by=rt]
    plot_dt$rt <- rank(plot_dt$rt)
    with(plot_dt, plot(rt, int, type="l", lwd=2, ylim=c(0, 1)))
    with(plot_dt, lines(rt, int+iqr_int))
    with(plot_dt, lines(rt, int-iqr_int))
  })
}
shinyApp(ui=ui, server = server, options = c(launch.browser=TRUE))
