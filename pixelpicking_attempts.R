
# Setup ----
library(tidyverse)
library(RaMS)
library(xcms)
library(data.table)
dataset_version <- "FT2040"
# dataset_version <- "MS3000"
# dataset_version <- "CultureData"
# dataset_version <- "Pttime"
output_folder <- paste0("made_data_", dataset_version, "/")
msnexp_filled <- readRDS(paste0(output_folder, "msnexp_filled.rds"))
feature_centers <- featureDefinitions(msnexp_filled) %>%
  as.data.frame() %>%
  select(mzmed, rtmed) %>%
  rownames_to_column("feat_id") %>%
  mutate(rtmed=rtmed/60)
msdata <- readRDS(paste0(output_folder, "msdata.rds"))



# Extract the raw data and coerce to pixel matrix ----
interp_dt <- pbapply::pbmapply(function(mzmed_i, rtmed_i, feat_id_i){
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
  bind_rows() %>%
  group_by(feature, filename) %>%
  mutate(int=int/max(int, na.rm = TRUE)) %>%
  as.data.table()



# Perform the PCA and check variance explained ----
interp_complete <- interp_dt[, .(rt=rank(rt), int), by=c("feature", "filename")] %>%
  complete(feature, filename, rt, fill=list(int=-1))
interp_mat <- interp_complete %>%
  pivot_wider(names_from=feature, values_from = int) %>%
  # select(which(colSums(is.na(.))==0)) %>%
  arrange(filename, rt) %>%
  select(-rt, -filename) %>%
  data.matrix()
pcaoutput <- prcomp(interp_mat)


# Visualize in two and three dimensions ----
pcaoutput$rotation %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>%
  ggplot() +
  geom_text(aes(x=PC1, y=PC2, label=feature))

pcaoutput$x[,"PC1"] %>%
  matrix(nrow=50, ncol=length(fileNames(msnexp_filled))) %>%
  as.data.frame() %>%
  mutate(rownum=row_number()) %>%
  pivot_longer(-rownum, names_to="rt", values_to="int") %>% 
  mutate(rt=factor(rt, levels=paste0("V", 1:length(fileNames(msnexp_filled))))) %>%
  ggplot() +
  geom_tile(aes(x=rownum, y=rt, fill=int))

library(plotly)
pcaoutput$rotation %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>%
  plot_ly(x=~PC1, y=~PC2, z=~PC3, mode="markers", type="scatter3d", opacity=0.5)



# Fact-check a single feature ----
row_data <- feature_centers %>% filter(feat_id=="FT0025")
msdata_gp <- msdata$EIC2[mz%between%pmppm(row_data$mzmed, 5)] %>%
  filter(rt%between%(row_data$rtmed+c(-1, 1))) %>%
  ggplot() +
  geom_line(aes(x=rt, y=int, group=filename)) +
  ggtitle(paste0(row_data$feat_id, ": ", round(row_data$mzmed, 5)))
pixel_gp <- interp_dt %>%
  filter(feature==row_data$feat_id) %>%
  mutate(samp_type=str_extract(filename, "Smp|Std|Blk|Poo")) %>%
  ggplot() +
  geom_tile(aes(x=rt, y=filename, fill=int)) +
  facet_wrap(~samp_type, ncol=1, strip.position = "left", scales = "free_y") +
  theme(axis.text.y=element_blank())
lmat <- matrix(c(1,2,2), ncol = 1)
plot(gridExtra::arrangeGrob(msdata_gp, pixel_gp, layout_matrix = lmat))

msdata$EIC2[mz%between%pmppm(row_data$mzmed, 5)] %>%
  filter(rt%between%(row_data$rtmed+c(-1, 1))) %>%
  mutate(plot_color=ifelse(str_detect(filename, "8501"), "Croco", "Other")) %>%
  ggplot() +
  geom_line(aes(x=rt, y=int, group=filename, color=plot_color)) +
  ggtitle(paste0(row_data$feat_id, ": ", round(row_data$mzmed, 5), " (Trimethylamine)"))

# Write a small shiny app to choose features ----
library(shiny)
library(shinyjs)
ui <- function(){fluidPage(
  useShinyjs(),
  extendShinyjs(text = "shinyjs.closeWindow = function() { window.close(); }", 
                functions = c("closeWindow")),
  sidebarLayout(
    sidebarPanel(
      h3("Group peakpicker"),
      h4("Settings"),
      numericInput("n_kmeans_groups", label="Number of k-means groups",
                   value = 4, min = 1, max = 10, step = 1),
      numericInput("n_kmeans_dims", label="Number of PCs to use for k-means",
                   value = 3, min = 1, max = 10, step = 1),
      actionButton("kmeans_click", label = "Rerun k-means"),
      p(" "),
      plotOutput(outputId = "pcaprops", height = "200px"), 
      p(" "),
      actionButton("chosen_good", label = "Flag selection as Good"),
      actionButton("chosen_bad", label = "Flag selection as Bad"),
      actionButton("endsession", label = "Return to R and write out"),
      width = 3
    ),
    mainPanel(
      fluidRow(
        column(width=8, plotlyOutput(outputId = "plotlypca")),
        column(width=4, plotOutput(outputId = "kmeans_avgpeak"))
      ),
      fluidRow(
        column(width = 6, plotOutput(outputId = "live_peak", height = "200px")),
        column(width = 6, plotOutput(outputId = "avg_selected_peak", height = "200px"))
      ),
      width=9
    )
  )
)}
plotpeak <- function(feat_ids){
  plot_dt <- interp_dt[feature%chin%feat_ids][order(rt)][
    , .(rt=1:.N, int), by=c("feature", "filename")][
      , .(int=mean(int, na.rm=TRUE), iqr_int=IQR(int, na.rm=TRUE)), by=rt]
  plot_dt$rt <- rank(plot_dt$rt)
  plot_title <- ifelse(length(feat_ids)==1, feat_ids, "Aggregate")
  par(mar=c(0.1, 0.1, 1.1, 0.1))
  with(plot_dt, plot(rt, int, type="l", lwd=2, ylim=c(0, 1), 
                     xlab="", ylab="", main=plot_title))
  with(plot_dt, lines(rt, int+iqr_int))
  with(plot_dt, lines(rt, int-iqr_int))
}
server <- function(input, output, session){
  init_par <- par(no.readonly = TRUE)
  on.exit(par(init_par))
  output$pcaprops <- renderPlot({
    par(mar=c(2.1, 4.1, 0.1, 0.1))
    perc_exp <- pcaoutput$sdev^2/sum(pcaoutput$sdev^2)
    barplot(head(perc_exp*100, 10), 
            ylab = "% variance explained", names.arg = paste0("PC", 1:10))
    exp_thresholds <- c(0.2, 0.5, 0.8)
    PCs_to_explain_perc <- sapply(exp_thresholds, function(exp_threshold){
      which(cumsum(perc_exp)>exp_threshold)[1]
    })
    legend_text <- paste(paste0(exp_thresholds*100, "%: "), PCs_to_explain_perc, "PCs")
    legend("topright", legend = legend_text, bty='n', bg="transparent")
  })
  kmeaned_df <- reactive({
    input$kmeans_click
    pcaoutput$rotation[,1:input$n_kmeans_dims] %>%
      as.data.frame() %>%
      mutate(cluster=factor(kmeans(., centers=input$n_kmeans_groups)$cluster)) %>%
      rownames_to_column("feature")
  })
  output$plotlypca <- renderPlotly({
    gp <- kmeaned_df() %>%
      ggplot(aes(x=PC1, y=PC2, label=feature, color=cluster, key=feature)) +
      geom_text()
    ggplotly(gp, source = "plotlypca") %>% layout(dragmode="lasso")
  })
  output$kmeans_avgpeak <- renderPlot({
    clustergroups <- kmeaned_df() %>% 
      split(.$cluster) %>%
      lapply(`[[`, "feature")
    plot_dt <- merge(interp_dt, kmeaned_df())[order(rt)][
      , .(rt=1:.N, int), by=c("cluster", "feature", "filename")][
        , .(int=mean(int, na.rm=TRUE), iqr_int=IQR(int, na.rm=TRUE)), 
        by=c("rt", "cluster")]
    plot_dt %>%
      ggplot(aes(x=rt, color=cluster)) +
      geom_line(aes(y=int), linewidth=1) +
      geom_line(aes(y=int+iqr_int)) +
      geom_line(aes(y=int-iqr_int)) +
      facet_wrap(~cluster) +
      coord_cartesian(ylim=c(0, 1), clip="on") +
      theme_bw() +
      theme(legend.position = "none")
  })
  output$live_peak <- renderPlot({
    ed_hover <- event_data(source = "plotlypca", event = c("plotly_hover"))
    req(ed_hover)
    plotpeak(ed_hover$key)
  })
  output$avg_selected_peak <- renderPlot({
    ed_selected <- event_data(source = "plotlypca", event = c("plotly_selected"))
    req(ed_selected)
    plotpeak(ed_selected$key)
  })
  observeEvent(input$chosen_good, {
    good_feats <<- event_data(source = "plotlypca", event = "plotly_selected")$key
  })
  observeEvent(input$chosen_bad, {
    bad_feats <<- event_data(source = "plotlypca", event = "plotly_selected")$key
  })
  observeEvent(input$endsession, {
    message(paste("Number of good features:", length(good_feats)))
    message(paste("Number of bad features:", length(bad_feats)))
    js$closeWindow()
    stopApp()
  })
  session$onSessionEnded(function() {
    stopApp()
  })
}
shinyApp(ui=ui, server = server, options = c(launch.browser=TRUE))

# Write out chosen features ----

feature_centers %>%
  mutate(feat_class=case_when(
    feat_id%in%good_feats ~ "Good",
    feat_id%in%bad_feats ~ "Bad",
    TRUE~"Unclassified"
  )) %>%
  write.csv(file = paste0(output_folder, "quickclass_feats.csv"), row.names = FALSE)
