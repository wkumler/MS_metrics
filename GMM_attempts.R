library(ClusterR)
vals <- c(rnorm(10000), rnorm(1000, mean = 8, sd = 1))
n_centers <- 2
modelout <- GMM(matrix(vals, ncol = 1), gaussian_comps = n_centers)
gauss <- data.frame(centers=modelout$centroids, cov=modelout$covariance_matrices)
x_space <- seq(min(vals), max(vals), 0.01)
hist(vals)
for(i in 1:nrow(gauss)){
  curve <- dnorm(x_space, mean = gauss$centers[i], sd=sqrt(gauss$cov[i]))
  scaled_curve <- curve*length(vals)/n_centers
  lines(x_space, scaled_curve)
}
cluster_assignments <- ifelse(modelout$Log_likelihood[,1]>modelout$Log_likelihood[,2], "C1", "C2")
table(cluster_assignments)

FT2040_features <- read.csv("made_data_FT2040/features_extracted.csv")

vals <- c(FT2040_features$med_cor, 2+FT2040_features$med_cor*-1)
n_centers <- 3
modelout <- GMM(matrix(vals, ncol = 1), gaussian_comps = n_centers)
gauss <- data.frame(centers=modelout$centroids, cov=modelout$covariance_matrices)
x_space <- seq(min(vals), max(vals), 0.01)
hist(vals, breaks = 200, xlim = c(min(vals), max(vals)))
for(i in 1:nrow(gauss)){
  curve <- dnorm(x_space, mean = gauss$centers[i], sd=sqrt(gauss$cov[i]))
  scaled_curve <- curve*length(vals)/n_centers/50
  lines(x_space, scaled_curve, lwd=2)
}

vals <- c(FT2040_features$med_cor, 2+FT2040_features$med_cor*-1)
hist(vals, breaks = 100)
n_centers <- 3
modelout <- GMM(matrix(vals, ncol = 1), gaussian_comps = n_centers)
gauss <- data.frame(centers=modelout$centroids, cov=modelout$covariance_matrices)
gauss$name <- c("C1", "C2", "C3")
gauss[order(gauss$centers),]
x_space <- seq(min(vals), 1, 0.01)
histdata <- hist(vals[vals<1], breaks = 100, xlim = c(min(vals), 1))
binwidth <- unique(round(diff(histdata$breaks), 10))
auc <- table(apply(modelout$Log_likelihood, 1, which.max))
for(i in 1:nrow(gauss)){
  curve <- dnorm(x_space, mean = gauss$centers[i], sd=sqrt(gauss$cov[i]))
  scaled_curve <- curve*auc[i]*binwidth
  lines(x_space, scaled_curve, lwd=2)
}
curveheights <- sapply(1:nrow(gauss), function(i){
  dnorm(x_space, mean = gauss$centers[i], sd=sqrt(gauss$cov[i]))*auc[i]*binwidth
})[length(x_space):1,1:2]
v <- as.data.frame(apply(curveheights, 2, cumsum))
names(v) <- c("AUC2", "AUC1")
v$ratio <- v$AUC1/(v$AUC2+v$AUC1)
v$x_space <- rev(x_space)
cut_point <- v[which.min(abs(v$ratio-0.05)),]
abline(v=cut_point$x_space, col="red", lwd=2)

table(pred_class=ifelse(FT2040_features$med_cor>cut_point$x_space, "Good", "Bad"),
      real_class=FT2040_features$feat_class)
FT2040_features %>%
  ggplot() +
  geom_histogram(aes(x=med_cor, fill=feat_class), bins=50) +
  geom_vline(xintercept = cut_point$x_space, color="red", linewidth=1) +
  facet_wrap(~feat_class, ncol=1)
