require(ggplot2)
require(gridExtra)
require(mp)

source("measures.R")

#automated.m <- function(D, labels) {
#  D.m <- D
#  for (label in unique(labels)) {
#    same.label <- labels == label
#    D.m[same.label, same.label] <- D[same.label, same.label] * 0.1
#  }
#
#  D.m
#}
automated.m <- function(Xs, labels) {
  n <- nrow(Xs)
  p <- ncol(Xs)
  Xs <- cbind(Xs, matrix(data=0, nrow=n, ncol=p))
  for (label in unique(labels)) {
    for (j in 1:p) {
      Xs[labels == label, j + p] <- mean(Xs[labels == label, j])
    }
  }

  dist(Xs)
}

color_scale.blue_orange <- function(name) {
  scale_colour_gradient(name = name, high = "#376092", low = "#e46c0a", space = "Lab")
}

color_scale.gradient2 <- function(name) {
  scale_colour_gradient2(name = name, mid = "#dddddd", space = "Lab")
}

test <- function(file, suffix, output.dir) {
  message("Testing dataset: ", file)
  dataset <- read.table(file)

  # Extract labels
  labels <- dataset[, ncol(dataset)]
  classes <- as.factor(labels)
  X <- dataset[, -ncol(dataset)]

  n <- nrow(X)

  # Calculate distances (X) and normalize
  message("\tCalculating dist(X)")
  Dx <- dist(X)
  Dx <- Dx / mean(Dx)
  Dx <- as.matrix(Dx)

  # Sample dataset
  sample.indices <- sample(n, 3*sqrt(n))
  classes.s <- as.factor(labels[sample.indices])

  # Automatic sample positioning
  message("\tCalculating Ys")
  Dx.s <- Dx[sample.indices, sample.indices]
  Ys   <- forceScheme(Dx.s)

  # LAMP
  message("\tCalculating Y")
  Y <- lamp(X, sample.indices, Ys)

  # Calculate distances (Y) and normalize
  message("\tCalculating dist(Y)")
  Dy <- dist(Y)
  Dy <- Dy / mean(Dy)
  Dy <- as.matrix(Dy)

  message("\tCalculating P and Q")
  prob <- d2p(Dx^2)
  P <- prob$P
  Q <- d2p.beta(Dy^2, prob$beta)

  # Calculate measures
  message("\tCalculating measures")
  np <- NP(Dx, Dy)
  silh <- silhouette(Dy, classes)
  precision <- klDivergence(Q, P)
  recall <- klDivergence(P, Q)

  measures <- rbind(data.frame(mean=mean(np), median=median(np), sd=sd(np)),
                    data.frame(mean=mean(silh), median=median(silh), sd=sd(silh)),
                    data.frame(mean=mean(precision), median=median(precision), sd=sd(precision)),
                    data.frame(mean=mean(recall), median=median(recall), sd=sd(recall)))
  write.table(measures, paste(output.dir, suffix, "-measures.csv", sep=""), row.names=F)

  if (!(all(is.finite(np)) &&
      all(is.finite(silh)) &&
      all(is.finite(precision)) &&
      all(is.finite(recall)))) {
    stop("Non-finite measures found")
  }

  # Plot results
  message("\tPlotting results")
  shape_scale <- scale_shape_manual(name = "Classe", values = 1:nlevels(classes))
  Ys <- cbind(as.data.frame(Ys), classes.s)
  Y  <- cbind(as.data.frame(Y), classes, np, silh, precision, recall)
  p.s <- ggplot(Ys) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes.s, colour = classes.s)) +
      shape_scale + scale_color_manual(name = "Classe", values = 1:nlevels(classes))
  ggsave(paste(output.dir, "subsample-", suffix, ".pdf", sep=""), p.s, width=5, height=5)

  p <- ggplot(Y) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = classes)) +
      shape_scale + scale_color_manual(name = "Classe", values = 1:nlevels(classes))
  ggsave(paste(output.dir, suffix, ".pdf", sep=""), p, width=5, height=5)

  p.np <- ggplot(Y) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = np)) +
      shape_scale + color_scale.blue_orange("NP")
  ggsave(paste(output.dir, "np-", suffix, ".pdf", sep=""), p.np, width=5, height=5)

  p.silh <- ggplot(Y) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = silh)) +
      shape_scale + color_scale.blue_orange("Silhueta")
  ggsave(paste(output.dir, "silh-", suffix, ".pdf", sep=""), p.silh, width=5, height=5)

  p.precision <- ggplot(Y) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = precision)) +
      shape_scale + color_scale.blue_orange("Precisão")
  ggsave(paste(output.dir, "precision-", suffix, ".pdf", sep=""), p.precision, width=5, height=5)

  p.recall <- ggplot(Y) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = recall)) +
      shape_scale + color_scale.blue_orange("Revocação")
  ggsave(paste(output.dir, "recall-", suffix, ".pdf", sep=""), p.recall, width=5, height=5)

  
  pdf(paste(output.dir, "all-", suffix, ".pdf", sep=""), width = 12, height = 16)
  grid.arrange(p.s, p, p.np, p.silh, p.precision, p.recall, ncol = 2)
  dev.off()

  # Perform manipulation
  message("\tCalculating Ys.m")
  Dx.m <- automated.m(X[sample.indices, ], labels[sample.indices])
  Ys.m <- forceScheme(Dx.m, Ys[, 1:2])

  # LAMP
  message("\tCalculating Y.m")
  Y.m <- lamp(X, sample.indices, Ys.m)

  # Calculate distances (Y.m) and normalize
  message("\tCalculating dist(Y.m)")
  Dy <- dist(Y.m)
  Dy <- Dy / mean(Dy)
  Dy <- as.matrix(Dy)

  message("\tCalculating Q")
  Q <- d2p.beta(Dy^2, prob$beta)

  # Calculate measures
  message("\tCalculating measures")
  np <- NP(Dx, Dy) - np
  silh <- silhouette(Dy, classes) - silh
  precision <- klDivergence(Q, P) - precision
  recall <- klDivergence(P, Q) - recall

  measures <- rbind(data.frame(mean=mean(np), median=median(np), sd=sd(np)),
                    data.frame(mean=mean(silh), median=median(silh), sd=sd(silh)),
                    data.frame(mean=mean(precision), median=median(precision), sd=sd(precision)),
                    data.frame(mean=mean(recall), median=median(recall), sd=sd(recall)))
  write.table(measures, paste(output.dir, suffix, "-manip-measures.csv", sep=""), row.names=F)

  if (!(all(is.finite(np)) &&
      all(is.finite(silh)) &&
      all(is.finite(precision)) &&
      all(is.finite(recall)))) {
    stop("Non-finite measures found")
  }

  # Plot results
  message("\tPlotting results")
  Ys.m <- cbind(as.data.frame(Ys.m), classes.s)
  Y.m  <- cbind(as.data.frame(Y.m), classes, np, silh, precision, recall)
  pm.s <- ggplot(Ys.m) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes.s, colour = classes.s)) +
      shape_scale + scale_color_manual(name = "Classe", values = 1:nlevels(classes))
  ggsave(paste(output.dir, "manip-subsample-", suffix, ".pdf", sep=""), pm.s, width=5, height=5)

  pm <- ggplot(Y.m) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = classes)) +
      shape_scale + scale_color_manual(name = "Classe", values = 1:nlevels(classes))
  ggsave(paste(output.dir, "manip-", suffix, ".pdf", sep=""), pm, width=5, height=5)

  p.np <- ggplot(Y.m) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = np)) +
      shape_scale + color_scale.gradient2("NP")
  ggsave(paste(output.dir, "manip-np-", suffix, ".pdf", sep=""), p.np, width=5, height=5)

  p.silh <- ggplot(Y.m) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = silh)) +
      shape_scale + color_scale.gradient2("Silhueta")
  ggsave(paste(output.dir, "manip-silh-", suffix, ".pdf", sep=""), p.silh, width=5, height=5)

  p.precision <- ggplot(Y.m) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = precision)) +
      shape_scale + color_scale.gradient2("Precisão")
  ggsave(paste(output.dir, "manip-precision-", suffix, ".pdf", sep=""), p.precision, width=5, height=5)

  p.recall <- ggplot(Y.m) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = recall)) +
      shape_scale + color_scale.gradient2("Revocação")
  ggsave(paste(output.dir, "manip-recall-", suffix, ".pdf", sep=""), p.recall, width=5, height=5)

  pdf(paste(output.dir, "original-manip-", suffix, ".pdf", sep=""), width = 10, height = 8)
  grid.arrange(p.s, p, pm.s, pm, ncol = 2)
  dev.off()

  pdf(paste(output.dir, "manip-measures-", suffix, ".pdf", sep=""), width = 10, height = 8)
  grid.arrange(p.np, p.silh, p.precision, p.recall, ncol = 2)
  dev.off()
  
  pdf(paste(output.dir, "manip-all-", suffix, ".pdf", sep=""), width = 12, height = 16)
  grid.arrange(pm.s, pm, p.np, p.silh, p.precision, p.recall, ncol = 2)
  dev.off()
}

#test(file = "datasets/iris.tbl", suffix = "iris", "tests/")
test(file = "datasets/wdbc.tbl", suffix = "wdbc", "tests/")
test(file = "datasets/segmentation.tbl", suffix = "segmentation", "tests/")
test(file = "datasets/glass.tbl", suffix = "glass", "tests/")
