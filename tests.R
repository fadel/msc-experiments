require(ggplot2)
require(gridExtra)
require(mp)

source("measures.R")

automated.m <- function(D, labels) {
  D.m <- D
  for (label in unique(labels)) {
    same.label <- labels == label
    D.m[same.label, same.label] <- D[same.label, same.label] * 0.1
  }

  D.m
}

test <- function(file, suffix, output.dir) {
  message("Testing dataset ", file)
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
  sigmas <- vector("numeric", n)
  sigmas[] <- 1
  P <- d2p(Dx, sigmas)
  Q <- d2p(Dy, sigmas)

  # Calculate measures
  message("\tCalculating measures")
  np <- NP(Dx, Dy)
  precision <- klDivergence(Q, P)
  recall <- klDivergence(P, Q)

  # Plot results
  #color_scale <- scale_colour_gradientn(colours = c("#e46c0a", "#dddddd", "#376092"))
  color_scale <- scale_colour_gradient2(mid = "#dddddd", space = "Lab")
  shape_scale <- scale_shape(solid = FALSE)
  Ys <- cbind(as.data.frame(Ys), classes.s)
  Y  <- cbind(as.data.frame(Y), classes, np, precision, recall)
  p.s <- ggplot(Ys) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes.s, colour = classes.s)) +
      shape_scale
  p <- ggplot(Y) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = classes)) +
      shape_scale
  p.np <- ggplot(Y) +
      theme_bw() +
      labs(x = "", y = "", title = "NP (9)") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = np)) +
      shape_scale + color_scale
  p.precision <- ggplot(Y) +
      theme_bw() +
      labs(x = "", y = "", title = "Precision") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = precision)) +
      shape_scale + color_scale
  p.recall <- ggplot(Y) +
      theme_bw() +
      labs(x = "", y = "", title = "Recall") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = recall)) +
      shape_scale + color_scale
  
  pdf(paste(output.dir, "original-", suffix, ".pdf", sep=""), width = 16, height = 8)
  grid.arrange(p.s, p, p.np, p.precision, p.recall, ncol=3)
  dev.off()

  # Perform manipulation
  message("\tCalculating Ys.m")
  Dx.m <- automated.m(Dx.s, labels[sample.indices])
  Ys.m <- forceScheme(Dx.m)

  # LAMP
  message("\tCalculating Y.m")
  Y.m <- lamp(X, sample.indices, Ys.m)

  # Calculate distances (Y.m) and normalize
  message("\tCalculating dist(Y.m)")
  Dy <- dist(Y.m)
  Dy <- Dy / mean(Dy)
  Dy <- as.matrix(Dy)

  message("\tCalculating Q")
  Q  <- d2p(Dy, sigmas)

  # Calculate measures
  message("\tCalculating measures")
  np <- np - NP(Dx, Dy)
  precision <- precision - klDivergence(Q, P)
  recall <- recall - klDivergence(P, Q)

  # Plot results
  Ys.m <- cbind(as.data.frame(Ys.m), classes.s)
  Y.m  <- cbind(as.data.frame(Y.m), classes, np, precision, recall)
  p.s <- ggplot(cbind(Ys.m, classes.s)) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes.s, colour = classes.s)) +
      scale_shape_identity() +
      shape_scale
  p <- ggplot(cbind(Y.m, classes)) +
      theme_bw() +
      labs(x = "", y = "") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = classes)) +
      shape_scale
  p.np <- ggplot(cbind(Y.m, np)) +
      theme_bw() +
      labs(x = "", y = "", title = "NP (9)") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = np)) +
      shape_scale + color_scale
  p.precision <- ggplot(cbind(Y.m, precision)) +
      theme_bw() +
      labs(x = "", y = "", title = "Precision") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = precision)) +
      shape_scale + color_scale
  p.recall <- ggplot(cbind(Y.m, recall)) +
      theme_bw() +
      labs(x = "", y = "", title = "Recall") +
      geom_point(aes(x = V1, y = V2, shape = classes, colour = recall)) +
      shape_scale + color_scale
  
  pdf(paste(output.dir, "manip-", suffix, ".pdf", sep=""), width = 16, height = 8)
  grid.arrange(p.s, p, p.np, p.precision, p.recall, ncol=3)
  dev.off()
}

test(file = "datasets/iris-std.tbl", suffix = "iris", "plots/")
test(file = "datasets/wdbc.tbl", suffix = "wdbc", "plots/")
test(file = "datasets/segmentation.tbl", suffix = "segmentation", "plots/")
test(file = "datasets/images.tbl", suffix = "images", "plots/")
