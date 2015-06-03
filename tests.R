require(ggplot2)
require(gridExtra)
require(mp)

source("measures.R")

automated.m <- function(D, labels) {
  D.m <- D
  for (label in unique(labels)) {
    same.label <- labels == label
    D.m[same.label, same.label] <- D[same.label, same.label] * 0.1
    #D.m[same.label, diff.label] <- D[same.label, diff.label] * 10
    #D.m[diff.label, same.label] <- D.m[same.label, diff.label]
  }

  D.m
}

xy.df <- function(M) {
  M <- as.data.frame(M)
  names(M) <- c("x", "y")

  M
}

test <- function(file, suffix, output.dir) {
  cat("Testing dataset ", file, "...\n")
  dataset <- read.table(file)

  # Extract labels
  labels  <- dataset[, ncol(dataset)]

  # Remove labels from dataset
  X <- dataset[, -ncol(dataset)]
  
  n <- nrow(X)

  # Calculate distances (X) and normalize
  Dx <- dist(X)
  Dx <- Dx / mean(Dx)
  Dx <- as.matrix(Dx)

  sample.indices <- sample(n, 3*sqrt(n))
  Dx.s <- Dx[sample.indices, sample.indices]
  Ys   <- forceScheme(Dx.s)
  Ys   <- xy.df(Ys)
  Y    <- lamp(X, sample.indices, Ys)
  Y    <- xy.df(Y)

  # Plot mapping
  classes <- as.factor(labels)
  classes.s <- as.factor(labels[sample.indices])
  p.s <- ggplot(cbind(Ys, classes.s), aes(x = x, y = y, colour = classes.s)) + geom_point()
  p   <- ggplot(cbind(Y, classes), aes(x = x, y = y, colour = classes)) + geom_point()
  pdf(paste(output.dir, "original-", suffix, ".pdf", sep=""), width = 10, height = 5)
  grid.arrange(p.s, p,
               widths = unit(rep_len(3, 2), "null"),
               heights = unit(rep_len(1, 2), "null"),
               ncol=2)
  dev.off()
  png(paste(output.dir, "original-", suffix, ".png", sep=""), width = 1200, height = 600)
  grid.arrange(p.s, p,
               widths = unit(rep_len(3, 2), "null"),
               heights = unit(rep_len(1, 2), "null"),
               ncol=2)
  dev.off()

  # Calculate distances (Y) and normalize
  Dy <- dist(Y)
  Dy <- Dy / mean(Dy)
  Dy <- as.matrix(Dy)

  # Calculate measures and plot
  sigmas <- vector("numeric", n)
  sigmas[] <- 1
  P <- d2p(Dx, sigmas)
  Q <- d2p(Dy, sigmas)
  np = NP(Dx, Dy)
  #stress = stress(Dx, Dy),
  precision <- klDivergence(Q, P)
  recall <- klDivergence(P, Q)
  p.np        <- ggplot(cbind(Y, np), aes(x = x, y = y, colour = np)) + geom_point() + labs(title = "NP (9)")
  p.precision <- ggplot(cbind(Y, precision), aes(x = x, y = y, colour = precision)) + geom_point() + labs(title = "Precision")
  p.recall    <- ggplot(cbind(Y, recall), aes(x = x, y = y, colour = recall)) + geom_point()    + labs(title = "Recall")
  pdf(paste(output.dir, "measures-original-", suffix, ".pdf", sep=""), width = 15, height = 5)
  grid.arrange(p.np, p.precision, p.recall,
               widths = unit(rep_len(3, 3), "null"),
               heights = unit(rep_len(1, 3), "null"),
               ncol=3)
  dev.off()
  png(paste(output.dir, "measures-original-", suffix, ".png", sep=""), width = 1800, height = 600)
  grid.arrange(p.np, p.precision, p.recall,
               widths = unit(rep_len(3, 3), "null"),
               heights = unit(rep_len(1, 3), "null"),
               ncol=3)
  dev.off()

  # Perform manipulation
  Dx.m <- automated.m(Dx.s, labels[sample.indices])
  Ys.m <- forceScheme(Dx.m)
  Ys.m <- xy.df(Ys.m)
  Y.m  <- lamp(X, sample.indices, Ys.m)
  Y.m  <- xy.df(Y.m)

  # Plot mapping
  p.s <- ggplot(cbind(Ys.m, classes.s), aes(x = x, y = y, colour = classes.s)) + geom_point()
  p   <- ggplot(cbind(Y.m, classes), aes(x = x, y = y, colour = classes)) + geom_point()
  pdf(paste(output.dir, "manip-", suffix, ".pdf", sep=""), width = 10, height = 5)
  grid.arrange(p.s, p,
               widths = unit(rep_len(3, 2), "null"),
               heights = unit(rep_len(1, 2), "null"),
               ncol=2)
  dev.off()
  png(paste(output.dir, "manip-", suffix, ".png", sep=""), width = 1200, height = 600)
  grid.arrange(p.s, p,
               widths = unit(rep_len(3, 2), "null"),
               heights = unit(rep_len(1, 2), "null"),
               ncol=2)
  dev.off()

  # Calculate distances (Y.m) and normalize
  Dy <- dist(Y.m)
  Dy <- Dy / mean(Dy)
  Dy <- as.matrix(Dy)
  Q <- d2p(Dy, sigmas)

  # Calculate measures and plot
  np = np - NP(Dx, Dy)
  #stress = stress(Dx, Dy),
  precision <- precision - klDivergence(Q, P)
  recall <- recall - klDivergence(P, Q)
  p.np        <- ggplot(cbind(Y.m, np), aes(x = x, y = y, colour = np)) + geom_point() + labs(title = "NP (9)")
  p.precision <- ggplot(cbind(Y.m, precision), aes(x = x, y = y, colour = precision)) + geom_point() + labs(title = "Precision")
  p.recall    <- ggplot(cbind(Y.m, recall), aes(x = x, y = y, colour = recall)) + geom_point()    + labs(title = "Recall")
  pdf(paste(output.dir, "measures-manip-", suffix, ".pdf", sep=""), width = 15, height = 5)
  grid.arrange(p.np, p.precision, p.recall,
               widths = unit(rep_len(3, 3), "null"),
               heights = unit(rep_len(1, 3), "null"),
               ncol=3)
  dev.off()
  png(paste(output.dir, "measures-manip-", suffix, ".png", sep=""), width = 1800, height = 600)
  grid.arrange(p.np, p.precision, p.recall,
               widths = unit(rep_len(3, 3), "null"),
               heights = unit(rep_len(1, 3), "null"),
               ncol=3)
  dev.off()
}

test(file = "datasets/iris-std.tbl", suffix = "iris", "plots/")
test(file = "datasets/wdbc.tbl", suffix = "wdbc", "plots/")
test(file = "datasets/segmentation.tbl", suffix = "segmentation", "plots/")
test(file = "datasets/images.tbl", suffix = "images", "plots/")
