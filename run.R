# run.R
#
# Main experiments script.

library(cluster)
library(logging)
library(MASS)
library(mp)
library(Rtsne)

source("measures.R")
source("util.R")

# Performs automated silhouette improvement manipulation, using a method
# inspired by Schaefer et al. (2013).
automated.silh <- function(Xs, labels) {
  n <- nrow(Xs)
  p <- ncol(Xs)
  Xs <- cbind(Xs, matrix(data=0, nrow=n, ncol=p))
  for (label in unique(labels)) {
    for (j in 1:p) {
      Xs[labels == label, j + p] <- mean(Xs[labels == label, j])
    }
  }

  Dx <- dist(Xs)
  # Dx <- Dx / mean(Dx)
  as.matrix(Dx)
}

# NOTE: This function requires the 'klmeasure' binary from:
# http://research.cs.aalto.fi/pml/software/dredviz/
nerv <- function(Dx, Y, lambda=0.1) {
  # Create SOM_PAK file for Dx
  Dx.fname <- tempfile()
  Dx.f <- file(Dx.fname, "w")
  cat(sprintf("%d\n", ncol(Dx)), file=Dx.f)
  write.table(Dx, Dx.f, col.names=F, row.names=F)
  close(Dx.f)

  # Create SOM_PAK file for Y
  Y.fname <- tempfile()
  Y.f <- file(Y.fname, "w")
  cat(sprintf("%d\n", ncol(Y)), file=Y.f)
  write.table(Y, Y.f, col.names=F, row.names=F)
  close(Y.f)

  # Run NeRV
  Ym.fname <- tempfile()
  system2("./nerv",
          stdout=F,
          stderr=F,
          args=c("--inputdist", Dx.fname,
                 "--outputfile", Ym.fname,
                 "--init", Y.fname,
                 "--lambda", sprintf("%.2f", lambda)))

  # Read results from generated file; remove file afterwards
  Ym <- read.table(Ym.fname, skip=1)
  file.remove(Dx.fname, Y.fname, Ym.fname)

  Ym
}

# Wrapper so that we can 'do.call' pekalska as we do with other techniques
pekalska.wrapper <- function(X, sample.indices, Ys) {
  pekalska(dist(X), sample.indices, Ys)
}

# Computes a random projection of a data matrix
random.projection <- function(X, k=2, fixed.seed=T) {
  if (fixed.seed) {
    set.seed(12345)
  }

  X <- as.matrix(X)

  # Not sure if factor is right for k > 2, but we use only k=2 for now
  factor <- sqrt(3) / sqrt(2)
  P <- matrix(sample(0:5, ncol(X)*k, replace=T), ncol=k)
  i.zeros <- P == 0
  i.ones  <- P == 1
  i.other <- P > 1
  P[i.zeros] <- factor
  P[i.ones]  <- -factor
  P[i.other] <- 0

  X %*% P
}

# Scales columns of projections so that all values are in [0, 1]
scale.Ys <- function(Ys) {
  for (j in 1:ncol(Ys)) {
    min.j <- min(Ys[, j])
    max.j <- max(Ys[, j])
    Ys[, j] <- (Ys[, j] - min.j) / (max.j - min.j)
  }

  Ys
}

# Extracts a "good" CP selection
extract.CPs <- function(Dx, k=-1) {
  if (k <= 0) {
    n <- nrow(Dx)
    k <- as.integer(sqrt(n)*3)
  }

  pam(Dx, k)$id.med
}

# Generates samples (one sample per iteration) and performs automated
# manipulation for all measures on each sample.
run.manipulation <- function(X, Dx, labels, k, ds, n.iter, output.dir) {
  n <- nrow(X)

  loginfo("Calculating all sample.indices and Ys")
  for (iter in 1:n.iter) {
    loginfo("Iteration: %02d", iter)
    # Sample dataset
    sample.indices <- sample(n, max(ncol(X), sqrt(n)*3))
    fname <- paste("sample-indices-", iter, ".tbl", sep="")
    write.table(sample.indices, file.path(output.dir, ds$name, fname), row.names=F, col.names=F)

    # Initial sample positioning
    loginfo("Calculating Ys")
    Dx.s  <- Dx[sample.indices, sample.indices]
    Ys    <- scale.Ys(cmdscale(Dx.s))
    fname <- paste("Ys-", iter, ".tbl", sep="")
    write.table(Ys, file.path(output.dir, ds$name, fname), row.names=F, col.names=F)

    # Perform manipulation
    loginfo("Running manipulation procedures")

    loginfo("Ys.m: Silhouette")
    Dx.m <- automated.silh(X[sample.indices, ], labels[sample.indices])
    Ys.silhouette <- scale.Ys(cmdscale(Dx.m))
    Ys.m <- Ys.silhouette
    fname <- paste("Ysm-silhouette-", iter, ".tbl", sep="")
    write.table(Ys.m, file.path(output.dir, ds$name, fname), row.names=F, col.names=F)

    loginfo("Ys.m: NP")
    Ys.np <- scale.Ys(Rtsne(X[sample.indices, ], perplexity=k)$Y)
    Ys.m <- Ys.np
    fname <- paste("Ysm-np-", iter, ".tbl", sep="")
    write.table(Ys.m, file.path(output.dir, ds$name, fname), row.names=F, col.names=F)

    loginfo("Ys.m: Stress")
    Ys.stress <- scale.Ys(sammon(Dx.s, Ys, tol=1e-20)$points)
    Ys.m <- Ys.stress
    fname <- paste("Ysm-stress-", iter, ".tbl", sep="")
    write.table(Ys.m, file.path(output.dir, ds$name, fname), row.names=F, col.names=F)

    loginfo("Ys.m: Precision")
    Ys.precision <- scale.Ys(nerv(Dx.s, Ys, 0.01))
    Ys.m <- Ys.precision
    fname <- paste("Ysm-precision-", iter, ".tbl", sep="")
    write.table(Ys.m, file.path(output.dir, ds$name, fname), row.names=F, col.names=F)

    loginfo("Ys.m: Recall")
    Ys.recall <- scale.Ys(nerv(Dx.s, Ys, 0.99))
    Ys.m <- Ys.recall
    fname <- paste("Ysm-recall-", iter, ".tbl", sep="")
    write.table(Ys.m, file.path(output.dir, ds$name, fname), row.names=F, col.names=F)
  }
}

run.technique <- function(X, Dx, labels, k, ds, n.iter, output.dir) {
  loginfo("Technique: %s", tech$name)
  dir.create.safe(file.path(output.dir, ds$name, tech$name))

  classes <- as.factor(labels)

  silhouette.Y <- c()
  np.Y         <- c()
  stress.Y     <- c()
  precision.Y  <- c()
  recall.Y     <- c()

  silhouette.Ym <- c()
  np.Ym         <- c()
  stress.Ym     <- c()
  precision.Ym  <- c()
  recall.Ym     <- c()

  for (iter in 1:n.iter) {
    loginfo("Iteration: %02d", iter)

    # Load sample indices...
    fname <- paste("sample-indices-", iter, ".tbl", sep="")
    sample.indices <- read.table(file.path(output.dir, ds$name, fname))$V1
    # ... and initial projection
    fname <- paste("Ys-", iter, ".tbl", sep="")
    Ys <- read.table(file.path(output.dir, ds$name, fname))

    loginfo("Calculating Y")
    Y     <- do.call(tech$fn, append(list(X=X, sample.indices=sample.indices, Ys=Ys), tech$args))
    fname <- paste("Y-", iter, ".tbl", sep="")
    write.table(Y, file.path(output.dir, ds$name, tech$name, fname), row.names=F, col.names=F)

    # Calculate distances (Y) and normalize
    loginfo("Calculating distances")
    Dy <- dist(Y)
    Dy <- Dy / mean(Dy)
    Dy <- as.matrix(Dy)

    # Calculate measures
    loginfo("Calculating measures")
    silhouette.Y <- c(silhouette.Y, mean(silhouette(Dy, classes)))
    np.Y         <- c(np.Y,         mean(NP(Dx, Dy, k)))
    stress.Y     <- c(stress.Y,     stress(Dx, Dy))
    precision.Y  <- c(precision.Y,  smoothed.pr(Dx, Dy, k)$s.precision)
    recall.Y     <- c(recall.Y,     smoothed.pr(Dx, Dy, k)$s.recall)

    # Testing manipulations
    loginfo("Projection using Ysm.silhouette")
    fname <- paste("Ysm-silhouette-", iter, ".tbl", sep="")
    Ys.m <- read.table(file.path(output.dir, ds$name, fname))
    Y.m  <- do.call(tech$fn, append(list(X=X, sample.indices=sample.indices, Ys=Ys.m), tech$args))
    fname <- paste("Ym-silhouette-", iter, ".tbl", sep="")
    write.table(Y.m, file.path(output.dir, ds$name, tech$name, fname), row.names=F, col.names=F)

    Dy <- dist(Y.m)
    Dy <- Dy / mean(Dy)
    Dy <- as.matrix(Dy)
    silhouette.Ym <- c(silhouette.Ym, mean(silhouette(Dy, classes)))


    loginfo("Projection using Ysm.np")
    fname <- paste("Ysm-np-", iter, ".tbl", sep="")
    Ys.m <- read.table(file.path(output.dir, ds$name, fname))
    Y.m <- do.call(tech$fn, append(list(X=X, sample.indices=sample.indices, Ys=Ys.m), tech$args))
    fname <- paste("Ym-np-", iter, ".tbl", sep="")
    write.table(Y.m, file.path(output.dir, ds$name, tech$name, fname), row.names=F, col.names=F)

    Dy <- dist(Y.m)
    Dy <- Dy / mean(Dy)
    Dy <- as.matrix(Dy)
    np.Ym <- c(np.Ym, mean(NP(Dx, Dy, k)))


    loginfo("Projection using Ysm.stress")
    fname <- paste("Ysm-stress-", iter, ".tbl", sep="")
    Ys.m <- read.table(file.path(output.dir, ds$name, fname))
    Y.m <- do.call(tech$fn, append(list(X=X, sample.indices=sample.indices, Ys=Ys.m), tech$args))
    fname <- paste("Ym-stress-", iter, ".tbl", sep="")
    write.table(Y.m, file.path(output.dir, ds$name, tech$name, fname), row.names=F, col.names=F)

    Dy <- dist(Y.m)
    Dy <- Dy / mean(Dy)
    Dy <- as.matrix(Dy)
    stress.Ym <- c(stress.Ym, stress(Dx, Dy))


    loginfo("Projection using Ysm.precision")
    fname <- paste("Ysm-precision-", iter, ".tbl", sep="")
    Ys.m <- read.table(file.path(output.dir, ds$name, fname))
    Y.m <- do.call(tech$fn, append(list(X=X, sample.indices=sample.indices, Ys=Ys.m), tech$args))
    fname <- paste("Ym-precision-", iter, ".tbl", sep="")
    write.table(Y.m, file.path(output.dir, ds$name, tech$name, fname), row.names=F, col.names=F)

    Dy <- dist(Y.m)
    Dy <- Dy / mean(Dy)
    Dy <- as.matrix(Dy)
    precision.Ym <- c(precision.Ym, smoothed.pr(Dx, Dy, k)$s.precision)


    loginfo("Projection using Ysm.recall")
    fname <- paste("Ysm-recall-", iter, ".tbl", sep="")
    Ys.m <- read.table(file.path(output.dir, ds$name, fname))
    Y.m <- do.call(tech$fn, append(list(X=X, sample.indices=sample.indices, Ys=Ys.m), tech$args))
    fname <- paste("Ym-recall-", iter, ".tbl", sep="")
    write.table(Y.m, file.path(output.dir, ds$name, tech$name, fname), row.names=F, col.names=F)

    Dy <- dist(Y.m)
    Dy <- Dy / mean(Dy)
    Dy <- as.matrix(Dy)
    recall.Ym <- c(recall.Ym, smoothed.pr(Dx, Dy, k)$s.recall)
  }

  write.table(silhouette.Y, file.path(output.dir, ds$name, tech$name, "silhouette-Y.tbl"), col.names=F, row.names=F)
  write.table(np.Y,         file.path(output.dir, ds$name, tech$name, "np-Y.tbl"),         col.names=F, row.names=F)
  write.table(stress.Y,     file.path(output.dir, ds$name, tech$name, "stress-Y.tbl"),     col.names=F, row.names=F)
  write.table(precision.Y,  file.path(output.dir, ds$name, tech$name, "precision-Y.tbl"),  col.names=F, row.names=F)
  write.table(recall.Y,     file.path(output.dir, ds$name, tech$name, "recall-Y.tbl"),     col.names=F, row.names=F)

  write.table(silhouette.Ym, file.path(output.dir, ds$name, tech$name, "silhouette-Ym.tbl"), col.names=F, row.names=F)
  write.table(np.Ym,         file.path(output.dir, ds$name, tech$name, "np-Ym.tbl"),         col.names=F, row.names=F)
  write.table(stress.Ym,     file.path(output.dir, ds$name, tech$name, "stress-Ym.tbl"),     col.names=F, row.names=F)
  write.table(precision.Ym,  file.path(output.dir, ds$name, tech$name, "precision-Ym.tbl"),  col.names=F, row.names=F)
  write.table(recall.Ym,     file.path(output.dir, ds$name, tech$name, "recall-Ym.tbl"),     col.names=F, row.names=F)
}

# The control points improvement experiment; n.iter sets of control points per
# dataset.
run <- function(datasets,
                techniques,
                output.dir,
                n.iter=30,
                intial.manipulation=T,
                kf=function(n) as.integer(min(sqrt(n), 0.05*n))) {
  dir.create.safe(output.dir)

  for (ds in datasets) {
    loginfo("Testing dataset: %s", ds$name)
    dir.create.safe(file.path(output.dir, ds$name))

    # Load and clean data by removing duplicates, center and scale
    X <- read.table(ds$data.file)
    if (!is.null(ds$labels.file)) {
      labels <- read.table(ds$labels.file)$V1
      labels <- labels[!duplicated(X)]
    }

    X <- unique(X)
    if (ds$scale) {
      X <- scale(X)
    }

    n <- nrow(X)
    k <- kf(n)

    # Calculate distances (X) and normalize
    loginfo("Calculating dist(X)")
    Dx <- dist(X)
    Dx <- Dx / mean(Dx)
    Dx <- as.matrix(Dx)

    # Generate samples, initial projections and all manipulations
    if (intial.manipulation) {
      run.manipulation(X, Dx, labels, k, ds, n.iter, output.dir)
    }

    # Test techniques
    for (tech in techniques) {
      run.technique(X, Dx, labels, k, ds, tech, n.iter, output.dir)
    }
  }
}

# Generates the base random CP projection and target manipulated projections for
# each measure.
run.manipulation.evo <- function(X, Dx, labels, sample.indices, k, ds, output.dir) {
  Dx.s <- Dx[sample.indices, sample.indices]

  # Initial sample positioning
  loginfo("Calculating Ys.i")
  Ys.i <- random.projection(X[sample.indices,], fixed.seed=T)
  Ys.i <- scale.Ys(Ys.i)
  write.table(Ys.i, file.path(output.dir, ds$name, "Ysi.tbl"), row.names=F, col.names=F)

  # Perform manipulation
  loginfo("Running manipulation procedures")

  loginfo("Ys.f: Silhouette")
  Dx.m <- automated.silh(X[sample.indices, ], labels[sample.indices])
  Ys.silhouette <- scale.Ys(cmdscale(Dx.m))
  Ys.m <- Ys.silhouette
  write.table(Ys.m, file.path(output.dir, ds$name, "Ysf-silhouette.tbl"), row.names=F, col.names=F)

  loginfo("Ys.f: NP")
  Ys.np <- scale.Ys(Rtsne(X[sample.indices, ], perplexity=k)$Y)
  Ys.m <- Ys.np
  write.table(Ys.m, file.path(output.dir, ds$name, "Ysf-np.tbl"), row.names=F, col.names=F)

  loginfo("Ys.f: Stress")
  Ys.stress <- scale.Ys(sammon(Dx.s, Ys.i, tol=1e-20)$points)
  Ys.m <- Ys.stress
  write.table(Ys.m, file.path(output.dir, ds$name, "Ysf-stress.tbl"), row.names=F, col.names=F)

  loginfo("Ys.f: Precision")
  Ys.precision <- scale.Ys(nerv(Dx.s, Ys.i, 0.01))
  Ys.m <- Ys.precision
  write.table(Ys.m, file.path(output.dir, ds$name, "Ysf-precision.tbl"), row.names=F, col.names=F)

  loginfo("Ys.f: Recall")
  Ys.recall <- scale.Ys(nerv(Dx.s, Ys.i, 0.99))
  Ys.m <- Ys.recall
  write.table(Ys.m, file.path(output.dir, ds$name, "Ysf-recall.tbl"), row.names=F, col.names=F)
}

# Produces Y for each interpolation step using the given technique and dataset.
run.technique.evo <- function(X, Dx, labels, k, ds, tech, n.samples, output.dir) {
  loginfo("Technique: %s", tech$name)
  dir.create.safe(file.path(output.dir, ds$name, tech$name))

  classes <- as.factor(labels)

  # Load sample indices...
  sample.indices <- read.table(file.path(output.dir, ds$name, "sample-indices.tbl"))$V1
  Dx.s <- Dx[sample.indices, sample.indices]
  # ... and initial projection
  Ys.i <- read.table(file.path(output.dir, ds$name, "Ysi.tbl"))

  alphas <- 0:(n.samples - 1)/(n.samples - 1)

  loginfo("Computing targets for measures")
  for (measure in measures) {
    if (is.null(ds$labels.file) && measure$name == "silhouette") {
      next
    }

    loginfo("Measure: %s", measure$name.pretty)
    fname <- paste("Ysf-", measure$name, ".tbl", sep="")
    Ys.f <- read.table(file.path(output.dir, ds$name, fname))

    for (iter in 1:n.samples) {
      loginfo("Calculating Y (%02d of %02d)", iter, n.samples)

      alpha <- alphas[iter]
      Ys    <- alpha * Ys.f + (1 - alpha) * Ys.i
      Y     <- do.call(tech$fn, append(list(X=X, sample.indices=sample.indices, Ys=Ys), tech$args))
      fname <- paste("Y-evo-", measure$name, "-", iter, ".tbl", sep="")
      write.table(Y, file.path(output.dir, ds$name, tech$name, fname), row.names=F, col.names=F)
    }
  }

  silhouette.Ys <- c()
  np.Ys         <- c()
  stress.Ys     <- c()
  precision.Ys  <- c()
  recall.Ys     <- c()

  silhouette.Y <- c()
  np.Y         <- c()
  stress.Y     <- c()
  precision.Y  <- c()
  recall.Y     <- c()

  loginfo("Computing measures")
  for (iter in 1:n.samples) {
    loginfo("Iteration: %02d", iter)
    alpha <- alphas[iter]

    if (!is.null(ds$labels.file)) {
      # Silhouette --------------------------------------------------------------
      Ys.f <- read.table(file.path(output.dir, ds$name, "Ysf-silhouette.tbl"))
      Ys <- alpha * Ys.f + (1 - alpha) * Ys.i

      # Calculate distances (Ys) and normalize
      loginfo("Calculating distances (Ys)")
      Dy.s <- dist(Ys)
      Dy.s <- Dy.s / mean(Dy.s)
      Dy.s <- as.matrix(Dy.s)

      fname <- paste("Y-evo-silhouette-", iter, ".tbl", sep="")
      Y <- read.table(file.path(output.dir, ds$name, tech$name, fname))

      # Calculate distances (Y) and normalize
      loginfo("Calculating distances (Y)")
      Dy <- dist(Y)
      Dy <- Dy / mean(Dy)
      Dy <- as.matrix(Dy)

      # Calculate measure
      loginfo("Calculating silhouette for Ys and Y")
      silhouette.Ys <- c(silhouette.Ys, mean(silhouette(Dy.s, classes[sample.indices])))
      silhouette.Y <- c(silhouette.Y, mean(silhouette(Dy, classes)))
    }

    # NP ----------------------------------------------------------------------
    Ys.f <- read.table(file.path(output.dir, ds$name, "Ysf-np.tbl"))
    Ys <- alpha * Ys.f + (1 - alpha) * Ys.i

    # Calculate distances (Ys) and normalize
    loginfo("Calculating distances (Ys)")
    Dy.s <- dist(Ys)
    Dy.s <- Dy.s / mean(Dy.s)
    Dy.s <- as.matrix(Dy.s)

    fname <- paste("Y-evo-np-", iter, ".tbl", sep="")
    Y <- read.table(file.path(output.dir, ds$name, tech$name, fname))

    # Calculate distances (Y) and normalize
    loginfo("Calculating distances (Y)")
    Dy <- dist(Y)
    Dy <- Dy / mean(Dy)
    Dy <- as.matrix(Dy)

    # Calculate measure
    loginfo("Calculating NP for Ys and Y")
    np.Ys <- c(np.Ys, mean(NP(Dx.s, Dy.s, k)))
    np.Y  <- c(np.Y,  mean(NP(Dx, Dy, k)))

    # Stress ------------------------------------------------------------------
    Ys.f <- read.table(file.path(output.dir, ds$name, "Ysf-stress.tbl"))
    Ys <- alpha * Ys.f + (1 - alpha) * Ys.i

    # Calculate distances (Ys) and normalize
    loginfo("Calculating distances (Ys)")
    Dy.s <- dist(Ys)
    Dy.s <- Dy.s / mean(Dy.s)
    Dy.s <- as.matrix(Dy.s)

    fname <- paste("Y-evo-stress-", iter, ".tbl", sep="")
    Y <- read.table(file.path(output.dir, ds$name, tech$name, fname))

    # Calculate distances (Y) and normalize
    loginfo("Calculating distances (Y)")
    Dy <- dist(Y)
    Dy <- Dy / mean(Dy)
    Dy <- as.matrix(Dy)

    # Calculate measure
    loginfo("Calculating stress for Ys and Y")
    stress.Ys <- c(stress.Ys, stress(Dx.s, Dy.s))
    stress.Y  <- c(stress.Y,  stress(Dx, Dy))

    # Precision ---------------------------------------------------------------
    Ys.f <- read.table(file.path(output.dir, ds$name, "Ysf-precision.tbl"))
    Ys <- alpha * Ys.f + (1 - alpha) * Ys.i

    # Calculate distances (Ys) and normalize
    loginfo("Calculating distances (Ys)")
    Dy.s <- dist(Ys)
    Dy.s <- Dy.s / mean(Dy.s)
    Dy.s <- as.matrix(Dy.s)

    fname <- paste("Y-evo-precision-", iter, ".tbl", sep="")
    Y <- read.table(file.path(output.dir, ds$name, tech$name, fname))

    # Calculate distances (Y) and normalize
    loginfo("Calculating distances (Y)")
    Dy <- dist(Y)
    Dy <- Dy / mean(Dy)
    Dy <- as.matrix(Dy)

    # Calculate measure
    loginfo("Calculating smoothed precision for Ys and Y")
    precision.Ys <- c(precision.Ys, smoothed.pr(Dx.s, Dy.s, k)$s.precision)
    precision.Y  <- c(precision.Y,  smoothed.pr(Dx, Dy, k)$s.precision)

    # Recall ------------------------------------------------------------------
    Ys.f <- read.table(file.path(output.dir, ds$name, "Ysf-recall.tbl"))
    Ys <- alpha * Ys.f + (1 - alpha) * Ys.i

    # Calculate distances (Ys) and normalize
    loginfo("Calculating distances (Ys)")
    Dy.s <- dist(Ys)
    Dy.s <- Dy.s / mean(Dy.s)
    Dy.s <- as.matrix(Dy.s)

    fname <- paste("Y-evo-recall-", iter, ".tbl", sep="")
    Y <- read.table(file.path(output.dir, ds$name, tech$name, fname))

    # Calculate distances (Y) and normalize
    loginfo("Calculating distances (Y)")
    Dy <- dist(Y)
    Dy <- Dy / mean(Dy)
    Dy <- as.matrix(Dy)

    # Calculate measure
    loginfo("Calculating smoothed recall for Ys and Y")
    recall.Ys <- c(recall.Ys, smoothed.pr(Dx.s, Dy.s, k)$s.recall)
    recall.Y  <- c(recall.Y,  smoothed.pr(Dx, Dy, k)$s.recall)
  }

  if (!is.null(ds$labels.file)) {
    write.table(silhouette.Ys, file.path(output.dir, ds$name, tech$name, "silhouette-Ys-evo.tbl"), col.names=F, row.names=F)
  }
  write.table(np.Ys,         file.path(output.dir, ds$name, tech$name, "np-Ys-evo.tbl"),         col.names=F, row.names=F)
  write.table(stress.Ys,     file.path(output.dir, ds$name, tech$name, "stress-Ys-evo.tbl"),     col.names=F, row.names=F)
  write.table(precision.Ys,  file.path(output.dir, ds$name, tech$name, "precision-Ys-evo.tbl"),  col.names=F, row.names=F)
  write.table(recall.Ys,     file.path(output.dir, ds$name, tech$name, "recall-Ys-evo.tbl"),     col.names=F, row.names=F)

  if (!is.null(ds$labels.file)) {
    write.table(silhouette.Y, file.path(output.dir, ds$name, tech$name, "silhouette-Y-evo.tbl"), col.names=F, row.names=F)
  }
  write.table(np.Y,         file.path(output.dir, ds$name, tech$name, "np-Y-evo.tbl"),         col.names=F, row.names=F)
  write.table(stress.Y,     file.path(output.dir, ds$name, tech$name, "stress-Y-evo.tbl"),     col.names=F, row.names=F)
  write.table(precision.Y,  file.path(output.dir, ds$name, tech$name, "precision-Y-evo.tbl"),  col.names=F, row.names=F)
  write.table(recall.Y,     file.path(output.dir, ds$name, tech$name, "recall-Y-evo.tbl"),     col.names=F, row.names=F)
}

# The control points improvement evolution experiment.
run.evo <- function(datasets,
                    techniques,
                    output.dir,
                    n.samples=30,
                    intial.manipulation=T,
                    kf=function(n) as.integer(min(sqrt(n), 0.05*n))) {
  dir.create.safe(output.dir)

  for (ds in datasets) {
    loginfo("Testing dataset: %s", ds$name)
    dir.create.safe(file.path(output.dir, ds$name))

    # Load and clean data by removing duplicates, center and scale
    X <- read.table(ds$data.file)
    if (!is.null(ds$labels.file)) {
      labels <- read.table(ds$labels.file)$V1
      labels <- labels[!duplicated(X)]
      classes <- as.factor(labels)
    }

    X <- unique(X)
    if (ds$scale) {
      X <- scale(X)
    }

    n <- nrow(X)
    k <- kf(n)

    # Calculate distances (X) and normalize
    loginfo("Calculating dist(X)")
    Dx <- dist(X)
    Dx <- Dx / mean(Dx)
    Dx <- as.matrix(Dx)

    loginfo("Extracting control points")
    sample.indices <- extract.CPs(Dx)
    write.table(sample.indices, file.path(output.dir, ds$name, "sample-indices.tbl"), row.names=F, col.names=F)

    # Computes each manipulation target
    run.manipulation.evo(X, Dx, labels, sample.indices, k, ds, output.dir)

    # Test techniques
    for (tech in techniques) {
      run.technique.evo(X, Dx, labels, k, ds, tech, n.samples, output.dir)
    }
  }
}


# Runs all techniques (and only the techniques) to generate all mappings from
# the original and manipulated samples.
run.Y <- function(datasets,
                  techniques,
                  output.dir,
                  n.iter=30,
                  kf=function(n) as.integer(min(sqrt(n), 0.05*n))) {
  for (ds in datasets) {
    loginfo("Testing dataset: %s", ds$name)

    # Load and clean data by removing duplicates, center and scale
    X <- read.table(ds$data.file)
    X <- unique(X)
    if (ds$scale) {
      X <- scale(X)
    }

    k <- kf(n)

    # Test techniques
    for (iter in 1:n.iter) {
      loginfo("Iteration: %d", iter)

      fname <- paste("sample-indices-", iter, ".tbl", sep="")
      sample.indices <- read.table(file.path(output.dir, ds$name, fname))$V1

      if (!is.null(ds$labels.file)) {
        fname <- paste("Ysm-silhouette-", iter, ".tbl", sep="")
        Ys.m <- read.table(file.path(output.dir, ds$name, fname))
        for (tech in techniques) {
          loginfo("Projection using Ysm.silhouette")
          Y.m <- do.call(tech$fn, append(list(X=X, sample.indices=sample.indices, Ys=Ys.m), tech$args))
          fname <- paste("Ym-silhouette-", iter, ".tbl", sep="")
          write.table(Y.m, file.path(output.dir, ds$name, tech$name, fname), row.names=F, col.names=F)
        }
      }

      fname <- paste("Ysm-np-", iter, ".tbl", sep="")
      Ys.m <- read.table(file.path(output.dir, ds$name, fname))
      for (tech in techniques) {
        loginfo("Projection using Ysm.np")
        Y.m <- do.call(tech$fn, append(list(X=X, sample.indices=sample.indices, Ys=Ys.m), tech$args))
        fname <- paste("Ym-np-", iter, ".tbl", sep="")
        write.table(Y.m, file.path(output.dir, ds$name, tech$name, fname), row.names=F, col.names=F)
      }

      fname <- paste("Ysm-stress-", iter, ".tbl", sep="")
      Ys.m <- read.table(file.path(output.dir, ds$name, fname))
      for (tech in techniques) {
        loginfo("Projection using Ysm.stress")
        Y.m <- do.call(tech$fn, append(list(X=X, sample.indices=sample.indices, Ys=Ys.m), tech$args))
        fname <- paste("Ym-stress-", iter, ".tbl", sep="")
        write.table(Y.m, file.path(output.dir, ds$name, tech$name, fname), row.names=F, col.names=F)
      }

      fname <- paste("Ysm-precision-", iter, ".tbl", sep="")
      Ys.m <- read.table(file.path(output.dir, ds$name, fname))
      for (tech in techniques) {
        loginfo("Projection using Ysm.precision")
        Y.m <- do.call(tech$fn, append(list(X=X, sample.indices=sample.indices, Ys=Ys.m), tech$args))
        fname <- paste("Ym-precision-", iter, ".tbl", sep="")
        write.table(Y.m, file.path(output.dir, ds$name, tech$name, fname), row.names=F, col.names=F)
      }

      fname <- paste("Ysm-recall-", iter, ".tbl", sep="")
      Ys.m <- read.table(file.path(output.dir, ds$name, fname))
      for (tech in techniques) {
        loginfo("Projection using Ysm.recall")
        Y.m <- do.call(tech$fn, append(list(X=X, sample.indices=sample.indices, Ys=Ys.m), tech$args))
        fname <- paste("Ym-recall-", iter, ".tbl", sep="")
        write.table(Y.m, file.path(output.dir, ds$name, tech$name, fname), row.names=F, col.names=F)
      }
    }
  }
}

# Computes confidence intervals for the difference in measures between
# manipulated and original samples.
confidence.intervals <- function(datasets, techniques, measures, output.dir, n.iter=30) {
  for (measure in measures) {
    measure.summary <- data.frame()
    for (tech in techniques) {
      for (ds in datasets) {
        if (is.null(ds$labels.file) && measure$name == "silhouette") {
          next
        }

        base.path <- file.path(output.dir, ds$name, tech$name)
        fname <- file.path(base.path, paste(measure$name, "Y.tbl", sep="-"))
        Y.measure  <- read.table(fname)$V1
        fname <- file.path(base.path, paste(measure$name, "Ym.tbl", sep="-"))
        Ym.measure <- read.table(fname)$V1
        measure.summary <- rbind(measure.summary, data.frame(tech=tech$name.pretty,
                                                            dataset=ds$name.pretty,
                                                            ci.fun(Ym.measure - Y.measure)))
      }
    }

    fname <- paste(measure$name, "-ci.tbl", sep="")
    write.table(measure.summary, file.path(output.dir, fname), col.names=T, row.names=F)
  }
}


# Experiment configuration
# Defines: datasets, techniques, output.dir
source("config.R")

args <- commandArgs(T)

# Logging setup
basicConfig()
addHandler(writeToFile,
           file=args[1],
           level="FINEST")

# The alpha and omega
run(datasets, techniques, output.dir=output.dir, initial.manipulation=F)
run.evo(datasets, techniques, output.dir=output.dir)

# Compute all confidence intervals
confidence.intervals(datasets, techniques, measures, output.dir)
