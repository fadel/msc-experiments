# run.R
#
# Main experiments script.

require(logging)
require(MASS)
require(mp)
require(Rtsne)

source("measures.R")

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

# Scales columns of projections so that all values are in [0, 1]
scale.Ys <- function(Ys) {
  for (j in 1:ncol(Ys)) {
    min.j <- min(Ys[, j])
    max.j <- max(Ys[, j])
    Ys[, j] <- (Ys[, j] - min.j) / (max.j - min.j)
  }

  Ys
}

# Creates a directory at given path, optionally logging the action.
# If it already exists, the directory is not created.
dir.create.safe <- function(path, log=T) {
  if (!dir.exists(path)) {
    if (log) {
      loginfo("Creating directory: %s", path)
    }

    dir.create(path)
  }
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

relative.improvements <- function(datasets,
                                  techniques,
                                  output.dir,
                                  n.iter=30,
                                  kf=function(n) as.integer(min(sqrt(n), 0.05*n))) {
  for (ds in datasets) {
    loginfo("Dataset: %s", ds$name)

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

    loginfo("Calculating dist(X)")
    Dx <- dist(X)
    Dx <- Dx / mean(Dx)
    Dx <- as.matrix(Dx)

    # Relative improvements per measure (for all iterations)
    rcp.silhouette <- c()
    rcp.np         <- c()
    rcp.stress     <- c()
    rcp.precision  <- c()
    rcp.recall     <- c()

    for (iter in 1:n.iter) {
      loginfo("Iteration: %d", iter)

      fname <- paste("sample-indices-", iter, ".tbl", sep="")
      sample.indices <- read.table(file.path(output.dir, ds$name, fname))$V1
      Dx.s <- Dx[sample.indices, sample.indices]
      fname <- paste("Ys-", iter, ".tbl", sep="")
      Ys <- read.table(file.path(output.dir, ds$name, fname))

      loginfo("Calculating dist(Ys)")
      Dy <- dist(Ys)
      Dy <- Dy / mean(Dy)
      Dy <- as.matrix(Dy)

      loginfo("Calculating measures for Ys")
      if (!is.null(ds$labels.file)) {
        silhouette.Ys <- mean(silhouette(Dy, classes[sample.indices]))
      }
      np.Ys         <- mean(NP(Dx.s, Dy, k))
      stress.Ys     <- stress(Dx.s, Dy)
      s.pr <- smoothed.pr(Dx.s, Dy, k)
      precision.Ys  <- s.pr$s.precision
      recall.Ys     <- s.pr$s.recall

      if (!is.null(ds$labels.file)) {
        fname <- paste("Ysm-silhouette-", iter, ".tbl", sep="")
        Ys.m <- read.table(file.path(output.dir, ds$name, fname))
        Dy <- dist(Ys.m)
        Dy <- Dy / mean(Dy)
        Dy <- as.matrix(Dy)
        silhouette.Ysm <- mean(silhouette(Dy, classes[sample.indices]))
        rcp.silhouette <- c(rcp.silhouette, silhouette.Ysm / silhouette.Ys)
      }

      fname <- paste("Ysm-np-", iter, ".tbl", sep="")
      Ys.m <- read.table(file.path(output.dir, ds$name, fname))
      Dy <- dist(Ys.m)
      Dy <- Dy / mean(Dy)
      Dy <- as.matrix(Dy)
      np.Ysm <- mean(NP(Dx.s, Dy, k))
      rcp.np <- c(rcp.np, np.Ysm / np.Ys)

      fname <- paste("Ysm-stress-", iter, ".tbl", sep="")
      Ys.m <- read.table(file.path(output.dir, ds$name, fname))
      Dy <- dist(Ys.m)
      Dy <- Dy / mean(Dy)
      Dy <- as.matrix(Dy)
      stress.Ysm <- stress(Dx.s, Dy)
      rcp.stress <- c(rcp.stress, stress.Ysm / stress.Ys)

      fname <- paste("Ysm-precision-", iter, ".tbl", sep="")
      Ys.m <- read.table(file.path(output.dir, ds$name, fname))
      Dy <- dist(Ys.m)
      Dy <- Dy / mean(Dy)
      Dy <- as.matrix(Dy)
      precision.Ysm <- smoothed.pr(Dx.s, Dy, k)$s.precision
      rcp.precision <- c(rcp.precision, precision.Ysm / precision.Ys)

      fname <- paste("Ysm-recall-", iter, ".tbl", sep="")
      Ys.m <- read.table(file.path(output.dir, ds$name, fname))
      Dy <- dist(Ys.m)
      Dy <- Dy / mean(Dy)
      Dy <- as.matrix(Dy)
      recall.Ysm <- smoothed.pr(Dx.s, Dy, k)$s.recall
      rcp.recall <- c(rcp.recall, precision.Ysm / precision.Ys)
    }

    write.table(rcp.silhouette, file.path(output.dir, ds$name, "r-cp-silhouette.tbl"), col.names=F, row.names=F)
    write.table(rcp.np,         file.path(output.dir, ds$name, "r-cp-np.tbl"),         col.names=F, row.names=F)
    write.table(rcp.stress,     file.path(output.dir, ds$name, "r-cp-stress.tbl"),     col.names=F, row.names=F)
    write.table(rcp.precision,  file.path(output.dir, ds$name, "r-cp-precision.tbl"),  col.names=F, row.names=F)
    write.table(rcp.recall,     file.path(output.dir, ds$name, "r-cp-recall.tbl"),     col.names=F, row.names=F)


    for (tech in techniques) {
      loginfo("Technique: %s", tech$name)

      r.silhouette <- c()
      r.np         <- c()
      r.stress     <- c()
      r.precision  <- c()
      r.recall     <- c()

      if (!is.null(ds$labels.file)) {
        silhouette.Y <- read.table(file.path(output.dir, ds$name, tech$name, "silhouette-Y.tbl"))$V1
      }
      np.Y         <- read.table(file.path(output.dir, ds$name, tech$name, "np-Y.tbl"))$V1
      stress.Y     <- read.table(file.path(output.dir, ds$name, tech$name, "stress-Y.tbl"))$V1
      precision.Y  <- read.table(file.path(output.dir, ds$name, tech$name, "precision-Y.tbl"))$V1
      recall.Y     <- read.table(file.path(output.dir, ds$name, tech$name, "recall-Y.tbl"))$V1

      if (!is.null(ds$labels.file)) {
        silhouette.Ym <- read.table(file.path(output.dir, ds$name, tech$name, "silhouette-Ym.tbl"))$V1
      }
      np.Ym         <- read.table(file.path(output.dir, ds$name, tech$name, "np-Ym.tbl"))$V1
      stress.Ym     <- read.table(file.path(output.dir, ds$name, tech$name, "stress-Ym.tbl"))$V1
      precision.Ym  <- read.table(file.path(output.dir, ds$name, tech$name, "precision-Ym.tbl"))$V1
      recall.Ym     <- read.table(file.path(output.dir, ds$name, tech$name, "recall-Ym.tbl"))$V1

      for (iter in 1:n.iter) {
        loginfo("Iteration: %d", iter)

        if (!is.null(ds$labels.file)) {
          r.silhouette <- c(r.silhouette, silhouette.Ym[iter] / silhouette.Y[iter])
        }
        r.np         <- c(r.np,        np.Ym[iter] / np.Y[iter])
        r.stress     <- c(r.stress,    stress.Ym[iter] / stress.Y[iter])
        r.precision  <- c(r.precision, precision.Ym[iter] / precision.Y[iter])
        r.recall     <- c(r.recall,    recall.Ym[iter] / recall.Y[iter])
      }

      if (!is.null(ds$labels.file)) {
        write.table(r.silhouette, file.path(output.dir, ds$name, tech$name, "r-silhouette.tbl"), col.names=F, row.names=F)
      }
      write.table(r.np,         file.path(output.dir, ds$name, tech$name, "r-np.tbl"),         col.names=F, row.names=F)
      write.table(r.stress,     file.path(output.dir, ds$name, tech$name, "r-stress.tbl"),     col.names=F, row.names=F)
      write.table(r.precision,  file.path(output.dir, ds$name, tech$name, "r-precision.tbl"),  col.names=F, row.names=F)
      write.table(r.recall,     file.path(output.dir, ds$name, tech$name, "r-recall.tbl"),     col.names=F, row.names=F)
    }
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

# Compute relative improvements for all datasets and techniques (and samples)
relative.improvements(datasets, techniques, output.dir)
