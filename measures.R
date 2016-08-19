# Measures used as manipulation targets.
# NOTE: This file is only a library, see run.R for more details.

NP <- function(Dx, Dy, k = 9) {
  if (any(Dx != t(Dx)) || any(Dy != t(Dy))) {
    stop("Dx and Dy must be symmetric")
  }

  Dx <- as.matrix(Dx)
  Dy <- as.matrix(Dy)
  if (nrow(Dx) != nrow(Dy)) {
    stop("Dx and Dy must have the same number of elements")
  }

  n <- nrow(Dx)
  if (k >= n) {
    stop("k must be smaller than the number of elements")
  }

  preservation <- vector("numeric", n)
  for (i in 1:n) {
    nx <- order(Dx[i, ])[1 + 1:k]
    ny <- order(Dy[i, ])[1 + 1:k]
    diff <- setdiff(nx, ny)
    preservation[i] <- (k - length(diff)) / k
  }

  preservation
}

#silhouette <- function(Dy, labels) {
#  if (any(t(Dy) != Dy)) {
#    stop("Dy must be symmetric")
#  }
#
#  Dy <- as.matrix(Dy)
#  n <- nrow(Dy)
#  cohesion <- vector("numeric", n)
#  separation <- vector("numeric", n)
#
#  for (i in 1:n) {
#    label <- labels[i]
#    separation[i] <- min(Dy[i, labels != label])
#    cohesion[i] <- mean(Dy[i, labels[-i] == label])
#  }
#
#  (separation - cohesion) / max(separation, cohesion)
#}
silhouette <- function(Dy, labels) {
  n <- nrow(Dy)
  if (n != length(labels)) {
    stop("Number of labels doesn't match number of points")
  }

  A_labels <- list()
  B_labels <- list()
  unique.labels <- unique(labels)
  for (l in unique.labels) {
    A_labels[[l]] <- labels == l
    B_labels[[l]] <- labels != l
  }

  # This factor excludes self comparisons when computing cohesion
  m.factor <- n / (n-1)
  s <- vector("numeric", n)
  for (i in 1:n) {
    label_i <- labels[i]

    a <- m.factor * mean(Dy[i, A_labels[[label_i]]]) # cohesion
    b <- Inf
    for (l in unique(labels)) {
      b <- min(mean(Dy[i, B_labels[[l]]]), b) # separation
    }

    s[i] <- (b - a) / max(b, a)
  }

  s
}

#stress <- function(Dx, Dy) {
#  if (any(Dx != t(Dx)) || any(Dy != t(Dy))) {
#    stop("Dx and Dy must be symmetric")
#  }
#
#  Dx <- as.matrix(Dx)
#  Dy <- as.matrix(Dy)
#  if (nrow(Dx) != nrow(Dy)) {
#    stop("Dx and Dy must have the same number of elements")
#  }
#
#  n <- nrow(Dx)
#  s <- vector("numeric", n)
#  for (i in 1:n) {
#    s[i] <- 0
#    for (j in 1:n) {
#      if (i == j) {
#        next
#      }
#
#      s[i] = s[i] + (Dx[i, j] - Dy[i, j])^2 / Dx[i, j]
#    }
#    s[i] = s[i] / sum(D[i, ])
#  }
#
#  s
#}
stress <- function(Dx, Dy) {
  n <- nrow(Dx)

  C <- 0
  D <- 0
  for (i in 1:(n-1)) {
    for (j in (i + 1):n) {
      C <- C + Dx[i, j]
      D <- D + (Dx[i, j] - Dy[i, j])^2 / Dx[i, j]
      if (is.nan(D)) {
        loginfo("%d, %d", i, j)
        loginfo("%f", Dx[i, j])
        stop("NaN")
      }
    }
  }

  D / C
}

# NOTE: This function requires the 'klmeasure' binary from:
# http://research.cs.aalto.fi/pml/software/dredviz/
smoothed.pr <- function(Dx, Dy, k) {
  # Create SOM_PAK file for Dx
  Dx.fname <- tempfile()
  Dx.f <- file(Dx.fname, "w")
  cat(sprintf("%d\n", ncol(Dx)), file=Dx.f)
  write.table(Dx, Dx.f, col.names=F, row.names=F)
  close(Dx.f)

  # Create SOM_PAK file for Dy
  Dy.fname <- tempfile()
  Dy.f <- file(Dy.fname, "w")
  cat(sprintf("%d\n", ncol(Dy)), file=Dy.f)
  write.table(Dy, Dy.f, col.names=F, row.names=F)
  close(Dy.f)

  output <- system2("./klmeasure",
                    stdout=T,
                    args=c("--datadist", Dx.fname,
                           "--projdist", Dy.fname,
                           "--neighbors", sprintf("%d", k)))
  output <- strsplit(output[2], " | ", fixed=T)
  output <- unlist(output)

  file.remove(Dx.fname, Dy.fname)

  list(s.precision=as.double(output[1]), s.recall=as.double(output[2]))
}
