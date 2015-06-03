stress <- function(Dx, Dy) {
  if (any(Dx != t(Dx)) || any(Dy != t(Dy))) {
    stop("Dx and Dy must be symmetric")
  }

  Dx <- as.matrix(Dx)
  Dy <- as.matrix(Dy)
  if (nrow(Dx) != nrow(Dy)) {
    stop("Dx and Dy must have the same number of elements")
  }

  n <- nrow(Dx)
  s <- vector("numeric", n)
  for (i in 1:n) {
    s[i] <- 0
    for (j in 1:n) {
      if (i == j) {
        next
      }

      s[i] = s[i] + (Dx[i, j] - Dy[i, j])^2 / Dx[i, j]
    }
    s[i] = s[i] / sum(D[i, ])
  }

  s
}

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

d2p <- function(D, sigmas) {
  if (any(D != t(D))) {
    stop("D must be symmetric")
  }

  D <- as.matrix(D)
  n <- nrow(D)
  P <- matrix(data=NA, nrow=nrow(D), ncol=ncol(D))

  for (i in 1:n) {
    denom <- sum(exp(-D[i, ] / sigmas))
    P[i, ] <- exp(-D[i, ] / sigmas) / denom
  }

  P
}

klDivergence <- function(P, Q) {
  if (nrow(P) != ncol(P) || nrow(Q) != ncol(Q)) {
    stop("P and Q must be square")
  }
  if (nrow(P) != nrow(Q)) {
    stop("P and Q must have the same number of elements")
  }

  n <- nrow(P)
  d <- vector("numeric", n)
  for (i in 1:n) {
    d[i] <- sum(P[i, ] * log(P[i, ] / Q[i, ]))
  }

  d
}
