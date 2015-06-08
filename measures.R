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

silhouette <- function(Dy, labels) {
  if (any(t(Dy) != Dy)) {
    stop("Dy must be symmetric")
  }

  Dy <- as.matrix(Dy)
  n <- nrow(Dy)
  cohesion <- vector("numeric", n)
  separation <- vector("numeric", n)

  for (i in 1:n) {
    label <- labels[i]
    separation[i] <- min(Dy[i, labels != label])
    cohesion[i] <- mean(Dy[i, labels[-i] == label])
  }

  silh <- (separation - cohesion) / max(separation, cohesion)
}

d2p <- function(D, perplexity = 9, tol = 1e-5, max.tries = 50) {
  if (any(D != t(D))) {
    stop("D must be symmetric")
  }

  D <- as.matrix(D)
  P <- matrix(data=0, nrow=nrow(D), ncol=ncol(D))
  n <- nrow(D)
  beta <- rep(1, n)
  logU <- log(perplexity)

  for (i in 1:n) {
    #denom <- sum(exp(-D[i, ] / sigmas))
    #P[i, ] <- exp(-D[i, ] / sigmas) / denom

    betaMin <- -Inf
    betaMax <-  Inf
    Di <- D[i, -i]

    tries <- 0
    repeat {
      Pi <- exp(-Di * beta[i])
      sumPi <- sum(Pi)
      H <- log(sumPi) + beta[i] * sum(Di * Pi) / sumPi
      Pi <- Pi / sumPi
      Hdiff <- H - logU

      if (abs(Hdiff) < tol || tries > max.tries) {
        break
      }

      if (Hdiff > 0) {
        betaMin <- beta[i]
        beta[i] <- if (is.finite(betaMax)) {
                     (beta[i] + betaMax) / 2
                   } else {
                     beta[i] * 2
                   }
      } else {
        betaMax <- beta[i]
        beta[i] <- if (is.finite(betaMin)) {
                     (beta[i] + betaMin) / 2
                   } else {
                     beta[i] / 2
                   }
      }

      tries <- tries + 1
    }

    P[i, -i] <- Pi
  }

  list(P = P, beta = beta) # sigmas = sqrt(1 / beta)
}

d2p.beta <- function(D, beta) {
  if (any(D != t(D))) {
    stop("D must be symmetric")
  }

  D <- as.matrix(D)
  n <- nrow(D)
  P <- matrix(data=0, nrow=nrow(D), ncol=ncol(D))
  for (i in 1:n) {
    P[i, -i] <- exp(-D[i, -i] * beta[i])
    P[i, -i] <- P[i, -i] / sum(P[i, -i])
  }

  P
}

klDivergence <- function(P, Q, eps = 1e-12) {
  if (nrow(P) != ncol(P) || nrow(Q) != ncol(Q)) {
    stop("P and Q must be square")
  }
  if (nrow(P) != nrow(Q)) {
    stop("P and Q must have the same number of elements")
  }

  P[P < eps] <- eps
  Q[Q < eps] <- eps
  n <- nrow(P)
  d <- vector("numeric", n)
  for (i in 1:n) {
    d[i] <- sum(P[i, -i] * log(P[i, -i] / Q[i, -i]))
  }

  d
}
