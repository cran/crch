## density
dcnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cdcnorm", x, mean, sd, left, right, log))
  if(is.matrix(x)) {
    rval <- matrix(rval, ncol = ncol(x), nrow = nrow(x))
    colnames(rval) <- colnames(x)
    rownames(rval) <- rownames(x)
  }
  return(rval)
}

## distribution function
pcnorm <- function(q, mean = 0, sd = 1, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  input <- data.frame(q = as.numeric(q), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cpcnorm", q, mean, sd, left, right, lower.tail, log.p))
  if(is.matrix(q)) {
    rval <- matrix(rval, ncol = ncol(q), nrow = nrow(q))
    colnames(rval) <- colnames(q)
    rownames(rval) <- rownames(q)
  }
  return(rval)
}

## random numbers
rcnorm <- function(n, mean = 0, sd = 1, left = -Inf, right = Inf) {
  rval <- rnorm(n) * sd + mean
  pmax(pmin(rval, right), left)
}

## quantiles
qcnorm <- function(p, mean = 0, sd = 1, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  rval <- qnorm(p, lower.tail = lower.tail, log.p = log.p) * sd + mean
  rval <- pmax(pmin(rval, right), left)
  if(is.matrix(p)) {
    rval <- matrix(rval, ncol = ncol(p), nrow = nrow(p))
    colnames(rval) <- colnames(p)
    rownames(rval) <- rownames(p)
  }
  return(rval)
}

## scores
scnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, which = NULL, drop = TRUE) {
  stopifnot(
    "parameter 'sd' must always be non-negative" = all(sd >= 0)
  )
  p <- c("mu", "sigma")
  if (is.null(which)) which <- p
  which <- match.arg(tolower(which), p, several.ok = TRUE)

  ## assure that all arguments are expanded to equal length
  n <- max(length(x), length(mean), length(sd), length(left), length(right))
  x <- rep_len(as.numeric(x), n)
  mean <- rep_len(as.numeric(mean), n)
  sd <- rep_len(as.numeric(sd), n)
  left <- rep_len(as.numeric(left), n)
  right <- rep_len(as.numeric(right), n)

  ## compute scores
  scr <- function(par) switch(par,
    "mu"  = .Call("scnorm_mu", x, mean, sd, left, right),
    "sigma" = .Call("scnorm_sigma", x, mean, sd, left, right))

  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    s <- scr(which)
  } else {
    s <- lapply(which, scr)
    s <- do.call("cbind", s)
    colnames(s) <- which
  }
  return(s)
}

## Hessian
hcnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, which = NULL, drop = TRUE, expected = FALSE) {
  if (expected) stop("only the observed hessian is available")

  ## available and selected parameters/combinations and mappings for symmetries
  p <- c("mu" = "mu", "sigma:mu" = "mu:sigma", "mu:sigma" = "mu:sigma", "sigma" = "sigma")
  if (is.null(which)) which <- names(p)
  
  ## which combinations need to be computed?
  which <- match.arg(which, names(p), several.ok = TRUE)
  w <- unique(p[which])

  ## sanity checks
  stopifnot(
    "parameter 'sd' must always be non-negative" = all(sd >= 0)
  )
  n <- max(length(x), length(mean), length(sd), length(left), length(right))
  x <- rep_len(as.numeric(x), n)
  mean <- rep_len(as.numeric(mean), n)
  sd <- rep_len(as.numeric(sd), n)
  left <- rep_len(as.numeric(left), n)
  right <- rep_len(as.numeric(right), n)

  ## function for computing Hessian elements (observed only)
  hess <- function(par) switch(par,
    "mu"    = .Call("hcnorm_mu", x, mean, sd, left, right),
    "sigma" = .Call("hcnorm_sigma", x, mean, sd, left, right),
    .Call("hcnorm_musigma", x, mean, sd, left, right))
  
  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    h <- hess(w)
  } else {
    h <- lapply(w, hess)
    h <- do.call("cbind", h)
    colnames(h) <- w
    if (!identical(w, which)) h <- h[, p[which], drop = FALSE]
    colnames(h) <- which
  }
  return(h)
}


## Expectation
ecnorm <- function(mean = 0, sd = 1, left = -Inf, right = Inf) {
  rmm <- (right-mean)/sd
  lmm <- (left-mean)/sd
  pncens <- pnorm(rmm)-pnorm(lmm)
  pncens*etnorm(mean = mean, sd = sd, left = left, right = right) + 
    pnorm(lmm)*left^(is.finite(left)) + 
    pnorm(rmm, lower.tail = FALSE)*right^(is.finite(right))
}


## Standard deviation
sdcnorm <- function(mean = 0, sd = 1, left = -Inf, right = Inf) {
    rmm <- (right - mean) / sd
    lmm <- (left - mean) / sd
    pl <- pnorm(lmm)
    pr <- pnorm(rmm, lower.tail = FALSE)
    pm <- 1 - pr - pl
    E <- etnorm(mean = mean, sd = sd, left = left, right = right)
    V <- sdtnorm(mean = mean, sd = sd, left = left, right = right)^2
    left <- left^is.finite(left)
    right <- right^is.finite(right)
    rval <- left^2 * pl * (1 - pl) + right^2 * pr * (1 - pr) - 2 * left * right * pl * pr +
        E^2 * pm * (1 - pm) - 2 * E * pm * (left * pl + right * pr) + V * pm
    sqrt(rval)
}
