## density
dtnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cdtnorm", x, mean, sd, left, right, log))
  if(is.matrix(x)) {
    rval <- matrix(rval, ncol = ncol(x), nrow = nrow(x))
    colnames(rval) <- colnames(x)
    rownames(rval) <- rownames(x)
  }
  return(rval)
}


## distribution function
ptnorm <- function(q, mean = 0, sd = 1, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  input <- data.frame(q = as.numeric(q), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cptnorm", q, mean, sd, left, right, lower.tail, log.p))
  if(is.matrix(q)) {
    rval <- matrix(rval, ncol = ncol(q), nrow = nrow(q))
    colnames(rval) <- colnames(q)
    rownames(rval) <- rownames(q)
  }
  return(rval)
}

## quantiles
qtnorm <- function(p, mean = 0, sd = 1, left = -Inf, right = Inf,
  lower.tail = TRUE, log.p = FALSE) {
  if(log.p) p <- exp(p) 
  lower <- if(lower.tail) left else right
  upper <- if(lower.tail) right else left
  p <- pnorm((lower-mean)/sd, lower.tail = lower.tail) * (1 - p) + 
    p*pnorm((upper - mean)/sd, lower.tail = lower.tail)
  rval <- qnorm(p, lower.tail = lower.tail)*sd + mean
  if(is.matrix(p)) {
    rval <- matrix(rval, ncol = ncol(p), nrow = nrow(p))
    colnames(rval) <- colnames(p)
    rownames(rval) <- rownames(p)
  }
  return(rval)
}

## random numbers
rtnorm <- function(n, mean = 0, sd = 1, left = -Inf, right = Inf) {
  qtnorm(runif(n), mean, sd, left = left, right = right)
}



## scores
stnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, which = NULL, drop = TRUE) {
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
    "mu"  = .Call("stnorm_mu", x, mean, sd, left, right),
    "sigma" = .Call("stnorm_sigma", x, mean, sd, left, right))

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
htnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, which = NULL, drop = TRUE, expected = FALSE) {
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
    "mu"    = .Call("htnorm_mu", x, mean, sd, left, right),
    "sigma" = .Call("htnorm_sigma", x, mean, sd, left, right),
    .Call("htnorm_musigma", x, mean, sd, left, right))
  
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

.erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
.erfc <- function(x) 2 * pnorm(x * sqrt(2), lower.tail = FALSE)
.erfcx <- function(x) 2 * pnorm(x * sqrt(2), lower.tail = FALSE) * exp(x^2)
.F1 <- function(x, y) {
    delta <- exp(x^2 - y^2)
    fx <- is.finite(x)
    fy <- is.finite(y)
    sx <- sign(x)
    sy <- sign(y)
    ifelse(fx & !fy, sy / .erfcx(sy * x),
    ifelse(!fx & fy, sx / .erfcx(sx * y),
    ifelse(abs(x) > y & y >= 0, (exp(-y^2) - exp(-x^2)) / (.erf(x) - .erf(y)),
    ifelse(x < 0 & y < 0, (1 - delta) / (delta * .erfcx(-y) - .erfcx(-x)),
    ifelse(x > 0 & y > 0, (1 - delta) / (.erfcx(x) - delta * .erfcx(y)),
    (1 - delta) * exp(-x^2) / (.erf(y) - .erf(x)))))))
}
.F2 <- function(x, y) {
    delta <- exp(x^2 - y^2)
    fx <- is.finite(x)
    fy <- is.finite(y)
    sx <- sign(x)
    sy <- sign(y)
    ifelse(fx & !fy, sy * x / .erfcx(sy * x),
    ifelse(!fx & fy, sx * y / .erfcx(sx * y),
    ifelse(abs(x) > y & y >= 0, (y * exp(-y^2) - x * exp(-x^2)) / (.erf(x) - .erf(y)),
    ifelse(x < 0 & y < 0, (x - y * delta) / (delta * .erfcx(-y) - .erfcx(-x)),
    ifelse(x > 0 & y > 0, (x - y * delta) / (.erfcx(x) - delta * .erfcx(y)),
    (x - y * delta) * exp(-x^2) / (.erf(y) - .erf(x)))))))
}


## Using the expressions in
## https://github.com/cossio/TruncatedNormal.jl/blob/fc904152f2da11a257e3ccdd3e49ef118b81d437/notes/normal.pdf
## to avoid catastrophic cancellation

etnorm <- function (mean = 0, sd = 1, left = -Inf, right = Inf) {
    rmm <- (right - mean) / sd / sqrt(2)
    lmm <- (left - mean) / sd / sqrt(2)
    ifelse(rmm == Inf & lmm == -Inf, mean,
           mean + sqrt(2 / pi) * .F1(lmm, rmm) * sd)
}

sdtnorm <- function (mean = 0, sd = 1, left = -Inf, right = Inf) {
    rmm <- (right - mean) / sd / sqrt(2)
    lmm <- (left - mean) / sd / sqrt(2)
    ifelse(rmm == Inf & lmm == -Inf, sd,
           sd * sqrt(1 + 2 / sqrt(pi) * .F2(lmm, rmm) - 2 / pi * (.F1(lmm, rmm))^2))
}

