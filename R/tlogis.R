## density
dtlogis <- function(x, location = 0, scale = 1, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), location = as.numeric(location), scale = as.numeric(scale), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cdtlogis", x, location, scale, left, right, log))
  if(is.matrix(x)) {
    rval <- matrix(rval, ncol = ncol(x), nrow = nrow(x))
    colnames(rval) <- colnames(x)
    rownames(rval) <- rownames(x)
  }
  return(rval)
}


## distribution function
ptlogis <- function(q, location = 0, scale = 1, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  input <- data.frame(q = as.numeric(q), location = as.numeric(location), scale = as.numeric(scale), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cptlogis", q, location, scale, left, right, lower.tail, log.p))
  if(is.matrix(q)) {
    rval <- matrix(rval, ncol = ncol(q), nrow = nrow(q))
    colnames(rval) <- colnames(q)
    rownames(rval) <- rownames(q)
  }
  return(rval)
}

## quantiles
qtlogis <- function(p, location = 0, scale = 1, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  if(log.p) p <- exp(p) 
  lower <- if(lower.tail) left else right
  upper <- if(lower.tail) right else left
  p <- plogis((lower-location)/scale, lower.tail = lower.tail) * (1 - p) + 
    p*plogis((upper - location)/scale, lower.tail = lower.tail)
  rval <- qlogis(p, lower.tail = lower.tail)*scale + location
  if(is.matrix(p)) {
    rval <- matrix(rval, ncol = ncol(p), nrow = nrow(p))
    colnames(rval) <- colnames(p)
    rownames(rval) <- rownames(p)
  }
  return(rval)
}

## random numbers
rtlogis <- function(n, location = 0, scale = 1, left = -Inf, right = Inf) {
  qtlogis(runif(n), location, scale, left = left, right = right)
}

## scores
stlogis <- function(x, location = 0, scale = 1, left = -Inf, right = Inf, which = NULL, drop = TRUE) {
  stopifnot(
    "parameter 'scale' must always be non-negative" = all(scale >= 0)
  )
  p <- c("mu", "sigma")
  if (is.null(which)) which <- p
  which <- match.arg(tolower(which), p, several.ok = TRUE)

  ## assure that all arguments are expanded to equal length
  n <- max(length(x), length(location), length(scale), length(left), length(right))
  x <- rep_len(as.numeric(x), n)
  location <- rep_len(as.numeric(location), n)
  scale <- rep_len(as.numeric(scale), n)
  left <- rep_len(as.numeric(left), n)
  right <- rep_len(as.numeric(right), n)

  ## compute scores
  scr <- function(par) switch(par,
    "mu"  = .Call("stlogis_mu", x, location, scale, left, right),
    "sigma" = .Call("stlogis_sigma", x, location, scale, left, right))

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
htlogis <- function(x, location = 0, scale = 1, left = -Inf, right = Inf, which = NULL, drop = TRUE, expected = FALSE) {
  if (expected) stop("only the observed hessian is available")

  ## available and selected parameters/combinations and mappings for symmetries
  p <- c("mu" = "mu", "sigma:mu" = "mu:sigma", "mu:sigma" = "mu:sigma", "sigma" = "sigma")
  if (is.null(which)) which <- names(p)
  
  ## which combinations need to be computed?
  which <- match.arg(which, names(p), several.ok = TRUE)
  w <- unique(p[which])

  ## sanity checks
  stopifnot(
    "parameter 'scale' must always be non-negative" = all(scale >= 0)
  )
  n <- max(length(x), length(location), length(scale), length(left), length(right))
  x <- rep_len(as.numeric(x), n)
  location <- rep_len(as.numeric(location), n)
  scale <- rep_len(as.numeric(scale), n)
  left <- rep_len(as.numeric(left), n)
  right <- rep_len(as.numeric(right), n)

  ## function for computing Hessian elements (observed only)
  hess <- function(par) switch(par,
    "mu"    = .Call("htlogis_mu", x, location, scale, left, right),
    "sigma" = .Call("htlogis_sigma", x, location, scale, left, right),
    .Call("htlogis_musigma", x, location, scale, left, right))
  
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
etlogis <- function(location = 0, scale = 1, left = -Inf, right = Inf) {
  ## scaled truncation points
  rmm <- (right - location)/scale
  lmm <- (left - location)/scale
  
  ## non-truncated probability
  pncens <- plogis(rmm) - plogis(lmm)

  ## effect on right and left truncation point
  rmm <- rmm * plogis(rmm) - log(1 + exp(rmm))
  rmm[!is.finite(rmm) | is.nan(rmm)] <- 0
  lmm <- lmm * plogis(lmm) - log(1 + exp(lmm))
  lmm[!is.finite(lmm) | is.nan(lmm)] <- 0

  ## mean of truncated variable
  rval <- location + scale * (rmm - lmm) / pncens
  rval
}
