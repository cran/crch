## density
dclogis <- function(x, location = 0, scale = 1, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), location = as.numeric(location), scale = as.numeric(scale), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cdclogis", x, location, scale, left, right, log))
  if(is.matrix(x)) {
    rval <- matrix(rval, ncol = ncol(x), nrow = nrow(x))
    colnames(rval) <- colnames(x)
    rownames(rval) <- rownames(x)
  }
  return(rval)
}

## distribution function
pclogis <- function(q, location = 0, scale = 1, left = -Inf, right = Inf,
  lower.tail = TRUE, log.p = FALSE) {
  input <- data.frame(q = as.numeric(q), location = as.numeric(location), scale = as.numeric(scale), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cpclogis", q, location, scale, left, right, lower.tail, log.p))
  if(is.matrix(q)) {
    rval <- matrix(rval, ncol = ncol(q), nrow = nrow(q))
    colnames(rval) <- colnames(q)
    rownames(rval) <- rownames(q)
  }
  return(rval)
}

## random numbers
rclogis <- function(n, location = 0, scale = 1, left = -Inf, right = Inf) {
  rval <- rlogis(n) * scale + location
  pmax(pmin(rval, right), left)
}

## quantiles
qclogis <- function(p, location = 0, scale = 1, left = -Inf, right = Inf,
  lower.tail = TRUE, log.p = FALSE) {
  rval <- qlogis(p, lower.tail = lower.tail, log.p = log.p) * scale + location
  rval <- pmax(pmin(rval, right), left)
  if(is.matrix(p)) {
    rval <- matrix(rval, ncol = ncol(p), nrow = nrow(p))
    colnames(rval) <- colnames(p)
    rownames(rval) <- rownames(p)
  }
  return(rval)
}

## scores
sclogis <- function(x, location = 0, scale = 1, left = -Inf, right = Inf, which = NULL, drop = TRUE) {
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
    "mu"  = .Call("sclogis_mu", x, location, scale, left, right),
    "sigma" = .Call("sclogis_sigma", x, location, scale, left, right))

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

# Hessian
hclogis <- function(x, location = 0, scale = 1, left = -Inf, right = Inf, which = NULL, drop = TRUE, expected = FALSE) {
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
    "mu"    = .Call("hclogis_mu", x, location, scale, left, right),
    "sigma" = .Call("hclogis_sigma", x, location, scale, left, right),
    .Call("hclogis_musigma", x, location, scale, left, right))
  
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
eclogis <- function(location = 0, scale = 1, left = -Inf, right = Inf) {
  rmm <- (right-location)/scale
  lmm <- (left-location)/scale
  pncens <- plogis(rmm)-plogis(lmm)
  pncens*etlogis(location = location, scale = scale, left = left, right = right) + 
    plogis(lmm)*left^(is.finite(left)) + 
    plogis(rmm, lower.tail = FALSE)*right^(is.finite(right))
}
