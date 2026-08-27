## density
dct <- function(x, location = 0, scale = 1, df, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), location = as.numeric(location), scale = as.numeric(scale), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cdct", x, location, scale, df, left, right, log))
  if(is.matrix(x)) {
    rval <- matrix(rval, ncol = ncol(x), nrow = nrow(x))
    colnames(rval) <- colnames(x)
    rownames(rval) <- rownames(x)
  }
  return(rval)
}

## distribution function
pct <- function(q, location = 0, scale = 1, df, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  input <- data.frame(q = as.numeric(q), location = as.numeric(location), scale = as.numeric(scale), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cpct", q, location, scale, df, left, right, lower.tail, log.p))
  if(is.matrix(q)) {
    rval <- matrix(rval, ncol = ncol(q), nrow = nrow(q))
    colnames(rval) <- colnames(q)
    rownames(rval) <- rownames(q)
  }
  return(rval)
}

## random numbers
rct <- function(n, location = 0, scale = 1, df, left = -Inf, right = Inf) {
  rval <- rt(n, df = df) * scale + location
  pmax(pmin(rval, right), left)
}

## quantiles
qct <- function(p, location = 0, scale = 1, df, left = -Inf, right = Inf,
  lower.tail = TRUE, log.p = FALSE) {
  rval <- qt(p, df = df, lower.tail = lower.tail, log.p = log.p) * scale + location
  rval <- pmax(pmin(rval, right), left)
  if(is.matrix(p)) {
    rval <- matrix(rval, ncol = ncol(p), nrow = nrow(p))
    colnames(rval) <- colnames(p)
    rownames(rval) <- rownames(p)
  }
  return(rval)
}


## scores
sct <- function(x, location = 0, scale = 1, df, left = -Inf, right = Inf, which = NULL, drop = TRUE, eps = .Machine$double.eps^(1/3)) {
  stopifnot(
    "parameter 'scale' must always be non-negative" = all(scale >= 0),
    "parameter 'df' must always be positive" = all(df > 0)
  )
  p <- c("df", "mu", "sigma")
  if (is.null(which)) which <- p[-1L]
  which <- match.arg(tolower(which), p, several.ok = TRUE)

  ## assure that all arguments are expanded to equal length
  n <- max(length(x), length(location), length(scale), length(df), length(left), length(right))
  x <- rep_len(as.numeric(x), n)
  location <- rep_len(as.numeric(location), n)
  scale <- rep_len(as.numeric(scale), n)
  df <- rep_len(as.numeric(df), n)
  left <- rep_len(as.numeric(left), n)
  right <- rep_len(as.numeric(right), n)

  ## compute scores
  scr <- function(par) switch(par,
    "df" = (dct(x, location = location, scale = scale, df = df + eps, left = left, right = right, log = TRUE) - 
            dct(x, location = location, scale = scale, df = df - eps, left = left, right = right, log = TRUE))/(2 * eps),
    "mu"  = .Call("sct_mu", x, location, scale, df, left, right),
    "sigma" = .Call("sct_sigma", x, location, scale, df, left, right))

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
hct <- function(x, location = 0, scale = 1, df, left = -Inf, right = Inf, which = NULL, drop = TRUE, expected = FALSE, eps = .Machine$double.eps^(1/4)) {
  if (expected) stop("only the observed hessian is available")

  ## available and selected parameters/combinations and mappings for symmetries
  p <- c("df" = "df", "mu:df" = "df:mu", "sigma:df" = "df:sigma",
    "df:mu" = "df:mu", "mu" = "mu", "sigma:mu" = "mu:sigma",
    "df:sigma" = "df:sigma", "mu:sigma" = "mu:sigma", "sigma" = "sigma")
  if (is.null(which)) which <- names(p)[-c(1, 2, 3, 4, 7)]
  
  ## which combinations need to be computed?
  which <- match.arg(which, names(p), several.ok = TRUE)
  w <- unique(p[which])

  ## sanity checks
  stopifnot(
    "parameter 'scale' must always be non-negative" = all(scale >= 0),
    "parameter 'df' must always be positive" = all(df > 0)
  )
  n <- max(length(x), length(location), length(scale), length(df), length(left), length(right))
  x <- rep_len(as.numeric(x), n)
  location <- rep_len(as.numeric(location), n)
  scale <- rep_len(as.numeric(scale), n)
  df <- rep_len(as.numeric(df), n)
  left <- rep_len(as.numeric(left), n)
  right <- rep_len(as.numeric(right), n)

  ## function for computing Hessian elements (observed only)
  scr <- function(par, eps) {
    par <- strsplit(par, ":", fixed = TRUE)[[1L]]
    par <- rep_len(par, 2L)
    if (par[2L] == "mu") {
      location <- location + eps
    } else if (par[2L] == "sigma") {
      scale <- scale + eps
    } else {
      df <- df + eps
    }
    sct(x, location = location, scale = scale, df = df, left = left, right = right, which = par[1L])
  }
  ## FIXME: check numerical stability of C functions at censoring points
  hess <- function(par) switch(par,
    "mu"    = .Call("hct_mu", x, location, scale, df, left, right),
    "sigma" = .Call("hct_sigma", x, location, scale, df, left, right),
    "mu:sigma" = .Call("hct_musigma", x, location, scale, df, left, right),
    (scr(par, eps) - scr(par, -eps))/(2 * eps))
  
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
ect <- function(location = 0, scale = 1, df, left = -Inf, right = Inf) {
  rmm <- (right-location)/scale
  lmm <- (left-location)/scale
  pncens <- pt(rmm, df = df)-pt(lmm, df = df)
  pncens*ett(location = location, scale = scale, df = df, left = left, right = right) + 
    pt(lmm, df = df)*left^(is.finite(left)) + 
    pt(rmm, df = df, lower.tail = FALSE)*right^(is.finite(right))
}
