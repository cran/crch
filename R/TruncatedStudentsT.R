TruncatedStudentsT <- function(df = numeric(), location = NULL, scale = NULL, left = NULL, right = NULL) {
  if (is.null(location)) location <- rep.int(0, length(df))
  if (is.null(scale))    scale    <- rep.int(1, length(df))
  if (is.null(left))     left     <- rep.int(-Inf, length(df))
  if (is.null(right))    right    <- rep.int(Inf, length(df))
  n <- c(length(df), length(location), length(scale), length(left), length(right))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" = all(n %in% c(1L, max(n))))
  d <- data.frame(df = df, location = location, scale = scale, left = left, right = right)
  class(d) <- c("TruncatedStudentsT", "distribution")
  d
}

mean.TruncatedStudentsT <- function(x, ...) {
  m <- ett(df = x$df, location = x$location, scale = x$scale, left = x$left, right = x$right)
  setNames(m, names(x))
}

random.TruncatedStudentsT <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rtt(n = at, df = d$df, location = d$location, scale = d$scale, left = d$left, right = d$right)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

pdf.TruncatedStudentsT <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dtt(x = at, df = d$df, location = d$location, scale = d$scale, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

log_pdf.TruncatedStudentsT <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dtt(x = at, df = d$df, location = d$location, scale = d$scale, left = d$left, right = d$right, log = TRUE)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

cdf.TruncatedStudentsT <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) ptt(q = at, df = d$df, location = d$location, scale = d$scale, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

quantile.TruncatedStudentsT <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) qtt(at, df = d$df, location = d$location, scale = d$scale, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

crps.TruncatedStudentsT <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) scoringRules::crps_tt(y = at, df = d$df, location = d$location, scale = d$scale, lower = d$left, upper = d$right)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

support.TruncatedStudentsT <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  distributions3::make_support(d$left, d$right, d, drop = drop)
}

is_discrete.TruncatedStudentsT <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

is_continuous.TruncatedStudentsT <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}

score.TruncatedStudentsT <- function(d, x, which = NULL, drop = TRUE, ...) {
  if (is.null(which)) which <- c("df", "mu", "sigma")
  which <- gsub("scale", "sigma", which, fixed = TRUE)
  which <- gsub("location", "mu", which, fixed = TRUE)
  s <- stt(x, location = d$location, scale = d$scale, df = d$df, left = d$left, right = d$right, which = which, drop = drop)
  if (!is.null(nam <- names(d))) {
    if (is.null(dim(s))) {
      names(s) <- nam
    } else {
      rownames(s) <- nam
    }
  }
  if (!is.null(dim(s))) {
    colnames(s) <- gsub("sigma", "scale", colnames(s), fixed = TRUE)
    colnames(s) <- gsub("mu", "location", colnames(s), fixed = TRUE)
  }
  return(s)
}

hessian.TruncatedStudentsT <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
  if (is.null(which)) which <- c("df", "mu:df", "sigma:df", "df:mu", "mu", "sigma:mu", "df:sigma", "mu:sigma", "sigma")
  which <- gsub("scale", "sigma", which, fixed = TRUE)
  which <- gsub("location", "mu", which, fixed = TRUE)
  h <- htt(x, location = d$location, scale = d$scale, df = d$df, left = d$left, right = d$right, which = which, drop = drop, expected = expected)
  if (!is.null(nam <- names(d))) {
    if (is.null(dim(h))) {
      names(h) <- nam
    } else {
      rownames(h) <- nam
    }
  }
  if (!is.null(dim(h))) {
    colnames(h) <- gsub("sigma", "scale", colnames(h), fixed = TRUE)
    colnames(h) <- gsub("mu", "location", colnames(h), fixed = TRUE)
  }
  return(h)
}
