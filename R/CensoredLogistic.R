CensoredLogistic <- function(location = 0, scale = 1, left = -Inf, right = Inf) {
  n <- c(length(location), length(scale), length(left), length(right))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" = all(n %in% c(1L, max(n))))
  d <- data.frame(location = location, scale = scale, left = left, right = right)
  class(d) <- c("CensoredLogistic", "distribution")
  d
}

mean.CensoredLogistic <- function(x, ...) {
  m <- eclogis(location = x$location, scale = x$scale, left = x$left, right = x$right)
  setNames(m, names(x))
}

random.CensoredLogistic <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rclogis(n = at, location = d$location, scale = d$scale, left = d$left, right = d$right)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

pdf.CensoredLogistic <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dclogis(x = at, location = d$location, scale = d$scale, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

log_pdf.CensoredLogistic <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dclogis(x = at, location = d$location, scale = d$scale, left = d$left, right = d$right, log = TRUE)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

cdf.CensoredLogistic <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) pclogis(q = at, location = d$location, scale = d$scale, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

quantile.CensoredLogistic <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) qclogis(at, location = d$location, scale = d$scale, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

crps.CensoredLogistic <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) scoringRules::crps_clogis(y = at, location = d$location, scale = d$scale, lower = d$left, upper = d$right)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

support.CensoredLogistic <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  distributions3::make_support(d$left, d$right, d, drop = drop)
}

is_discrete.CensoredLogistic <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

is_continuous.CensoredLogistic <- function(d, ...) {
  setNames(!is.finite(d$left) & !is.finite(d$right), names(d))
}

score.CensoredLogistic <- function(d, x, which = NULL, drop = TRUE, ...) {
  if (is.null(which)) which <- c("mu", "sigma")
  which <- gsub("scale", "sigma", which, fixed = TRUE)
  which <- gsub("location", "mu", which, fixed = TRUE)
  s <- sclogis(x, location = d$location, scale = d$scale, left = d$left, right = d$right, which = which, drop = drop)
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

hessian.CensoredLogistic <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
  if (is.null(which)) which <- c("mu", "sigma:mu", "mu:sigma", "sigma")
  which <- gsub("scale", "sigma", which, fixed = TRUE)
  which <- gsub("location", "mu", which, fixed = TRUE)
  h <- hclogis(x, location = d$location, scale = d$scale, left = d$left, right = d$right, which = which, drop = drop, expected = expected)
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
