CensoredStudentsT <- function(df = numeric(), location = NULL, scale = NULL, left = NULL, right = NULL) {
  if (is.null(location)) location <- rep.int(0, length(df))
  if (is.null(scale))    scale    <- rep.int(1, length(df))
  if (is.null(left))     left     <- rep.int(-Inf, length(df))
  if (is.null(right))    right    <- rep.int(Inf, length(df))
  n <- c(length(df), length(location), length(scale), length(left), length(right))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" = all(n %in% c(1L, max(n))))
  d <- data.frame(df = df, location = location, scale = scale, left = left, right = right)
  class(d) <- c("CensoredStudentsT", "distribution")
  d
}

mean.CensoredStudentsT <- function(x, ...) {
  m <- ect(df = x$df, location = x$location, scale = x$scale, left = x$left, right = x$right)
  setNames(m, names(x))
}

random.CensoredStudentsT <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rct(n = at, df = d$df, location = d$location, scale = d$scale, left = d$left, right = d$right)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

pdf.CensoredStudentsT <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dct(x = at, df = d$df, location = d$location, scale = d$scale, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

log_pdf.CensoredStudentsT <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dct(x = at, df = d$df, location = d$location, scale = d$scale, left = d$left, right = d$right, log = TRUE)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

cdf.CensoredStudentsT <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) pct(q = at, df = d$df, location = d$location, scale = d$scale, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

quantile.CensoredStudentsT <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) qct(at, df = d$df, location = d$location, scale = d$scale, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

crps.CensoredStudentsT <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) scoringRules::crps_ct(y = at, df = d$df, location = d$location, scale = d$scale, lower = d$left, upper = d$right)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

support.CensoredStudentsT <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  distributions3::make_support(d$left, d$right, d, drop = drop)
}

is_discrete.CensoredStudentsT <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

is_continuous.CensoredStudentsT <- function(d, ...) {
  setNames(!is.finite(d$left) & !is.finite(d$right), names(d))
}

score.CensoredStudentsT <- function(d, x, which = NULL, drop = TRUE, ...) {
  if (is.null(which)) which <- c("df", "mu", "sigma")
  which <- gsub("scale", "sigma", which, fixed = TRUE)
  which <- gsub("location", "mu", which, fixed = TRUE)
  s <- sct(x, location = d$location, scale = d$scale, df = d$df, left = d$left, right = d$right, which = which, drop = drop)
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

## FIXME: check numerical stability of C functions at censoring points
## hessian.CensoredStudentsT <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
##   if (is.null(which)) which <- c("df", "mu:df", "sigma:df", "df:mu", "mu", "sigma:mu", "df:sigma", "mu:sigma", "sigma")
##   which <- gsub("scale", "sigma", which, fixed = TRUE)
##   which <- gsub("location", "mu", which, fixed = TRUE)
##   h <- hct(x, location = d$location, scale = d$scale, df = d$df, left = d$left, right = d$right, which = which, drop = drop, expected = expected)
##   if (!is.null(nam <- names(d))) {
##     if (is.null(dim(h))) {
##       names(h) <- nam
##     } else {
##       rownames(h) <- nam
##     }
##   }
##   if (!is.null(dim(h))) {
##     colnames(h) <- gsub("sigma", "scale", colnames(h), fixed = TRUE)
##     colnames(h) <- gsub("mu", "location", colnames(h), fixed = TRUE)
##   }
##   return(h)
## }
hessian.CensoredStudentsT <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
  if (is.null(which)) which <- c("df", "location:df", "scale:df", "df:location", "location", "scale:location", "df:scale", "location:scale", "scale")
  NextMethod()
}
