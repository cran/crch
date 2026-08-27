## packages
library("crch")
library("distributions3")

## parameters
p <- expand.grid(x = 0:3/3, location = 0:2, scale = 1:3, df = 5, left = 0, right = 1)
x <- p$x
p <- p[-1]
tol <- .Machine$double.eps^(1/3)

## compare C-based analytic solution from crch with R-based numeric solution from distributions3
d <- do.call("CensoredStudentsT", p)
w <- c("df", "location", "scale")
all.equal(score(d, x, which = w), distributions3:::score.distribution(d, x, which = w), tolerance = tol)
w <- c("df", "location", "scale", "df:location", "df:scale", "location:scale")
all.equal(hessian(d, x, which = w), distributions3:::hessian.distribution(d, x, which = w), tolerance = tol)

d <- do.call("TruncatedStudentsT", p)
w <- c("df", "location", "scale")
all.equal(score(d, x, which = w), distributions3:::score.distribution(d, x, which = w), tolerance = tol)
w <- c("df", "location", "scale", "df:location", "df:scale", "location:scale")
all.equal(hessian(d, x, which = w), distributions3:::hessian.distribution(d, x, which = w), tolerance = tol)

p$df <- NULL
d <- do.call("CensoredLogistic", p)
w <- c("location", "scale")
all.equal(score(d, x, which = w), distributions3:::score.distribution(d, x, which = w), tolerance = tol)
w <- c("location", "scale:location", "location:scale", "scale")
all.equal(hessian(d, x, which = w), distributions3:::hessian.distribution(d, x, which = w), tolerance = tol)

d <- do.call("TruncatedLogistic", p)
w <- c("location", "scale")
all.equal(score(d, x, which = w), distributions3:::score.distribution(d, x, which = w), tolerance = tol)
w <- c("location", "scale:location", "location:scale", "scale")
all.equal(hessian(d, x, which = w), distributions3:::hessian.distribution(d, x, which = w), tolerance = tol)

names(p)[1:2] <- c("mu", "sigma")
d <- do.call("CensoredNormal", p)
w <- c("mu", "sigma")
all.equal(score(d, x, which = w), distributions3:::score.distribution(d, x, which = w), tolerance = tol)
w <- c("mu", "sigma:mu", "mu:sigma", "sigma")
all.equal(hessian(d, x, which = w), distributions3:::hessian.distribution(d, x, which = w), tolerance = tol)

d <- do.call("TruncatedNormal", p)
w <- c("mu", "sigma")
all.equal(score(d, x, which = w), distributions3:::score.distribution(d, x, which = w), tolerance = tol)
w <- c("mu", "sigma:mu", "mu:sigma", "sigma")
all.equal(hessian(d, x, which = w), distributions3:::hessian.distribution(d, x, which = w), tolerance = tol)
