\name{CensoredLogistic}
\alias{CensoredLogistic}
\alias{mean.CensoredLogistic}
\alias{variance.CensoredLogistic}
\alias{skewness.CensoredLogistic}
\alias{kurtosis.CensoredLogistic}
\alias{pdf.CensoredLogistic}
\alias{log_pdf.CensoredLogistic}
\alias{cdf.CensoredLogistic}
\alias{quantile.CensoredLogistic}
\alias{random.CensoredLogistic}
\alias{crps.CensoredLogistic}
\alias{support.CensoredLogistic}
\alias{is_discrete.CensoredLogistic}
\alias{is_continuous.CensoredLogistic}

\title{Create a Censored Logistic Distribution}

\description{
Class and methods for left-, right-, and interval-censored logistic
distributions using the workflow from the \pkg{distributions3} package.
}

\usage{
CensoredLogistic(location = 0, scale = 1, left = -Inf, right = Inf)
}

\arguments{
\item{location}{numeric. The location parameter of the underlying \emph{uncensored}
  logistic distribution, typically written \eqn{\mu} in textbooks.
  Can be any real number, defaults to \code{0}.}
\item{scale}{numeric. The scale parameter (standard deviation) of
  the underlying \emph{uncensored} logistic distribution,
  typically written \eqn{\sigma} in textbooks.
  Can be any positive number, defaults to \code{1}.}
\item{left}{numeric. The left censoring point. Can be any real number,
  defaults to \code{-Inf} (uncensored). If set to a finite value, the
  distribution has a point mass at \code{left} whose probability corresponds
  to the cumulative probability function of the \emph{uncensored}
  logistic distribution at this point.}
\item{right}{numeric. The right censoring point. Can be any real number,
  defaults to \code{Inf} (uncensored). If set to a finite value, the
  distribution has a point mass at \code{right} whose probability corresponds
  to 1 minus the cumulative probability function of the \emph{uncensored}
  logistic distribution at this point.}
}

\details{
The constructor function \code{CensoredLogistic} sets up a distribution
object, representing the censored logistic probability distribution by the
corresponding parameters: the latent mean \code{location} = \eqn{\mu} and
latent standard deviation \code{scale} = \eqn{\sigma} (i.e., the parameters
of the underlying \emph{uncensored} logistic variable), the \code{left} censoring
point (with \code{-Inf} corresponding to uncensored), and the
\code{right} censoring point (with \code{Inf} corresponding to uncensored). 

The censored logistic distribution has probability density function (PDF) \eqn{f(x)}:

\tabular{ll}{
\eqn{\Lambda((left - \mu)/\sigma)} \tab if \eqn{x \le left}\cr
\eqn{1 - \Lambda((right - \mu)/\sigma)} \tab if \eqn{x \ge right}\cr
\eqn{\lambda((x - \mu)/\sigma)/\sigma} \tab otherwise
}
     
where \eqn{\Lambda} and \eqn{\lambda} are the cumulative distribution function
and probability density function of the standard logistic distribution,
respectively.

All parameters can also be vectors, so that it is possible to define a vector
of censored logistic distributions with potentially different parameters.
All parameters need to have the same length or must be scalars (i.e.,
of length 1) which are then recycled to the length of the other parameters.

For the \code{CensoredLogistic} distribution objects there is a wide range
of standard methods available to the generics provided in the \pkg{distributions3}
package: \code{\link[distributions3]{pdf}} and \code{\link[distributions3]{log_pdf}}
for the (log-)density (PDF), \code{\link[distributions3]{cdf}} for the probability
from the cumulative distribution function (CDF), \code{quantile} for quantiles,
\code{\link[distributions3]{random}} for simulating random variables,
\code{\link[scoringRules]{crps}} for the continuous ranked probability score
(CRPS), and \code{\link[distributions3]{support}} for the support interval
(minimum and maximum). Internally, these methods rely on the usual d/p/q/r
functions provided for the censored logistic distributions in the \pkg{crch}
package, see \code{\link{dclogis}}, and the \code{\link[scoringRules]{crps_clogis}}
function from the \pkg{scoringRules} package.
The methods \code{\link[distributions3]{is_discrete}} and \code{\link[distributions3]{is_continuous}}
can be used to query whether the distributions are discrete on the entire support
(always \code{FALSE}) or continuous on the entire support (only \code{TRUE} if
there is no censoring, i.e., if both \code{left} and \code{right} are infinite).

See the examples below for an illustration of the workflow for the class and methods.
}

\value{
A \code{CensoredLogistic} distribution object.
}

\seealso{\code{\link{dclogis}}, \code{\link[distributions3]{Logistic}}, \code{\link{TruncatedLogistic}},
  \code{\link{CensoredNormal}}, \code{\link{CensoredStudentsT}}}

\examples{
\dontshow{ if(!requireNamespace("distributions3")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }
}
## package and random seed
library("distributions3")
set.seed(6020)

## three censored logistic distributions:
## - uncensored standard logistic
## - left-censored at zero with latent location = 1 and scale = 1
## - interval-censored in [0, 5] with latent location = 2 and scale = 1
X <- CensoredLogistic(
  location = c(   0,   1, 2),
  scale    = c(   1,   1, 1),
  left     = c(-Inf,   0, 0),
  right    = c( Inf, Inf, 5)
)
\donttest{
X
}

## compute mean of the censored distribution
mean(X)
## higher moments (variance, skewness, kurtosis) are not implemented yet

## support interval (minimum and maximum)
support(X)

## simulate random variables
random(X, 5)

## histograms of 1,000 simulated observations
x <- random(X, 1000)
hist(x[1, ], main = "uncensored")
hist(x[2, ], main = "left-censored at zero")
hist(x[3, ], main = "interval-censored in [0, 5]")

## probability density function (PDF) and log-density (or log-likelihood)
x <- c(0, 0, 1)
pdf(X, x)
pdf(X, x, log = TRUE)
log_pdf(X, x)

## cumulative distribution function (CDF)
cdf(X, x)

## quantiles
quantile(X, 0.5)

## cdf() and quantile() are inverses (except at censoring points)
cdf(X, quantile(X, 0.5))
quantile(X, cdf(X, 1))

## all methods above can either be applied elementwise or for
## all combinations of X and x, if length(X) = length(x),
## also the result can be assured to be a matrix via drop = FALSE
p <- c(0.05, 0.5, 0.95)
quantile(X, p, elementwise = FALSE)
quantile(X, p, elementwise = TRUE)
quantile(X, p, elementwise = TRUE, drop = FALSE)

## compare theoretical and empirical mean from 1,000 simulated observations
cbind(
  "theoretical" = mean(X),
  "empirical" = rowMeans(random(X, 1000))
)

## evaluate continuous ranked probability score (CRPS) using scoringRules
library("scoringRules")
crps(X, x)
}
