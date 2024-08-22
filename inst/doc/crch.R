## ----preliminaries, include = FALSE-------------------------------------------
knitr::opts_chunk$set(
  engine = "R",
  collapse = TRUE,
  comment = "##",
  message = FALSE,
  warning = FALSE,
  echo = TRUE
)
options(width = 70, prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE, digits = 5)
library("crch")


## -----------------------------------------------------------------------------
#| eval: false
## crch(formula, data, subset, na.action, weights, offset, link.scale = "log",
##   dist = "gaussian", df = NULL, left = -Inf, right = Inf, truncated = FALSE,
##   type = "ml", control = crch.control(...),
##   model = TRUE, x = FALSE, y = FALSE, ...)


## -----------------------------------------------------------------------------
library("crch")
data("RainIbk", package = "crch")


## -----------------------------------------------------------------------------
RainIbk <- sqrt(RainIbk)
RainIbk$ensmean <- apply(RainIbk[,grep('^rainfc',names(RainIbk))], 1, mean)
RainIbk$enssd <- apply(RainIbk[,grep('^rainfc',names(RainIbk))], 1, sd)
RainIbk <- subset(RainIbk, enssd > 0)


## -----------------------------------------------------------------------------
#| eval: false
## plot(rain ~ ensmean, data = RainIbk, pch = 19, col = gray(0, alpha = 0.2))
## abline(0, 1, col = 2)


## -----------------------------------------------------------------------------
CRCH <- crch(rain ~ ensmean | log(enssd), data = RainIbk, left = 0, 
  dist = "logistic")
summary(CRCH)


## -----------------------------------------------------------------------------
#| eval: false
## abline(coef(CRCH)[1:2], col = 4)


## ----scatter------------------------------------------------------------------
#| fig-width: 9
#| fig-height: 5
#| out-width: 100%
#| echo: false
#| label: fig-scatter
#| fig-cap: "Square rooted precipitation amount against ensemble mean forecasts. A line with intercept 0 and slope 1 is shown in red and the censored regression fit in blue."
plot(rain~ensmean, data = RainIbk, pch = 19, col = gray(0, alpha = 0.2))
abline(0,1, col = 2, lwd = 2)
abline(coef(CRCH)[1:2], col = 4, lwd = 2)


## -----------------------------------------------------------------------------
CR <- crch(rain ~ ensmean, data = RainIbk, left = 0, dist = "logistic")
cbind(AIC(CR, CRCH), BIC = BIC(CR, CRCH)[,2])


## -----------------------------------------------------------------------------
CRCHgau <- crch(rain ~ ensmean | log(enssd), data = RainIbk, left = 0, 
  dist = "gaussian")
CRCHstud <- crch(rain ~ ensmean | log(enssd), data = RainIbk, left = 0, 
  dist = "student")
AIC(CRCH, CRCHgau, CRCHstud)


## ----dist---------------------------------------------------------------------
#| fig-width: 9
#| fig-height: 5
#| out-width: 100%
#| echo: false
#| label: fig-dist
#| fig-cap: "Probability density functions of a student-t distribution with 9.56 degrees of freedom, a logistic, and a normal distribution. The densities of the logistic and normal distribution are scaled to facilitate comparison."
a <- seq(-5, 5, 0.01)
df <- exp(tail(coef(CRCHstud), 1))
#par(mfrow = c(1,2))
plot(a, dt(a, df), type = "l", xlab = "", ylab = "density", main = "Probability density function", lwd = 2)
lines(a, dlogis(a*dt(0, df)/dlogis(0))*dt(0, df)/dlogis(0), col = 2, lwd = 2, lty = 2)
lines(a, dnorm(a*dt(0, df)/dnorm(0))*dt(0, df)/dnorm(0), col = 4, lwd = 1, lty = 1)

#plot(a, pt(a, df), type = "l", xlab = "", ylab = "density", main = "Cumulative distribution function", lwd = 2)
#lines(a, plogis(a*dt(0, df)/dlogis(0)), col = 2, lwd = 2, lty = 2)
legend("topright", lwd = c(2,2,1), lty = c(1,2,1), col = c(1,2,4), c("student-t", "scaled logistic", "scaled normal"), bty = "n")


## -----------------------------------------------------------------------------
library("glmx")
BIN <- hetglm(I(rain > 0) ~ ensmean | log(enssd), data = RainIbk,
  family = binomial(link = "logit"))
TRCH <- crch(rain~ensmean | log(enssd), data = RainIbk, subset = rain > 0, 
  left = 0, dist = "logistic", truncated = TRUE)


## -----------------------------------------------------------------------------
cbind("CRCH" = c(coef(CRCH, "location")/exp(coef(CRCH, "scale"))[1], 
    coef(CRCH, "scale")[2]), 
  "BIN" = coef(BIN), 
  "TRCH" = c(coef(TRCH, "location")/exp(coef(TRCH, "scale"))[1], 
    coef(TRCH, "scale")[2]))


## -----------------------------------------------------------------------------
loglik <- c("Censored" = logLik(CRCH), "Two-Part" = logLik(BIN) + logLik(TRCH))
df <- c(4, 7)
aic <- -2 * loglik + 2 * df
bic <- -2 * loglik + log(nrow(RainIbk)) * df
cbind(df, AIC = aic, BIC = bic)


## -----------------------------------------------------------------------------
newdata <- data.frame(ensmean = 1.8, enssd = 0.9)
predict(CRCH, newdata, type = "quantile", at = 0.5)^2


## -----------------------------------------------------------------------------
p <- predict(BIN, newdata)
predict(TRCH, newdata, type = "quantile", at = (p - 0.5)/p)^2


## -----------------------------------------------------------------------------
mu <- predict(CRCH, newdata, type = "location")
sigma <- predict(CRCH, newdata, type = "scale")
pclogis(sqrt(5), mu, sigma, lower.tail = FALSE, left = 0)

mu <- predict(TRCH, newdata, type = "location")
sigma <- predict(TRCH, newdata, type = "scale")
p * ptlogis(sqrt(5), mu, sigma, lower.tail = FALSE, left = 0)

