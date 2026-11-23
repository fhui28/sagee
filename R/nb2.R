#' @title Family object for negative binomial distribution.
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' A specified family object for the negative binomial distribution using the log link function and a quadratic mean-variance relationship.
#'
#' @return An object of class "family".
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>
#'
#' @export
#' @md
nb2 <- function() {
    link <- "log"
    linkfun <- function(mu)
        return(log(mu))
    linkinv <- function(eta)
        return(pmax(exp(eta), .Machine$double.eps))
    mu.eta <- function(eta)
        return(pmax(exp(eta), .Machine$double.eps))
    validmu <- function(mu)
        return(all(is.finite(mu)) && all(mu > 0))
    variance <- function(mu, phi = 0)
        return(pmax(mu+phi*mu^2, .Machine$double.eps))

    structure(list(family = "negative_binomial",
                   link = link,
                   linkfun = linkfun,
                   linkinv = linkinv,
                   mu.eta = mu.eta,
                   validmu = validmu,
                   variance = variance), class = "family")
    }
