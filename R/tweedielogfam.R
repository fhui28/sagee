#' @title Family object for Tweedie distribution
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' A specified family object for the Tweedie distribution using the log link function and a power mean-variance relationship
#'
#' @return An object of class "family".
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>
#'
#' @export
#' @md
tweedielogfam <- function() {
     link <- "log"
     linkfun <- function(mu)
          return(log(mu))
     linkinv <- function(eta)
          return(pmax(exp(eta), .Machine$double.eps))
     mu.eta <- function(eta)
          return(pmax(exp(eta), .Machine$double.eps))
     validmu <- function(mu)
          return(all(is.finite(mu)) && all(mu > 0))
     variance <- function(mu, power, phi)
          return(pmax(phi * mu^power, .Machine$double.eps))

     structure(list(family = "tweedie",
                    link = link,
                    linkfun = linkfun,
                    linkinv = linkinv,
                    mu.eta = mu.eta,
                    validmu = validmu,
                    variance = variance),
               class = "family")
     }

