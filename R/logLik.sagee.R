#' @title Pseudo log-likelihood from a fitted SAGEE model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Returns the Pearson residual pseudo log-likelihood of a fitted `"sagee"` object, evaluated at the estimated parameters and species groupings, along with the number of estimated parameters.
#'
#' @param object An object of class `"sagee"`.
#' @param ... Not used.
#'
#' @details
#' Since SAGEE only specifies the marginal mean, the marginal variance, and a working correlation, there is no true likelihood function. Instead, the function returns the so-called "pseudo log-likelihood", defined as
#' \deqn{
#' Q(\Psi) = -\frac{1}{2} \sum_{i} \left\{ r_i^\top R_{g_i}^{-1} r_i + \sum_{s} \log v_{si} + \log |R_{g_i}| \right\},
#' }
#' where \eqn{r_i} denotes the vector of Pearson residuals for species \eqn{i}, \eqn{v_{si}} is its fitted marginal variance at site \eqn{s}, and \eqn{R_{g_i}} is the working correlation matrix of the corresponding assigned archetype. The constant involving \eqn{2\pi} is omitted. Note the above is the same criterion used to assign species to archetypes as part of the SAGEE model fitting process.
#'
#' The `"df"` attribute counts the number of estimated parameters as follows:
#' \itemize{
#'   \item The number of archetypal and species-specific coefficients, i.e., `length(object$archetypal_coefficients) + length(object$species_coefficients)`.
#'   \item Any estimated species-specific dispersion and power parameters.
#'   \item The number of working correlation parameters, i.e., the spatial range and spatial proportion when `cor_type = "matern"` or `cor_type = "NNGP"`, either once if `common = TRUE` or once per archetype if `common = FALSE`. No parameters are counted for `cor_type = "independence"` or `cor_type = "supplied"`.
#'   }
#'
#' @note
#' Because the result is an object of class `"logLik"` with a `"df"` attribute, then functions such as [stats::AIC()] can be applied to a fitted SAGEE model. However, these **should be taken with a large grain of salt** as standard information criteria are not formally justified for a pseudo log-likelihood.
#'
#' @return An object of class `"logLik"`.
#'
#' @seealso
#' [SAGEE()]
#'
#' @author
#' Francis K. C. Hui \email{fhui28@gmail.com}
#'
#' @importFrom stats logLik
#'
#' @export
#' @md
logLik.sagee <- function(object, ...) {
    if(!inherits(object, "sagee"))
        stop("object must be of class 'sagee'.")

    num_species <- ncol(object$species_coefficients)
    num_dispersion <- switch(object$family$family[1],
                             gaussian = num_species,
                             Gamma = num_species,
                             inverse.gaussian = num_species,
                             negative_binomial = num_species,
                             tweedie = 2 * num_species,
                             0)

    num_correlation <- sum(!is.na(unlist(object$working_correlation_params)))
    df <- length(object$archetypal_coefficients) + length(object$species_coefficients) + num_dispersion + num_correlation

    out <- object$pseudo_logLik
    attr(out, "df") <- df
    class(out) <- "logLik"
    return(out)
    }
