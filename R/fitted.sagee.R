#' @title Fitted values from a fitted SAGEE model
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Returns the fitted marginal means for each species from a fitted `"sagee"` object.
#'
#' @param object An object of class `"sagee"`.
#' @param ... Not used.
#'
#' @details
#' For species \eqn{i} assigned to archetype \eqn{g_i}, the fitted value is \eqn{g^{-1}(X_{\alpha}\widehat{\alpha}_{g_i} + X_{\beta}\widehat{\beta}_i)}, where \eqn{\alpha_{g_i}} denotes the vector of coefficients for the assigned species archetype, and \eqn{\beta_i} is the coefficient vector specific to species \eqn{i}. The species-specific component should generally contain the intercept and any covariates whose effects are not intended to characterize the species archetypes. Finally, \eqn{g^{-1}} denotes the inverse link function of the assumed `family`.
#'
#' @return A matrix with the same dimensions as the response matrix `y` inputted into the fitted SAGEE model, where column \eqn{i} contains the fitted marginal means for species \eqn{i}.
#'
#' @seealso
#' [SAGEE()], [summary.sagee()]
#'
#' @author
#' Francis K. C. Hui \email{fhui28@gmail.com}
#'
#' @importFrom stats fitted
#'
#' @export
#' @md
fitted.sagee <- function(object, ...) {
    if(!inherits(object, "sagee"))
        stop("object must be of class 'sagee'.")

    return(object$species_fitted_values)
    }
