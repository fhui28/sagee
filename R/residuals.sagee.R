#' @title Residuals from a fitted SAGEE model
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Extracts Pearson residuals or Dunn-Smyth (randomized quantile) residuals from a fitted `"sagee"` object, for each site-species combination.
#'
#' @param object An object of class `"sagee"`.
#' @param type Either `"pearson"` (the default) or `"dunnsmyth"` residuals to return.
#' @param seed An optional seed for the randomization used in Dunn-Smyth residuals for discrete responses, for reproducibility.
#' @param ... Not used.
#'
#' @details
#' Pearson residuals are the raw residuals divided by the square root of the fitted marginal variance. They rely only on the marginal mean and variance models, i.e., the same assumptions made by SAGEE itself.
#'
#' Dunn-Smyth residuals (Dunn and Smyth, 1996) are computed as \eqn{\Phi^{-1}(u)}, where \eqn{\Phi} denotes the standard normal distribution function, and \eqn{u = F(y)} where \eqn{F} is the fitted marginal distribution function. For discrete responses, some "jittering" is added by drawing \eqn{u} uniformly between \eqn{F(y^-)} (the CDF at the lower limit approaching from the left) and \eqn{F(y)}. If the fitted marginal distribution is correct, the Dunn-Smyth residuals are exactly standard normal, which makes them much better suited to normal probability plots than Pearson residuals for discrete data. On the other hand, they come with a big caveat as follows: since SAGEE is a generalized estimating equation approach, meaning it only assumes a model for the marginal mean and the marginal variance of the responses (along with a marginal working correlation structure across sites), then **it does not assume a full, explicit distribution for the actual response**. Dunn-Smyth residuals however require a full marginal distribution. Therefore, a compromise is made in this function by **choosing a distribution for the response that arises as the "natural" choice implied by each family's mean and mean-variance relationship**; see [mvabund::residuals.manyglm()] for an analogous assumption.
#'
#' Put simply, Dunn-Smyth residuals assess an additional assumption that SAGEE itself does not make. This means a poor-looking normal probability plot of Dunn-Smyth residuals, say, may therefore reflect an incorrect (but superfluous) distributional assumption, even when the mean and variance models, and hence the SAGEE estimates, are perfectly adequate. Pearson residuals remain the more direct check of the assumptions that SAGEE actually makes, and Dunn-Smyth residuals are best treated as a complementary, somewhat stricter, residual choice.
#'
#' @return A matrix with the same dimensions as the response matrix used in the fitted SAGEE model, containing the residuals for each site-species combination.
#'
#' @references
#' Dunn, P. K., and Smyth, G. K. (1996).
#' Randomized quantile residuals.
#' \emph{Journal of Computational and Graphical Statistics}, 5, 236--244.
#' \doi{10.1080/10618600.1996.10474708}
#'
#' @seealso
#' [SAGEE()]
#'
#' @author
#' Francis K. C. Hui \email{fhui28@gmail.com}
#'
#' @importFrom stats residuals pbinom ppois pnbinom pnorm pgamma qnorm runif
#' @importFrom statmod pinvgauss
#' @importFrom tweedie ptweedie
#'
#' @export
#' @md
residuals.sagee <- function(object,
                            type = "pearson",
                            seed = NULL,
                            ...) {

    ##---------------------------------
    #' # Checks and Pearson residuals
    ##---------------------------------
    if(!inherits(object, "sagee"))
        stop("object must be of class 'sagee'.")
    type <- match.arg(type, choices = c("pearson", "dunnsmyth"))

    if(type == "pearson")
        return(object$pearson_residuals)


    ##---------------------------------
    #' # Dunn-Smyth residuals
    ##---------------------------------
    warning("Dunn-Smyth residuals assume a full marginal distribution for the responses, whereas SAGEE only assumes a marginal mean and variance. A poor fit of these residuals may therefore reflect an incorrect distributional shape rather than a problem with the SAGEE fit itself. Please use with caution!")

    if(!is.null(seed))
        set.seed(seed)

    orig_resp <- object$residuals + object$species_fitted_values
    if(object$family$family[1] %in% c("binomial", "poisson", "negative_binomial")) {
        # Since the response is reconstructed as residuals + fitted values, then remove floating point error before evaluating discrete distribution functions
        orig_resp <- round(orig_resp)
        }

    num_sites <- nrow(object$species_fitted_values)
    num_species <- ncol(object$species_fitted_values)
    disp_param_mat <- matrix(object$dispersion_parameters, nrow = num_sites, ncol = num_species, byrow = TRUE)
    pow_param_mat <- matrix(object$power_parameters, nrow = num_sites, ncol = num_species, byrow = TRUE)


    #' ## Compute u = F(y), randomized over the jump for discrete responses
    .randomize <- function(lower, upper)
        lower + stats::runif(length(lower)) * (upper - lower)

    if(object$family$family[1] == "binomial")
        u <- .randomize(stats::pbinom(orig_resp - 1,
                                      size = 1,
                                      prob = object$species_fitted_values),
                        stats::pbinom(orig_resp,
                                      size = 1,
                                      prob = object$species_fitted_values))

    if(object$family$family[1] == "poisson")
        u <- .randomize(stats::ppois(orig_resp - 1,
                                     lambda = object$species_fitted_values),
                        stats::ppois(orig_resp, lambda = object$species_fitted_values))

    if(object$family$family[1] == "negative_binomial")
        u <- .randomize(stats::pnbinom(orig_resp - 1,
                                       size = 1 / disp_param_mat,
                                       mu = object$species_fitted_values),
                        stats::pnbinom(orig_resp,
                                       size = 1 / disp_param_mat,
                                       mu = object$species_fitted_values))

    if(object$family$family[1] == "gaussian")
        u <- stats::pnorm(orig_resp,
                          mean = object$species_fitted_values,
                          sd = sqrt(disp_param_mat))

    if(object$family$family[1] == "Gamma")
        u <- stats::pgamma(orig_resp,
                           shape = 1 / disp_param_mat,
                           rate = 1 / (disp_param_mat * object$species_fitted_values))

    if(object$family$family[1] == "inverse.gaussian")
        u <- statmod::pinvgauss(orig_resp,
                                mean = object$species_fitted_values,
                                dispersion = disp_param_mat)

    if(object$family$family[1] == "tweedie") {
        u <- matrix(NA, nrow = num_sites, ncol = num_species)
        for(j in 1:num_species) {
            F_y <- tweedie::ptweedie(orig_resp[, j],
                                     power = pow_param_mat[1, j],
                                     mu = object$species_fitted_values[, j],
                                     phi = disp_param_mat[1, j])
            is_zero <- orig_resp[, j] <= 0
            u[, j] <- F_y
            u[is_zero, j] <- .randomize(rep(0, sum(is_zero)), F_y[is_zero])
            }
        }


    #' ## Transform to the standard normal scale
    u <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps) # Guard against infinite residuals from extreme tails
    out <- matrix(stats::qnorm(u), nrow = num_sites, ncol = num_species)
    dimnames(out) <- dimnames(object$pearson_residuals)

    return(out)
    }
