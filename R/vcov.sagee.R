#' @title Variance-covariance matrices for the coefficients from a fitted SAGEE model
#'
#' @description
#'
#' Returns the estimated variance-covariance matrices of the archetypal coefficients \eqn{\alpha_g} and the species-specific coefficients \eqn{\beta_i} from a fitted `"sagee"` object, as a list with one matrix per archetype.
#'
#' @param object An object of class `"sagee"`.
#' @param type The type of covariance matrix to return. Either `"model-based"` (the default), `"sandwich"`, or `"sandwich_adjusted"` (the sandwich estimator with the degrees-of-freedom adjustment of MacKinnon and White, 1985). These match the third dimension of `object$se_archetypal_coefficients`.
#' @param ... Not used.
#'
#' @details
#' Note each archetype in a SAGEE model is fitted as a separate GEE conditional on the estimated species groupings, and so the **covariance between coefficients belonging to different archetypes is zero**. This is why the variance-covariance matrices are returned as a list with the number of elements equal to the number of archetypes. The covariance matrix does not account for uncertainty in the species groupings themselves.
#'
#' For `type = "model-based"`, the matrix for archetype \eqn{g} covers the archetypal coefficients \eqn{\alpha_g}, followed by the species-specific coefficients of each species assigned to archetype \eqn{g}, in the order those species appear in the columns of `y`. Rows and columns are named by covariate for \eqn{\alpha_g} (e.g., `"temp"`), and by species and covariate for the \eqn{\beta_i} (e.g., `"species3:(Intercept)"`).
#'
#' For `type = "sandwich"` and `"sandwich_adjusted"`, the variance-covariance matrix corresponding to \eqn{\alpha_g} are returned only, as the species-specific coefficients are informed by only a single species (response) and hence lack replication for conventional sandwich inference; see [GEEfn()].
#'
#'
#' @return A list of length `num_archetypes`, named `"archetype1"`, `"archetype2"`, and so on, where each element is a symmetric covariance matrix.
#'
#'
#' @references
#' MacKinnon, J. and White, H. (1985).
#' Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties.
#' \emph{Journal of Econometrics}.
#' \doi{10.1016/0304-4076(85)90158-7}
#'
#' @seealso
#' [SAGEE()], [GEEfn()]
#'
#' @author
#' Francis K. C. Hui \email{fhui28@gmail.com}
#'
#' @importFrom stats vcov
#'
#' @export
#' @md
vcov.sagee <- function(object,
                       type = "model-based",
                       ...) {

    ##---------------------------------
    #' # Checks and set up
    ##---------------------------------
    if(!inherits(object, "sagee"))
        stop("object must be of class 'sagee'.")
    if(is.null(object$vcov_by_archetype))
        stop("object does not contain covariance information. Please refit the model with the current version of SAGEE().")
    type <- match.arg(type, choices = c("model-based", "sandwich", "sandwich_adjusted"))

    alpha <- object$archetypal_coefficients
    beta <- object$species_coefficients
    num_archetypes <- nrow(alpha)
    species_groupings <- as.integer(object$species_groupings)

    ##---------------------------------
    #' # Extract matrices and return
    ##---------------------------------
    out <- vector("list", num_archetypes)
    names(out) <- rownames(alpha)

    for(k0 in 1:num_archetypes) {
        cw_vcov <- object$vcov_by_archetype[[k0]]
        if(is.null(cw_vcov))
            next;

        if(type == "model-based") {
            #' GEEfn() orders its parameters as alpha, then the species in that archetype in increasing column order, matching y[, Ind] in SAGEE()
            species_names <- colnames(beta)[species_groupings == k0]
            coef_names <- c(colnames(alpha), paste0(rep(species_names, each = nrow(beta)), ":", rep(rownames(beta), times = length(species_names))))
            V <- cw_vcov$vcov_model
            # if(!identical(dim(V), rep(length(coef_names), 2)))
            #     stop("The stored covariance matrix for archetype ", k0, " does not match its species membership. Please refit the model.")
            }
        if(type != "model-based") {
            coef_names <- colnames(alpha)
            V <- if(type == "sandwich") cw_vcov$vcov_sandwich_alpha else cw_vcov$vcov_sandwich_alpha_adj
            }

        dimnames(V) <- list(coef_names, coef_names)
        out[[k0]] <- V
        }

    return(out)
    }
