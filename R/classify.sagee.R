#' @title Classify Species Using a Fitted SAGEE Model
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Classifies species (response variables) into archetypes based on a fitted species archetype generalized estimating equation (SAGEE) model. For each species and candidate archetype, the fitted archetypal regression coefficients are held fixed while species-specific parameters are re-estimated using [GEEfn()]. The species is then assigned to the archetype based on the same hard-classification step in
#' [SAGEE()].
#'
#' @param object A fitted object of class `"sagee"`.
#' @param new_y A new matrix of responses to classify. Rows represent sites or observational units and columns represent species or response variables.
#' @param newdata A data frame containing the variables referenced by `Xalpha_formula` and `Xbeta_formula`.
#' @param new_spatial_coords A numeric matrix with `nrow(y)` rows and two columns containing the spatial coordinates of each site. This is used when `cor_type != "independence"` was used in the fitted SAGEE `object`.
#' @param trace If `TRUE`, then print progress messages to the console.
#' @param ... Additional arguments passed to [GEEfn()].
#'
#' @details
#' Because the function has to estimate the species-specific coefficients/dispersion parameters for each species in `new_y`, then it is not true "marginal" hard-classification. Put another way, the resulting classifications are conditional on nuisance parameters estimated from the same response data, and so will be slightly optimistic.
#'
#' @return An object of class `"classify.sagee"` containing:
#' \itemize{
#'   \item \code{call} The matched function call.
#'   \item \code{species_groupings} Integer vector containing the assigned archetype for each species.
#'   \item \code{classification_criterion} A matrix containing the hard-classification criterion (effectively the Pearson residual pseudo-log-likelihood) for every species-archetype combination.
#' }
#'
#' @seealso
#' [SAGEE()], [GEEfn()]
#'
#'
#' @method classify sagee
#'
#' @importFrom stats model.matrix cov2cor
#' @importFrom fields stationary.cov
#' @export
#' @md
classify.sagee <- function(object,
                           new_y,
                           newdata,
                           new_spatial_coords = NULL,
                           trace = FALSE,
                           ...) {

    ##---------------------------------
    #' # Initial setup
    ##---------------------------------
    call <- match.call()
    if(!inherits(object, "sagee"))
        stop("object must inherit from class \"sagee\".")

    if(!is.matrix(new_y))
        new_y <- as.matrix(new_y)
    if(is.null(rownames(new_y)))
        rownames(new_y) <- paste0("site", 1:nrow(new_y))
    if(is.null(colnames(new_y)))
        colnames(new_y) <- paste0("species", 1:ncol(new_y))

    if(!is.data.frame(newdata))
        newdata <- as.data.frame(newdata)
    if(nrow(newdata) != nrow(new_y))
        stop("newdata and new_y must have the same number of rows.")

    num_sites <- nrow(new_y)
    num_species <- ncol(new_y)

    Xmat_alpha <- stats::model.matrix(object$Xalpha_formula, data = newdata)
    Xmat_beta <- stats::model.matrix(object$Xbeta_formula, data = newdata)
    alpha <- as.matrix(object$archetypal_coefficients)
    num_archetypes <- nrow(alpha)
    p_specific <- ncol(Xmat_beta)
    if(is.null(rownames(alpha)))
        rownames(alpha) <- paste0("archetype", 1:num_archetypes)


    ##---------------------------------
    #' # Construct working correlations
    ##---------------------------------
    make_spatial_correlation <- function(correlation_parameters,
                                         cor_type) {
        coordinates <- as.matrix(new_spatial_coords)
        spatial_proportion <- unname(correlation_parameters[2])

        if(cor_type == "matern") {
            spatial_component <- fields::stationary.cov(x1 = coordinates,
                                                        Covariance = "Matern",
                                                        smoothness = object$working_correlation_control$matern_smoothness,
                                                        aRange = unname(correlation_parameters[1]))

            R <- spatial_proportion * spatial_component + (1 - spatial_proportion) * diag(num_sites)
            }

        if(cor_type == "NNGP") {
            covfun <- switch(as.character(object$working_correlation_control$matern_smoothness),
                             "0.5" = "exponential_isotropic",
                             "1.5" = "matern15_isotropic",
                             "2.5" = "matern25_isotropic",
                             "matern_isotropic")

            .working_cor_from_par <- function(par,
                                              covfun,
                                              mat_smooth,
                                              locs) {
                cp <- if(covfun == "matern_isotropic") c(1, par[1], mat_smooth, 1/par[2] - 1) else c(1, par[1], 1/par[2] - 1)
                return(stats::cov2cor(getExportedValue("GpGp", covfun)(cp, locs)))
                }
            R <- .working_cor_from_par(par = correlation_parameters,
                                       covfun = covfun,
                                       mat_smooth = object$working_correlation_control$matern_smoothness,
                                       locs = coordinates)
            }

        R <- as.matrix(R)
        R <- 0.5 * (R + t(R))
        tryCatch(chol(R), error = function(e) { stop("Constructed working correlation matrix is not positive definite.") })

        return(R)
        }

    if(object$working_correlation_control$cor_type == "independence")
        working_correlations <- replicate(num_archetypes, diag(nrow = num_sites), simplify = FALSE)
    if(object$working_correlation_control$cor_type %in% c("matern", "NNGP")) {
        if(is.null(new_spatial_coords))
            stop("new_spatial_coords must be provided when the fitted SAGEE model used a non-independence working correlation.")
        if(!is.matrix(new_spatial_coords) || ncol(new_spatial_coords) != 2 || nrow(new_spatial_coords) != nrow(new_y))
            stop("new_spatial_coords must be a numeric matrix with two columns and the same number of rows as new_y.")

        if(object$working_correlation_control$common) {
            common_R <- make_spatial_correlation(object$working_correlation_params,
                                                 cor_type = object$working_correlation_control$cor_type)
            working_correlations <- replicate(num_archetypes, common_R, simplify = FALSE)
            }
        if(!object$working_correlation_control$common) {
            correlation_parameters <- object$working_correlation_params

            working_correlations <- lapply(1:num_archetypes,  function(archetype) {
                make_spatial_correlation(correlation_parameters[[archetype]],
                                         cor_type = object$working_correlation_control$cor_type)
                })
            }
        }

    names(working_correlations) <- rownames(alpha)
    workcor_chol <- lapply(working_correlations, function(R) { t(chol(R)) })
    workcor_logdet <- sapply(workcor_chol, function(Rchol) { x = 2 * sum(log(diag(Rchol))); return(as.numeric(x)) })


    ##---------------------------------
    #' # Conditional candidate fits and hard classification
    ##---------------------------------
    criterion <- matrix(Inf,
                        nrow = num_species,
                        ncol = num_archetypes,
                        dimnames = list(colnames(new_y), rownames(alpha)))

    candidate_fits <- vector("list", num_archetypes)
    names(candidate_fits) <- rownames(alpha)

    for(archetype in 1:num_archetypes) {
        if(trace)
            message(paste0("Fitting candidate model for archetype ", rownames(alpha)[archetype], " (", archetype, "/", num_archetypes, ")..."))

        #' Use the IEE estimates of the archetypal coefficients to get a starting value for the species-specific coefficients.
        start_fit <- GEEfn(y = new_y,
                           Xmat_alpha = Xmat_alpha,
                           Xmat_beta = Xmat_beta,
                           penalty_alpha = NULL,
                           workingcorrelation = diag(nrow = num_sites),
                           family = object$family,
                           start_alpha = alpha[archetype, ],
                           fix_alpha = TRUE,
                           se = FALSE,
                           ...)

        if(object$working_correlation_control$cor_type == "independence") {
            candidate_fit <- start_fit
            rm(start_fit)
            }

        if(object$working_correlation_control$cor_type != "independence") {
            candidate_fit <- GEEfn(y = new_y,
                                   Xmat_alpha = Xmat_alpha,
                                   Xmat_beta = Xmat_beta,
                                   penalty_alpha = NULL,
                                   workingcorrelation = working_correlations[[archetype]],
                                   family = object$family,
                                   start_alpha = alpha[archetype, ],
                                   start_beta = start_fit$beta,
                                   fix_alpha = TRUE,
                                   se = FALSE,
                                   ...)
            }

        candidate_pearson <- candidate_fit$residuals / sqrt(candidate_fit$variance)
        whitened_pearsonres <- forwardsolve(workcor_chol[[archetype]], candidate_pearson)
        criterion[, archetype] <- colSums(whitened_pearsonres^2) + colSums(log(candidate_fit$variance)) + workcor_logdet[archetype]

        candidate_fits[[archetype]] <- candidate_fit
        }


    classification <- max.col(-criterion, ties.method = "first")
    names(classification) <- colnames(new_y)

    #' ## Final wrap up
    out <- list(call = call,
                species_groupings = classification,
                archetype_sizes = stats::setNames(tabulate(classification, nbins = num_archetypes), rownames(alpha)),
                criterion = criterion)

    class(out) <- "classify.sagee"
    return(out)
    }


#' @title Classify Using a Fitted Model
#'
#' @description
#' A generic function for classification using a fitted model.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments passed to a class-specific method.
#'
#' @return The value returned by the applicable method.
#'
#' @export
classify <- function(object, ...) {
    UseMethod("classify")
    }
