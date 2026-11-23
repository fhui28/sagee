#' @title Compute a working correlation structure from a fitted SAGEE model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Takes a fitted `"sagee"` object and estimates a working correlation matrix (or a list of matrices, one per archetype) from its Pearson residuals, given the newly supplied arguments in `working_correlation_control`. The main use case of this function is when the original SAGEE fit assumed independence, as the resulting spatial working correlation matrix/list of matrices estimated from this function can be passed back into [SAGEE()] via `cor_type = "supplied"`.
#'
#' @param object An object of class `"sagee"`.
#' @param working_correlation_control Working correlation structure among the rows of `y`, within each species, to estimate. A named list with the following components:
#' \itemize{
#'  \item \code{cor_type} A character string specifying the type of working correlation structure. Currently available options are `"independence"`, `"matern"`, and `"NNGP"`.
#'  \item \code{matern_smoothness} A positive numeric value specifying the smoothness parameter of the Matern correlation function. This is used when `cor_type == "matern"` or `cor_type == "NNGP"`. Note for `cor_type = "NNGP"`, setting `matern_smoothness` to be 0.5, 1.5, or 2.5 generally gives faster computation.
#'  \item \code{num_neighbors} A positive integer specifying the number of nearest neighbors to use when constructing the NNGP working correlation structure. This is used when `cor_type = "NNGP"`.
#'  \item \code{spatial_coords} A numeric matrix with `nrow(object$data)` rows and two columns containing the spatial coordinates of each site. This is used when `cor_type == "matern"` or `cor_type == "NNGP"`.
#'  \item \code{common} When set to `TRUE`, the same working correlation structure is assumed for all species. When set to `FALSE`, a separate working correlation structure is estimated for each archetype.
#'   }
#' @param num_cores `NULL` for sequential computation, or a positive integer specifying the number of cores to use for parallel computation.
#'
#' @details
#' Note this is a **one-step estimate**: the regression coefficients and species groupings in `object` are not updated.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{working_correlation} If `working_correlation_control$common = TRUE`, then a single estimated working correlation matrix with `nrow(y)` rows and columns. If `working_correlation_control$common = FALSE`, then a list of length `num_archetypes`, where each element is an estimated working correlation matrix with `nrow(y)` rows and columns.
#'   \item \code{working_correlation_params} The estimated parameters defining the working correlation structure in the case where `cor_type = "matern"` or `cor_type = "NNGP"`. If `working_correlation_control$common = TRUE`, then it is a vector where the first element is the estimated spatial range parameter, and the second element is the estimated proportion of the marginal variance due to spatial correlation. A value closer to one indicates a stronger spatial correlation, while the limiting case of zero corresponds to a pure nugget effect. If `working_correlation_control$common = FALSE`, then it is a list of length `num_archetypes`, where each element is a 2-vector as described.
#'   }
#'
#' @seealso
#' [SAGEE()]
#'
#' @author
#' Francis K. C. Hui \email{fhui28@gmail.com}
#'
#' @importFrom fields stationary.cov spatialProcess
#' @importFrom GpGp order_maxmin find_ordered_nn vecchia_Linv
#' @importFrom Matrix sparseMatrix
#' @importFrom stats median dist cov2cor nlminb
#' @importFrom utils modifyList
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#'
#' @export
#' @md
compute_working_correlation <- function(object,
                                working_correlation_control = list(cor_type = "independence",
                                                                   matern_smoothness = 0.5,
                                                                   num_neighbors = 30,
                                                                   spatial_coords = NULL,
                                                                   common = TRUE),
                                num_cores = NULL) {

    ##---------------------------------
    #' # Checks and set up
    ##---------------------------------
    if(!inherits(object, "sagee"))
        stop("object must be of class 'sagee'.")

    if(!is.null(num_cores)) {
        if(!is.numeric(num_cores) || length(num_cores) != 1 || num_cores < 1 || num_cores != as.integer(num_cores))
            stop("num_cores must be a positive integer.")
        doParallel::registerDoParallel(cores = num_cores)
        on.exit(doParallel::stopImplicitCluster(), add = TRUE)
    }

    num_sites <- nrow(object$species_fitted_values)
    num_species <- ncol(object$species_fitted_values)
    num_archetypes <- length(object$archetype_sizes)
    species_groupings <- as.integer(object$species_groupings)

    working_correlation_control_defaults <- list(cor_type = "independence",
                                                 matern_smoothness = 0.5,
                                                 num_neighbors = 30,
                                                 spatial_coords = NULL,
                                                 common = TRUE)
    working_correlation_control <- utils::modifyList(working_correlation_control_defaults, working_correlation_control)

    working_correlation_control$cor_type <- match.arg(working_correlation_control$cor_type, c("independence", "matern", "NNGP"))
    if(working_correlation_control$cor_type %in% c("matern", "NNGP") && is.null(working_correlation_control$spatial_coords))
        stop("working_correlation_control$spatial_coords must be provided when working_correlation_control$cor_type is not set to independence.")
    if(working_correlation_control$cor_type %in% c("matern", "NNGP") && (!is.numeric(working_correlation_control$matern_smoothness) || length(working_correlation_control$matern_smoothness) != 1 || working_correlation_control$matern_smoothness <= 0))
        stop("working_correlation_control$matern_smoothness must be a positive numeric value.")
    if(working_correlation_control$cor_type == "NNGP" && (!is.numeric(working_correlation_control$num_neighbors) || length(working_correlation_control$num_neighbors) != 1 || working_correlation_control$num_neighbors <= 0))
        stop("working_correlation_control$num_neighbors must be a positive numeric value.")
    if(!is.null(working_correlation_control$spatial_coords)) {
        working_correlation_control$spatial_coords <- as.matrix(working_correlation_control$spatial_coords)
        if(!is.matrix(working_correlation_control$spatial_coords) || ncol(working_correlation_control$spatial_coords) != 2 || nrow(working_correlation_control$spatial_coords) != num_sites)
            stop("working_correlation_control$spatial_coords must be a numeric matrix with two columns and the same number of rows as object$data.")
    }

    working_cor <- NULL
    working_cor_params <- NULL

    ##---------------------------------
    #' # Independence
    ##---------------------------------
    if(working_correlation_control$cor_type == "independence") {
        if(working_correlation_control$common)
            working_cor <- diag(x = 1, nrow = num_sites)
        if(!working_correlation_control$common)
            working_cor <- replicate(num_archetypes, diag(x = 1, nrow = num_sites), simplify = FALSE)

        new_working_cor_params <- NULL
        }

    ##---------------------------------
    #' # Matern and NNGP
    ##---------------------------------
    if(working_correlation_control$cor_type %in% c("matern", "NNGP")) {
        marginal_mean <- object$species_fitted_values
        if(object$family$family[1] %in% c("poisson", "binomial", "gaussian", "Gamma", "inverse.gaussian"))
            get_variance <- object$family$variance(marginal_mean) * matrix(object$dispersion_parameters, nrow = num_sites, ncol = num_species, byrow = TRUE)
        if(object$family$family[1] == "negative_binomial")
            get_variance <- object$family$variance(marginal_mean, phi = matrix(object$dispersion_parameters, nrow = num_sites, ncol = num_species, byrow = TRUE))
        if(object$family$family[1] == "tweedie")
            get_variance <- object$family$variance(marginal_mean, phi = matrix(object$dispersion_parameters, nrow = num_sites, ncol = num_species, byrow = TRUE), power = matrix(object$power_parameters, nrow = num_sites, ncol = num_species, byrow = TRUE))
        get_variance  <- pmax(get_variance, 1e-6) ## protect against variances too close to zero
        pearson_residuals <- object$residuals / sqrt(get_variance)
        pearson_residuals <- scale(pearson_residuals, center = FALSE, scale = TRUE)
        rm(marginal_mean, get_variance)
        }

    if(working_correlation_control$cor_type == "matern") {
        if(working_correlation_control$common) {
            spatial_fit <- fields::spatialProcess(x = working_correlation_control$spatial_coords,
                                                  y = pearson_residuals,
                                                  cov.function = "stationary.cov",
                                                  cov.args = list(Covariance = "Matern", smoothness = working_correlation_control$matern_smoothness),
                                                  #cov.params.start = list(aRange = new_working_cor_params),
                                                  mKrig.args = list(m = 1, find.trA = FALSE))
            #sigma2_hat <- unname(spatial_fit$MLESummary["sigma2"])
            lambda_hat <- unname(spatial_fit$MLESummary["lambda"])
            range_hat <- unname(spatial_fit$MLESummary["aRange"])
            new_working_cor_params <- c(range_hat, 1/(1+lambda_hat))
            names(new_working_cor_params) <- c("spatial_range", "spatial_proportion")

            working_cor <- fields::stationary.cov(x1 = working_correlation_control$spatial_coords,
                                                  Covariance = "Matern",
                                                  smoothness = working_correlation_control$matern_smoothness,
                                                  aRange = range_hat)
            working_cor <- new_working_cor_params[2] * (working_cor + lambda_hat * diag(num_sites))
            working_cor <- 0.5 * (working_cor + t(working_cor))
            #print(paste0("Updated spatial decay parameters: ", round(new_working_cor_params, 4)))
            rm(spatial_fit, lambda_hat, range_hat)
            }
        if(!working_correlation_control$common) {
            new_working_cor_params <- vector("list", num_archetypes)
            working_cor <- vector("list", num_archetypes)

            update_correlation_fn <- function(j) {
                Ind <- which(species_groupings == j)
                if(length(Ind) < 1)
                    return(NULL)

                spatial_fit <- fields::spatialProcess(x = working_correlation_control$spatial_coords,
                                                      y = pearson_residuals[, Ind, drop = FALSE],
                                                      cov.function = "stationary.cov",
                                                      cov.args = list(Covariance = "Matern", smoothness = working_correlation_control$matern_smoothness),
                                                      mKrig.args = list(m = 1, find.trA = FALSE))
                #sigma2_hat <- unname(spatial_fit$MLESummary["sigma2"])
                lambda_hat <- unname(spatial_fit$MLESummary["lambda"])
                range_hat <- unname(spatial_fit$MLESummary["aRange"])

                new_working_cor_params <- c(range_hat, 1/(1+lambda_hat))

                working_cor <- fields::stationary.cov(x1 = working_correlation_control$spatial_coords,
                                                      Covariance = "Matern",
                                                      smoothness = working_correlation_control$matern_smoothness,
                                                      aRange = range_hat)
                working_cor <- new_working_cor_params[2] * (working_cor + lambda_hat * diag(num_sites))
                working_cor <- 0.5 * (working_cor + t(working_cor))

                return(list(new_working_cor_params = new_working_cor_params, working_cor = working_cor))
                }

            if(is.null(num_cores))
                update_results <- lapply(1:num_archetypes, update_correlation_fn)
            if(!is.null(num_cores))
                update_results <- foreach(j = 1:num_archetypes,
                                          .export = c("species_groupings",
                                                      "pearson_residuals",
                                                      "working_correlation_control",
                                                      "num_sites"),
                                          .errorhandling = "stop") %dopar% update_correlation_fn(j)

            for(j in 1:num_archetypes) {
                if(is.null(update_results[[j]]))
                    next;

                new_working_cor_params[[j]] <- update_results[[j]]$new_working_cor_params
                working_cor[[j]] <- update_results[[j]]$working_cor
                }
            rm(update_results, update_correlation_fn)
            }
        }
    if(working_correlation_control$cor_type == "NNGP") {
        nngp_config <- list(ord = GpGp::order_maxmin(working_correlation_control$spatial_coords))
        nngp_config$sp_coords_ord <- working_correlation_control$spatial_coords[nngp_config$ord, , drop = FALSE]
        nngp_config$NNarray <- GpGp::find_ordered_nn(nngp_config$sp_coords_ord, m = working_correlation_control$num_neighbors)
        nngp_config$NNarray_lowertri <- !is.na(nngp_config$NNarray)
        nngp_config$B <- Matrix::sparseMatrix(i = row(nngp_config$NNarray)[nngp_config$NNarray_lowertri],
                                              j = nngp_config$NNarray[nngp_config$NNarray_lowertri],
                                              x = 1:sum(nngp_config$NNarray_lowertri),
                                              dims = c(num_sites, num_sites))
        nngp_config$perm <- as.integer(nngp_config$B@x) # maps Linv[ok] into B's internal (column-major) order
        nngp_config$covfun <- switch(as.character(working_correlation_control$matern_smoothness),
                                     "0.5" = "exponential_isotropic",
                                     "1.5" = "matern15_isotropic",
                                     "2.5" = "matern25_isotropic",
                                     "matern_isotropic")

        .nngp_fit_fn <- function(covfun,
                                 resp_ord,
                                 sp_coords_ord,
                                 mat_smooth,
                                 B,
                                 NNarray,
                                 NNarray_lowertri,
                                 perm,
                                 starting_params) {

            n_sites <- nrow(resp_ord)
            n_resp <- ncol(resp_ord)

            negloglik <- function(th) {
                expanded_cp <- if(covfun == "matern_isotropic") c(1, exp(th[1]), mat_smooth, exp(th[2])) else c(1, exp(th[1]), exp(th[2])) # marginal variance, spatial range, smoothness if required, "nugget". Note the marginal variance here is set to one as the actual marginal variance is profiled out.
                Linv <- GpGp::vecchia_Linv(covparms = expanded_cp,
                                           covfun_name = covfun,
                                           locs = sp_coords_ord,
                                           NNarray = NNarray)
                B@x <- Linv[NNarray_lowertri][perm]
                s2 <- sum((B %*% resp_ord)^2) / (n_sites * n_resp)
                logdetR <- -2 * sum(log(Linv[, 1])) # diag of B is Linv[, 1]

                out <- list(nll = 0.5 * (n_sites * n_resp * log(2*pi*s2) + n_resp * logdetR + n_sites * n_resp),
                            s2 = s2)
                return(out)
                }

            nngp_param_start <- log(starting_params)
            opt <- stats::nlminb(start = nngp_param_start,
                                 objective = function(x) negloglik(x)$nll,
                                 control = list(iter.max = 100, trace = 0))

            return(opt)
            }

        .working_cor_from_par <- function(par,
                                          covfun,
                                          mat_smooth,
                                          locs) {
            cp <- if(covfun == "matern_isotropic") c(1, par[1], mat_smooth, 1/par[2] - 1) else c(1, par[1], 1/par[2] - 1)
            return(stats::cov2cor(getExportedValue("GpGp", covfun)(cp, locs)))
            }

        initial_cor_params <- c(spatial_range = stats::median(stats::dist(working_correlation_control$spatial_coords[sample(1:num_sites, size = min(100, num_sites)), , drop = FALSE])) / 4,
                                spatial_proportion = 0.1)

        if(working_correlation_control$common) {
            update_nngpfit <- .nngp_fit_fn(covfun = nngp_config$covfun,
                                           resp_ord = pearson_residuals[nngp_config$ord, , drop = FALSE],
                                           sp_coords_ord = nngp_config$sp_coords_ord,
                                           mat_smooth = working_correlation_control$matern_smoothness,
                                           B = nngp_config$B,
                                           NNarray = nngp_config$NNarray,
                                           NNarray_lowertri = nngp_config$NNarray_lowertri,
                                           perm = nngp_config$perm,
                                           starting_params = c(initial_cor_params[1], 1/initial_cor_params[2] - 1))

            new_working_cor_params <- c(exp(update_nngpfit$par[1]), 1/(1+exp(update_nngpfit$par[2])))
            names(new_working_cor_params) <- c("spatial_range", "spatial_proportion")

            working_cor <- .working_cor_from_par(c(new_working_cor_params[1], new_working_cor_params[2]),
                                                 covfun = nngp_config$covfun,
                                                 mat_smooth = working_correlation_control$matern_smoothness,
                                                 locs = working_correlation_control$spatial_coords) # The estimated marginal variance is ignored since we want a working correlation matrix, not a covariance matrix
            rm(update_nngpfit)
            }
        if(!working_correlation_control$common) {
            new_working_cor_params <- vector("list", num_archetypes)
            working_cor <- vector("list", num_archetypes)

            update_correlation_fn <- function(j) {
                Ind <- which(species_groupings == j)
                if(length(Ind) < 1)
                    return(NULL)

                update_nngpfit <- .nngp_fit_fn(covfun = nngp_config$covfun,
                                               resp_ord = pearson_residuals[nngp_config$ord, Ind, drop = FALSE],
                                               sp_coords_ord = nngp_config$sp_coords_ord,
                                               mat_smooth = working_correlation_control$matern_smoothness,
                                               B = nngp_config$B,
                                               NNarray = nngp_config$NNarray,
                                               NNarray_lowertri = nngp_config$NNarray_lowertri,
                                               perm = nngp_config$perm,
                                               starting_params = c(initial_cor_params[1], 1/initial_cor_params[2] - 1))

                cw_working_cor_params <- c(exp(update_nngpfit$par[1]), 1/(1+exp(update_nngpfit$par[2])))

                cw_working_cor <- .working_cor_from_par(c(cw_working_cor_params[1], cw_working_cor_params[2]),
                                                        covfun = nngp_config$covfun,
                                                        mat_smooth = working_correlation_control$matern_smoothness,
                                                        locs = working_correlation_control$spatial_coords) # The estimated marginal variance is ignored since we want a working correlation matrix, not a covariance matrix
                rm(update_nngpfit)

                return(list(new_working_cor_params = cw_working_cor_params, working_cor = cw_working_cor))
                }

            if(is.null(num_cores))
                update_results <- lapply(1:num_archetypes, update_correlation_fn)
            if(!is.null(num_cores))
                update_results <- foreach(j = 1:num_archetypes,
                                          .export = c("nngp_config",
                                                      "initial_cor_params",
                                                      "pearson_residuals",
                                                      ".nngp_fit_fn",
                                                      ".working_cor_from_par",
                                                      "num_sites",
                                                      "working_correlation_control"),
                                          .errorhandling = "stop") %dopar% update_correlation_fn(j)

            for(j in 1:num_archetypes) {
                if(is.null(update_results[[j]]))
                    next;

                new_working_cor_params[[j]] <- update_results[[j]]$new_working_cor_params
                working_cor[[j]] <- update_results[[j]]$working_cor
                }
            rm(update_results, update_correlation_fn)
            }
        }

    ##---------------------------------
    #' # Wrap up
    ##---------------------------------
    if(working_correlation_control$common && working_correlation_control$cor_type %in% c("matern", "NNGP")) {
        rownames(working_cor) <- colnames(working_cor) <- rownames(object$species_fitted_values)
        names(new_working_cor_params) <- c("spatial_range", "spatial_proportion")
        }
    if(!working_correlation_control$common && working_correlation_control$cor_type %in% c("matern", "NNGP")) {
        names(working_cor) <- paste0("archetype", 1:num_archetypes)
        names(new_working_cor_params) <- paste0("archetype", 1:num_archetypes)

        for(j in seq_len(num_archetypes)) {
            rownames(working_cor[[j]]) <- colnames(working_cor[[j]]) <- rownames(object$species_fitted_values)
            names(new_working_cor_params[[j]]) <- c("spatial_range", "spatial_proportion")
            }
        }

    out <- list(working_correlation = working_cor,
                working_correlation_params = new_working_cor_params)
    return(out)
    }


