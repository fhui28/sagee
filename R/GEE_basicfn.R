#' @title Fit a Multivariate Generalized Estimating Equation Model
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Fits a generalized estimating equation (GEE) model to a matrix of multivariate responses e.g., multi-species abundance data in community ecology. All response variables share a common vector of regression coefficients, while selected response families may have column-specific regression/dispersion and power parameters. The marginal linear predictor for response variable \eqn{j} is
#'
#' \deqn{
#' \eta_j = X_{\alpha}\alpha + X_{\beta}\beta_j,
#' }
#'
#' where \eqn{\alpha} denotes the vector of regression coefficients shared across all response variables, and \eqn{\beta_j} is a vector of coefficients specific to response \eqn{j}. The shared regression coefficients may be optionally penalized using one or more penalty matrices, e.g., for penalized smooth terms.
#'
#' @param y A matrix of responses. Rows represent sites or observational units and columns represent species or response variables.
#' @param Xmat_alpha A model matrix for the shared regression coefficients, with `nrow(y)` rows. The corresponding coefficient vector is shared across all response variables.
#' @param Xmat_beta A model matrix for the response-specific regression coefficients, with `nrow(y)` rows. A separate coefficient vector is estimated for each response variable.
#' @param penalty_alpha An optional positive semi-definite penalty matrix or list of such matrices for the shared regression coefficients. If supplied, the penalized working least-squares problems are solved using [mgcv::pcls()].
#' @param lambda_alpha An optional vector of non-negative smoothing or penalty parameters, with one value corresponding to each matrix in `penalty_alpha`.
#' @param penalty_offset An integer vector giving the *zero-based* starting location of each penalty matrix within the shared coefficient vector \eqn{\alpha}. An offset of zero means that the penalty begins at the first element of \eqn{\alpha}. The length of `penalty_offset` must equal the length of `penalty_alpha`.
#' @param weights A numeric vector of positive weights with length equal to `ncol(y)`.
#' @param workingcorrelation A symmetric, positive-definite numeric matrix describing the working correlation among rows within each response column, with `nrow(y)` rows and columns.
#' @param family A `family` object describing the assumed link function and marginal mean-variance relationship. Currently supported families include Gaussian, binary, Poisson, gamma, inverse Gaussian, along with the negative binomial distribution with a log link function courtesy of the custom [nb2()] function, and the Tweedie distribution with a log link function courtesy of the custom [tweedielogfam()] function.
#' @param start_alpha An optional numeric vector of starting values for the shared regression coefficients, with length equal to `ncol(Xmat_alpha)`.
#' @param start_beta An optional numeric matrix of starting values for the response-specific regression coefficients, with `ncol(Xmat_beta)` rows and `ncol(y)` columns. Column \eqn{j} contains the starting values for response variable \eqn{j}.
#' @param fix_alpha If `TRUE`, then `start_alpha` will be held fixed at the supplied values and not estimated. This is useful when conditionally fitting a response variable against a previously estimated set of shared coefficients. Note if `TRUE`, then `start_alpha` must be supplied.
#' @param se A logical value for whether standard errors should be computed.
#' @param return_vcov A logical value for whether the covariance matrices associated with the model-based and sandwich standard errors should be returned.
#' @param control A list of control parameters for the iterative fitting algorithm:
#'   \itemize{
#'     \item \code{maxit} Maximum number of GEE iterations.
#'     \item \code{tol} Convergence tolerance for the relative change in the regression coefficients.
#'     \item \code{step_halving_maxit} Maximum number of step-halving attempts made during each coefficient update.
#'     \item \code{trace} Logical value indicating whether iteration progress should be printed.
#'     }
#'
#' @details
#' See the accompanying manuscript, currently in preparation, for further methodological details. See also Ito and Sugasawa (2022) for the broader grouped GEE context from which the SAGEE model was developed.
#' 
#' Note the function was designed to be a workhorse function for the main [SAGEE()] function, and so it is not "equipped" with user-friendly features such as formula interfaces, S3 methods etc...
#'
#'
#' @return A list with the following components:
#' \itemize{
#' \item \code{alpha} Estimated regression coefficients shared across all response variables.
#' \item \code{beta} Matrix with `ncol(Xmat_beta)` rows and `ncol(y)` columns, where column \eqn{j} contains the estimated response-specific regression coefficients for response variable \eqn{j}.
#' \item \code{disp_param} Estimated dispersion parameters, with one value for each response variable. For binomial and Poisson responses, these values remain equal to one.
#' \item \code{pow_param} Estimated power parameters, with one value for each response variable. For responses other than Tweedie, these values remain equal to 1.6.
#' \item \code{iterations} Number of GEE iterations performed.
#' \item \code{converged} Logical value indicating whether the relative change in the coefficients fell below `control$tol`.
#' \item \code{solver} Which solver was used to update the regression coefficients. If no penalty is supplied, then a direct solve is used. If a penalty is supplied, then the penalized least squares solver from [mgcv::pcls()] is used.
#' \item \code{penalty_value} The value of the penalty term.
#' \item \code{linear_predictor} A matrix of fitted linear predictors, with the same dimensions as `y`.
#' \item \code{fitted_values} A matrix of fitted marginal means, with the same dimensions as `y`.
#' \item \code{residuals} Matrix of raw residuals, calculated as the observed response matrix minus the fitted marginal means.
#' \item \code{variance} Matrix of fitted marginal variances, with the same dimensions as `y`.
#' \item \code{pearson_residuals} Matrix of Pearson residuals, calculated by dividing the raw residuals by the square root of the fitted marginal variances.
#' \item \code{se_model} If `se = TRUE`, then a list containing the model-based (naive) standard errors for the shared regression coefficients (note these are `NA` when `fix_alpha = TRUE`) and the response-specific regression coefficients. Note if a penalty is included for the shared regression coefficients, the model-based standard errors are computed using the penalized Hessian matrix.
#' \item \code{se_sandwich_alpha} If `se = TRUE`, then a numeric vector of sandwich standard errors for the shared regression coefficients. The standard errors are `NA` when `fix_alpha = TRUE`.
#' \item \code{se_sandwich_alpha_adj} If `se = TRUE`, then a numeric vector of sandwich standard errors for the shared regression coefficients, adjusted for degrees of freedom in the style suggested by MacKinnon (1985). Note the adjustment does not attempt to account for when a penalty is applied to the shared regression coefficients. The standard errors are `NA` when `fix_alpha = TRUE`.
#' \item \code{vcov_model/vcov_sandwich_alpha/vcov_sandwich_alpha_adj} If `return_vcov = TRUE`, then these are the covariance matrices corresponding to `se_model`, `se_sandwich_alpha`, and `se_sandwich_alpha_adj`. `vcov_model` covers all coefficients, ordered as \eqn{\alpha} followed by \eqn{\beta_1, \beta_2, \ldots}. The sandwich matrices cover \eqn{\alpha} only. When `fix_alpha = TRUE`, the \eqn{\alpha} rows and columns of `vcov_model` and both sandwich matrices are `NA`.
#'   }
#'
#' @note
#' Note sandwich standard errors are not computed for the response-specific regression coefficients, as they are informed by only a single response variable and hence lack the response-level replication for conventional sandwich inference. Standard errors for response-specific coefficients are therefore model-based and rely on the specified marginal variance and working correlation model.
#'
#' @references
#' Ito, T. and Sugasawa, S. (2022).
#' Grouped generalized estimating equations for longitudinal data analysis.
#' \emph{Biometrics}.
#' \doi{10.1111/biom.13718}
#'
#' MacKinnon, J. and White, H. (1985).
#' Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties.
#' \emph{Journal of Econometrics}.
#' \doi{10.1016/0304-4076(85)90158-7}
#'
#'
#' @seealso
#' [SAGEE()]
#'
#' @author
#' Francis K. C. Hui \email{fhui28@gmail.com}
#' 
#' @examples
#' \donttest{
#' ##------------------------------
#' ## Example 1: Binary responses
#' ##------------------------------
#' ## Simulate multi-species presence-absence data from a marginal species archetype model,
#' ## with archetype-specific marginal spatial correlations
#' set.seed(062027)
#' num_sites <- 500
#' num_species <- 60
#' num_archetypes <- 3
#' p <- 8
#' covariate_dat <- mvtnorm::rmvnorm(num_sites, sigma = 0.5^abs(outer(1:p, 1:p, "-")))
#' colnames(covariate_dat) <- paste0("X", 1:p)
#' sp_coords <- cbind(runif(num_sites), runif(num_sites))
#'
#' true_archetypal_coefficients <- rbind(c(1, -1, 0, 0, 0.5, -0.5, 0.25, -0.25),
#'                                       c(-1, 1, -1, 0.5, 0, 0, 0.25, 0),
#'                                       c(1, 0.5, 0, 0, -1, 0, -0.5, -0.25))
#' colnames(true_archetypal_coefficients) <- colnames(covariate_dat)
#' true_spp_intercepts <- runif(num_species, -2.5, 0)
#'
#' ## Marginal spatial correlations for each archetype: exponential plus a nugget
#' true_range <- c(0.25, 0.35, 0.55)
#' true_spvarprop <- c(0.6, 0.5, 0.4)
#' true_R <- lapply(1:num_archetypes, function(k0) {
#'     true_spvarprop[k0] * fields::Exp.cov(x1 = sp_coords, aRange = true_range[k0]) +
#'         (1 - true_spvarprop[k0]) * diag(num_sites)
#'     })
#'
#' simdat <- create_sageelife(Xalpha_formula = ~ . - 1,
#'                            Xbeta_formula = ~ 1,
#'                            data = as.data.frame(covariate_dat),
#'                            family = binomial(),
#'                            mixing_probabilities = c(0.5, 0.25, 0.25),
#'                            archetypal_coefficients = true_archetypal_coefficients,
#'                            species_coefficients = matrix(true_spp_intercepts, nrow = 1),
#'                            common_correlation = FALSE,
#'                            marginal_correlation = true_R,
#'                            seed = 062027)
#' summary(colMeans(simdat$Y)) # Species prevalences
#'
#'
#' ## Fit a GEE separately to the species belonging to each archetype, using the true
#' ## archetype memberships and the true marginal correlation as the working correlation.
#' ## Starting values for shared coefficients are obtained using species-specific GLMs.
#' ## This is used purely to demonstrate how the function can be used in practice...
#' fit_binary <- lapply(1:num_archetypes, function(k0) {
#'     get_resps <- simdat$Y[, simdat$species_groupings == k0, drop = FALSE]
#'     get_start <- sapply(1:ncol(get_resps), function(j) {
#'         stats::glm(get_resps[,j] ~ covariate_dat, family = binomial())$coefficients[-1]
#'         })
#'         
#'     GEEfn(y = get_resps,
#'           Xmat_alpha = covariate_dat,
#'           Xmat_beta = matrix(1, nrow = num_sites, ncol = 1),
#'           start_alpha = rowMeans(get_start),
#'           workingcorrelation = true_R[[k0]],
#'           family = binomial(),
#'           control = list(trace = TRUE))
#'     })
#'
#'
#' ## Compare the true and estimated archetypal coefficients, along with sandwich standard errors
#' lapply(1:num_archetypes, function(k0) {
#'     rbind(true = true_archetypal_coefficients[k0,],
#'           estimate = fit_binary[[k0]]$alpha,
#'           se_sandwich = fit_binary[[k0]]$se_sandwich_alpha,
#'           se_sandwich_adj = fit_binary[[k0]]$se_sandwich_alpha_adj)
#'     })
#'
#'
#' ##------------------------------
#' ## Example 2: Negative binomial responses
#' ##------------------------------
#' ## Simulate multi-species overdispersed count data, reusing the covariates, 
#' ## spatial coordinates, and archetypal coefficients from Example 1. 
#' true_disp_param <- runif(num_species, 0.5, 5)
#' true_range <- c(0.25, 0.35, 0.55) * 0.5
#' true_R <- lapply(1:num_archetypes, function(k0) {
#'     true_spvarprop[k0] * fields::Exp.cov(x1 = sp_coords, aRange = true_range[k0]) +
#'         (1 - true_spvarprop[k0]) * diag(num_sites)
#'     })
#'
#' simdat <- create_sageelife(Xalpha_formula = ~ . - 1,
#'                            Xbeta_formula = ~ 1,
#'                            data = as.data.frame(covariate_dat),
#'                            family = nb2(),
#'                            mixing_probabilities = c(0.5, 0.25, 0.25),
#'                            archetypal_coefficients = true_archetypal_coefficients,
#'                            species_coefficients = matrix(true_spp_intercepts, nrow = 1),
#'                            dispersion_parameters = true_disp_param,
#'                            common_correlation = FALSE,
#'                            marginal_correlation = true_R,
#'                            seed = 062027)
#' summary(colMeans(simdat$Y > 0)) # Proportion of non-zero counts per species
#'
#'
#' fit_nb <- lapply(1:num_archetypes, function(k0) {
#'     get_resps <- simdat$Y[, simdat$species_groupings == k0, drop = FALSE]
#'     get_start <- sapply(1:ncol(get_resps), function(j) {
#'         mgcv::gam(get_resps[,j] ~ covariate_dat, family = mgcv::nb())$coefficients[-1]
#'         })
#'
#'     GEEfn(y = get_resps,
#'           Xmat_alpha = covariate_dat,
#'           Xmat_beta = matrix(1, nrow = num_sites, ncol = 1),
#'           start_alpha = rowMeans(get_start),
#'           workingcorrelation = true_R[[k0]],
#'           family = nb2(),
#'           control = list(trace = TRUE))
#'     })
#'
#'
#' ## Compare the true and estimated archetypal coefficients, along with sandwich standard errors
#' lapply(1:num_archetypes, function(k0) {
#'     rbind(true = true_archetypal_coefficients[k0,],
#'           estimate = fit_nb[[k0]]$alpha,
#'           se_sandwich = fit_nb[[k0]]$se_sandwich_alpha,
#'           se_sandwich_adj = fit_nb[[k0]]$se_sandwich_alpha_adj)
#'     })
#'     
#' }
#'
#' @importFrom mgcv pcls
#' @importFrom stats binomial optimize lm.fit
#' @importFrom utils modifyList
#'
#' @export
#' @md
GEEfn <- function(y,
                  Xmat_alpha,
                  Xmat_beta,
                  penalty_alpha = NULL,
                  lambda_alpha = NULL,
                  penalty_offset = NULL,
                  weights = rep(1, ncol(y)),
                  workingcorrelation,
                  family = binomial(),
                  start_alpha = NULL,
                  start_beta = NULL,
                  fix_alpha = FALSE,
                  se = TRUE,
                  return_vcov = FALSE,
                  control = list(maxit = 30, tol = 1e-5, step_halving_maxit = 30, trace = FALSE)) {

    ##---------------------------------
    #' # Initial set up
    ##---------------------------------
    if(!is.matrix(y))
        y <- as.matrix(y)

    if(!is.matrix(Xmat_alpha))
        Xmat_alpha <- as.matrix(Xmat_alpha)
    if(ncol(Xmat_alpha) < 1)
        stop("Xmat_alpha must have at least one column.")

    if(!is.matrix(Xmat_beta))
        Xmat_beta <- as.matrix(Xmat_beta)
    if(ncol(Xmat_beta) < 1)
        stop("Xmat_beta must have at least one column.")

    if(nrow(Xmat_alpha) != nrow(y))
        stop("Xmat_alpha and y must have the same number of rows.")
    if(nrow(Xmat_beta) != nrow(y))
        stop("Xmat_beta and y must have the same number of rows.")

    has_penalty <- !is.null(penalty_alpha) && length(penalty_alpha) > 0
    if(!has_penalty) {
        if(!is.null(lambda_alpha))
            warning("lambda_alpha will be ignored since penalty_alpha is not supplied.")
        if(!is.null(penalty_offset))
            warning("penalty_offset will not be used since penalty_alpha is not supplied.")

        penalty_alpha <- list()
        lambda_alpha <- numeric(0)
        penalty_offset <- integer(0)
        }
    if(has_penalty) {
        if(!is.list(penalty_alpha))
            penalty_alpha <- list(as.matrix(penalty_alpha))
        if(is.list(penalty_alpha))
            penalty_alpha <- lapply(penalty_alpha, as.matrix)

        if(is.null(lambda_alpha))
            stop("lambda_alpha must be supplied when penalty_alpha is supplied.")
        if(length(lambda_alpha) != length(penalty_alpha))
            stop("lambda_alpha must contain one value for each penalty matrix.")
        if(is.null(penalty_offset))
            stop("penalty_offset must be supplied when penalty_alpha is supplied.")
        if(length(penalty_offset) != length(penalty_alpha))
            stop("penalty_offset must contain one value for each penalty matrix.")


        for(r in 1:length(penalty_alpha)) {
            S_r <- penalty_alpha[[r]]

            if(nrow(S_r) != ncol(S_r))
                stop("penalty_alpha[[", r, "]] must be a square numeric matrix.")

            alpha_idx_r <- penalty_offset[r] + 1:(ncol(S_r))
            if(any(alpha_idx_r < 1) || any(alpha_idx_r > ncol(Xmat_alpha)))
                stop("penalty_offset[", r, "] is invalid for penalty_alpha[[", r, "]].")
            }
        }

    num_penalties <- length(penalty_alpha)
    lambda_alpha <- as.numeric(lambda_alpha)
    penalty_offset <- as.integer(penalty_offset)

    if(!is.matrix(workingcorrelation))
        workingcorrelation <- as.matrix(workingcorrelation)
    if(!all(dim(workingcorrelation) == c(nrow(y), nrow(y))))
        stop("workingcorrelation must be nrow(y) by nrow(y).")
    # if(any(abs(workingcorrelation - t(workingcorrelation)) > sqrt(.Machine$double.eps)))
    #     stop("workingcorrelation must be symmetric.")

    if(!is.numeric(weights) || length(weights) != ncol(y) || any(!is.finite(weights)) || any(weights <= 0))
        stop("weights must be a numeric vector of length ncol(y) and contain only positive values.")

    if(!is.null(start_alpha) && length(start_alpha) != ncol(Xmat_alpha))
        stop("start_alpha must have length ncol(Xmat_alpha).")
    if(!is.null(start_beta)) {
        start_beta <- as.matrix(start_beta)
        if(nrow(start_beta) != ncol(Xmat_beta) || ncol(start_beta) != ncol(y))
            stop("start_beta must be matrix with the number of rows equal to ncol(Xmat_beta) and the number of columns equal to ncol(y).")
        }

    if(fix_alpha && is.null(start_alpha))
        stop("start_alpha must be supplied when fix_alpha = TRUE.")

    if(!(family$family %in% c("binomial", "poisson", "gaussian", "negative_binomial", "Gamma", "inverse.gaussian", "tweedie")))
        stop("Proposed family is not supported...sorry!")


    R_chol <- tryCatch(chol(workingcorrelation), error = function(e) { stop("workingcorrelation must be positive definite.") })
    num_sites <- nrow(y)
    num_species <- ncol(y)
    p_shared <- ncol(Xmat_alpha)
    p_specific <- ncol(Xmat_beta)
    p_total <- p_shared + p_specific * num_species

    if(is.null(start_alpha))
        new_alpha <- numeric(p_shared)
    if(!is.null(start_alpha))
        new_alpha <- as.vector(start_alpha)
    if(is.null(start_beta))
        new_beta <-  matrix(0, nrow = p_specific, ncol = num_species, byrow = FALSE)
    if(!is.null(start_beta))
        new_beta <- as.matrix(start_beta)
    new_disp_param <- rep(1, num_species)
    new_pow_param <- rep(1.6, num_species)

    control_defaults <- list(maxit = 30, tol = 1e-5, step_halving_maxit = 30, trace = FALSE)
    control <- utils::modifyList(control_defaults, control)

    if(!se && return_vcov)
        stop("return_vcov = TRUE is not allowed when se = FALSE.")


    ##---------------------------------
    #' # Define a fast inner function to compute U, M, and Q, where the last is a quadratic root-finding heuristic, and no penalties are included.
    #' When penalties are supplied, then switch to a penalized least squares (pcls) approach using the mgcv::pcls() function.
    #' Also define function to build the penalty matrix for the shared regression coefficients
    ##---------------------------------
    build_penalty_matrix <- function(penalty_alpha,
                                     penalty_offset,
                                     lambda_alpha,
                                     p_shared,
                                     p_total,
                                     fix_alpha) {

        P_alpha <- matrix(0, nrow = p_shared, ncol = p_shared)

        if(!fix_alpha && length(penalty_alpha) > 0) {
            for(r in seq_along(penalty_alpha)) {
                S_r <- penalty_alpha[[r]]
                idx_r <- penalty_offset[r] + 1:ncol(S_r)
                P_alpha[idx_r, idx_r] <- P_alpha[idx_r, idx_r] + lambda_alpha[r] * S_r
                }
            }

            P_full <- matrix(0, nrow = p_total, ncol = p_total)
            P_full[1:p_shared, 1:p_shared] <- P_alpha
            return(list(alpha = P_alpha, full = P_full))
        }

    .directsolve_fn <- function(theta,
                                disp_param,
                                pow_param = NULL,
                                already_all_eta = NULL,
                                return_species_U = FALSE) {

        if(is.null(already_all_eta)) {
            alpha <- theta[1:p_shared]
            Beta <- matrix(theta[-(1:p_shared)], nrow = p_specific, ncol = num_species, byrow = FALSE)

            all_eta <- Xmat_beta %*% Beta
            all_eta <- all_eta + as.vector(Xmat_alpha %*% alpha)
            }
        if(!is.null(already_all_eta))
            all_eta <- already_all_eta

        all_mu <- family$linkinv(all_eta)
        all_residual <- y - all_mu
        U <- numeric(p_total)
        M <- matrix(0, nrow = p_total, ncol = p_total)

        if(return_species_U)
            species_U <- matrix(0, nrow = p_total, ncol = num_species)

        for(j in 1:num_species) {
            specific_idx <- p_shared + (j - 1) * p_specific + 1:p_specific

            if(family$family[1] %in% c("poisson", "binomial", "gaussian", "Gamma", "inverse.gaussian"))
                variance_j <- family$variance(all_mu[,j]) * disp_param[j]
            if(family$family[1] == "negative_binomial")
                variance_j <- family$variance(all_mu[,j], phi = disp_param[j])
            if(family$family[1] == "tweedie")
                variance_j <- family$variance(all_mu[,j], power = pow_param[j], phi = disp_param[j])

            variance_j <- pmax(variance_j, 1e-6) # Prevent numerical issues with very small variances
            dmudeta_j <- family$mu.eta(all_eta[,j])
            Ajinv_sqrt_dmudeta_j <- dmudeta_j / sqrt(variance_j)

            ## Derivative matrix for response j
            ## Columns: [alpha | beta_1 | beta_2 | ... | beta_m]
            D_scaled <- cbind(Xmat_alpha * Ajinv_sqrt_dmudeta_j, Xmat_beta * Ajinv_sqrt_dmudeta_j)
            Rinv_all <- Rinv %*% cbind(D_scaled, all_residual[,j] / sqrt(variance_j))

            ## Accumulate first and second derivatives for NR update
            cw_U <- weights[j] * drop(crossprod(D_scaled, Rinv_all[, ncol(Rinv_all)]))
            U[1:p_shared] <- U[1:p_shared] + cw_U[1:p_shared]
            U[specific_idx] <- U[specific_idx] + cw_U[-(1:p_shared)]
            if(return_species_U) {
                species_U[1:p_shared, j] <- cw_U[1:p_shared]
                species_U[specific_idx, j] <- cw_U[-(1:p_shared)]
                }

            cw_M <- as.matrix(weights[j] * crossprod(D_scaled, Rinv_all[, 1:ncol(D_scaled)]))
            all_index <- c(1:p_shared, specific_idx)
            M[all_index, all_index] <- M[all_index, all_index] + cw_M
            }

        if(fix_alpha) {
            U_update <- U[-(1:p_shared)]
            M_update <- M[-(1:p_shared), -(1:p_shared)]
            delta_update <- solve(M_update, U_update)
            delta <- numeric(p_total)
            delta[-(1:p_shared)] <- delta_update

            Q <- as.vector(crossprod(U_update, delta_update))
            }
        if(!fix_alpha) {
            delta <- solve(M, U)
            Q <- as.vector(crossprod(U, delta))
            }

        out <- list(U = U,
                    M = M,
                    U_pen = U,
                    M_pen = M,
                    delta = as.vector(delta),
                    Q = Q)
        if(return_species_U)
            out$species_U <- species_U
        return(out)
        }

    .pcls_fn <- function(theta,
                         disp_param,
                         pow_param = NULL,
                         already_all_eta = NULL,
                         return_species_U = FALSE) {

        alpha <- theta[1:p_shared]
        Beta <- matrix(theta[-(1:p_shared)], nrow = p_specific, ncol = num_species, byrow = FALSE)

        if(is.null(already_all_eta)) {
            all_eta <- Xmat_beta %*% Beta
            all_eta <- all_eta + as.vector(Xmat_alpha %*% alpha)
            }
        if(!is.null(already_all_eta))
            all_eta <- already_all_eta

        all_mu <- family$linkinv(all_eta)

        ## Store species-level blocks because they are also needed for sandwich covariance calculations
        x_blocks <- vector("list", num_species)
        z_blocks <- vector("list", num_species)
        R_chol_t <- t(R_chol)
        if(fix_alpha)
            offset_j <- drop(Xmat_alpha %*% alpha)

        for(j in 1:num_species) {
            specific_idx <- p_shared + (j - 1) * p_specific + 1:p_specific

            if(family$family[1] %in% c("poisson", "binomial", "gaussian", "Gamma", "inverse.gaussian"))
                variance_j <- family$variance(all_mu[,j]) * disp_param[j]
            if(family$family[1] == "negative_binomial")
                variance_j <- family$variance(all_mu[,j], phi = disp_param[j])
           if(family$family[1] == "tweedie")
                variance_j <- family$variance(all_mu[,j], power = pow_param[j], phi = disp_param[j])
            variance_j <- pmax(variance_j, 1e-6) # Prevent numerical issues with very small variances
            dmudeta_j <- family$mu.eta(all_eta[,j])

            z_j <- all_eta[,j] + (y[, j] - all_mu[,j]) / dmudeta_j # Working response: z = g(mu) + (y - mu) * g'(mu)
            D_scaled_j <- dmudeta_j / sqrt(variance_j) # dmu/deta * A^{-1/2}

            if(fix_alpha) {
                X_j <- matrix(0, nrow = num_sites, ncol = p_specific * num_species)
                X_j[, (j - 1) * p_specific + 1:p_specific] <- Xmat_beta
                target_j <- z_j - offset_j
                }
            if(!fix_alpha) {
                X_j <- matrix(0, nrow = num_sites, ncol = p_total)
                X_j[, 1:p_shared] <- Xmat_alpha
                X_j[, specific_idx] <- Xmat_beta
                target_j <- z_j
                }

            X_scaled_j <- X_j * D_scaled_j
            target_scaled_j <- target_j * D_scaled_j

            X_whitened_j <- forwardsolve(R_chol_t, X_scaled_j) # Final working model matrix
            target_whitened_j <- forwardsolve(R_chol_t, target_scaled_j) # Final working response vector

            X_whitened_j <- sqrt(weights[j]) * X_whitened_j
            target_whitened_j <- sqrt(weights[j]) * target_whitened_j

            x_blocks[[j]] <- X_whitened_j
            z_blocks[[j]] <- target_whitened_j
            }

        X_work <- do.call(rbind, x_blocks)
        z_work <- as.vector(unlist(z_blocks))
        if(fix_alpha) {
            fit <- stats::lm.fit(x = X_work, y = z_work, singular.ok = FALSE)
            beta_new_vector <- as.vector(fit$coefficients)
            theta_new <- c(alpha, beta_new_vector)
            }
        if(!fix_alpha) {
            pcls_problem <- list(y = z_work,
                                 X = X_work,
                                 w = rep(1, length(z_work)),
                                 S = penalty_alpha,
                                 off = penalty_offset,
                                 sp = lambda_alpha,
                                 p = as.vector(theta),
                                 C = matrix(0, nrow = 0, ncol = 0),
                                 Ain = matrix(0, nrow = 1, ncol = p_total),
                                 bin = -1)

            theta_new <- as.vector(mgcv::pcls(pcls_problem))
            fit <- list(coefficients = theta_new, pcls_problem = pcls_problem)
            }

        delta <- theta_new - theta

        ## Recover score and sensitivity quantities for compatibility with the rest of the GEE code.
        if(fix_alpha) {
            free_idx <- p_shared + 1:(p_specific * num_species)
            M <- matrix(0, nrow = p_total, ncol = p_total)
            M[free_idx, free_idx] <- crossprod(X_work)
            M_pen <- M

            U <- numeric(p_total)
            beta_current_vector <- theta[-(1:p_shared)]
            U[free_idx] <- as.vector(crossprod(X_work, z_work - drop(X_work %*% beta_current_vector)))
            U_pen <- U

            Q <- drop(crossprod(U, delta))
            }
        if(!fix_alpha) {
            penalty_full <- penalty_object$full
            M <- crossprod(X_work)
            M_pen <- M + penalty_full

            ## This is the score at the input theta, not at the newly fitted theta
            U <- as.vector(crossprod(X_work, z_work - drop(X_work %*% theta)))
            U_pen <- U - as.vector(penalty_full %*% theta)

            Q <- as.vector(crossprod(U_pen, delta))
            }

        out <- list(U = U,
                    M = M,
                    U_pen = U_pen,
                    M_pen = M_pen,
                    delta = delta,
                    theta_new = theta_new,
                    Q = Q)


        if(return_species_U) {
            species_U <- matrix(0, nrow = p_total, ncol = num_species)

            for(j in 1:num_species) {
                if(fix_alpha) {
                    beta_current_vector <- theta[-(1:p_shared)]
                    free_idx <- p_shared + 1:(p_specific * num_species)
                    residual_j <- z_blocks[[j]] - as.vector(x_blocks[[j]] %*% beta_current_vector)
                    species_U[free_idx, j] <- as.vector(crossprod(x_blocks[[j]], residual_j))
                    }
                if(!fix_alpha) {
                    residual_j <- z_blocks[[j]] - as.vector(x_blocks[[j]] %*% theta)
                    species_U[, j] <- as.vector(crossprod(x_blocks[[j]], residual_j))
                    }
                }

            out$species_U <- species_U
            }
        return(out)
        }

    if(has_penalty) {
        selected_inner_fn <- .pcls_fn
        selected_solver <- "pcls"
        }
    if(!has_penalty) {
        selected_inner_fn <- .directsolve_fn
        selected_solver <-  "direct"
        }


    ##---------------------------------
    #' # GEE iterations
    ##---------------------------------
    converged <- FALSE
    err <- Inf
    Q_new <- Inf
    Rinv <- solve(workingcorrelation)
    penalty_object <- build_penalty_matrix(penalty_alpha = penalty_alpha,
                                           penalty_offset = penalty_offset,
                                           lambda_alpha = lambda_alpha,
                                           p_shared = p_shared,
                                           p_total = p_total,
                                           fix_alpha = fix_alpha)

    for(k in 1:control$maxit) {
        cw_theta <- c(new_alpha, as.vector(new_beta))
        cw_disp_param <- new_disp_param
        cw_pow_param <- new_pow_param

        #' ## Update coefficients
        #' STEP HALVING: theta_try = theta + step_size * Delta, where step_size = 1, 1/2, 1/4, ... Accept when Q(theta_try) <= Q(theta)
        current <- selected_inner_fn(theta = cw_theta, disp_param = cw_disp_param, pow_param = cw_pow_param)
        Q_old <- current$Q
        delta <- current$delta

        #accepted <- FALSE
        #' Note if the step size fails to reduce Q even after step_halving_maxit attempts, the algorithm will still accept the last step size and move on to the next iteration!
        step_factor <- 1
        for(ihalf in 0:control$step_halving_maxit) {
            new_theta_try <- cw_theta + step_factor * delta
            trial <- try(selected_inner_fn(theta = new_theta_try, disp_param = cw_disp_param, pow_param = cw_pow_param), silent = TRUE)

            if(inherits(trial, "try-error")) {
                step_factor <- 0.5 * step_factor
                next;
                }

            if(trial$Q <= Q_old) {
                #accepted <- TRUE
                break;
                }
            step_factor <- 0.5 * step_factor
            }

        new_theta <- new_theta_try
        if(fix_alpha)
            new_theta[1:p_shared] <- as.vector(start_alpha)
        new_alpha <- new_theta[1:p_shared]
        new_beta <- matrix(new_theta[-(1:p_shared)], nrow = p_specific, ncol = num_species, byrow = FALSE)
        Q_new <- trial$Q
        new_linear_predictor <- as.vector(Xmat_alpha %*% new_alpha) + Xmat_beta %*% new_beta
        new_fitted_value <- family$linkinv(new_linear_predictor)

        #' ## Update species-specific dispersion parameters
        if(family$family[1] == "gaussian") {
            rawres <- y - new_fitted_value
            new_disp_param <- as.vector(colMeans(rawres^2))
            }
        if(family$family[1] == "Gamma") {
            rawres <- y - new_fitted_value
            new_disp_param <- as.vector(colMeans(rawres^2 / new_fitted_value^2))
            }
        if(family$family[1] == "inverse.gaussian") {
            rawres <- y - new_fitted_value
            new_disp_param <- as.vector(colMeans(rawres^2 / new_fitted_value^3))
            }
        if(family$family[1] %in% c("negative_binomial")) {
            rawres <- y - new_fitted_value
            for(k0 in 1:num_species) {
                profilelogLdisp_param <- function(x) {
                    varvals <- new_fitted_value[,k0] + x*new_fitted_value[,k0]^2
                    out <- sum(-0.5*log(varvals) - 0.5*rawres[,k0]^2/varvals)
                    return(out)
                    }

                update_disp_param <- optimize(f = profilelogLdisp_param, lower = 1e-6, upper = 1e4, maximum = TRUE)
                new_disp_param[k0] <- update_disp_param$maximum
                }
            }
        if(family$family[1] %in% c("tweedie")) {
            rawres <- y - new_fitted_value
            for(k0 in 1:num_species) {
                r2_k <- rawres[, k0]^2
                mu_k <- pmax(new_fitted_value[, k0], 1e-6)
                .phi_given_pow <- function(pow)
                    mean(r2_k / mu_k^pow)
                .profile_pow <- function(pow) {
                    v <- .phi_given_pow(pow) * mu_k^pow
                    sum(-0.5 * log(v) - 0.5 * r2_k / v)
                    }

                update_disppow_param <- optimize(.profile_pow, lower = 1 + 1e-4, upper = 2 - 1e-4, maximum = TRUE)
                new_pow_param[k0] <- update_disppow_param$maximum
                new_disp_param[k0] <- .phi_given_pow(update_disppow_param$maximum)
                }
            }

        if(fix_alpha)
            convergence_indices <- -(1:p_shared)
        if(!fix_alpha)
            convergence_indices <- 1:p_total

        err <- sum(abs(new_theta[convergence_indices] - cw_theta[convergence_indices])) / sum(abs(cw_theta[convergence_indices]) + 1e-5) ## Relative change in regression coefficients
        if(control$trace)
            message("Iteration: ", k, "\t Relative change in regression coefficients: ", round(err, 4))
        if(err < control$tol) {
            converged <- TRUE
            break;
            }

        }

    ##---------------------------------
    #' # Wrap up
    ##---------------------------------
    out <- list(alpha = new_alpha,
                beta = new_beta,
                disp_param = new_disp_param,
                pow_param = new_pow_param,
                iterations = k,
                converged = converged,
                solver = selected_solver,
                penalty_value = 0)

    if(has_penalty) {
        out$penalty_value <- 0
        for(r in 1:num_penalties) {
            idx_r <- penalty_offset[r] + 1:ncol(penalty_alpha[[r]])
            out$penalty_value <- out$penalty_value + 0.5 * lambda_alpha[r] * as.vector(crossprod(new_alpha[idx_r], penalty_alpha[[r]] %*% new_alpha[idx_r]))
            }
        rm(idx_r)
        }


    out$linear_predictor <- new_linear_predictor
    out$fitted_values <- family$linkinv(out$linear_predictor)
    out$residuals <- y - out$fitted_values
    if(family$family[1] %in% c("poisson", "binomial", "gaussian", "Gamma", "inverse.gaussian"))
        out$variance <- family$variance(out$fitted_values) * matrix(out$disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE)
    if(family$family[1] == "negative_binomial")
        out$variance <- family$variance(out$fitted_values, phi = matrix(out$disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE))
    if(family$family[1] == "tweedie")
        out$variance <- family$variance(out$fitted_values,
                                        power = matrix(out$pow_param, nrow = num_sites, ncol = num_species, byrow = TRUE),
                                        phi = matrix(out$disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE))
    out$variance <- pmax(out$variance, 1e-6)
    out$pearson_residuals <- out$residuals / sqrt(out$variance)


    #' ## Standard errors
    if(se) {
        final_gee <- selected_inner_fn(theta = c(new_alpha, as.vector(new_beta)),
                                       disp_param = new_disp_param,
                                       pow_param = new_pow_param,
                                       already_all_eta = new_linear_predictor,
                                       return_species_U = TRUE)
        if(fix_alpha) {
            beta_bread <- final_gee$M_pen[-(1:p_shared), -(1:p_shared), drop = FALSE]
            beta_inv <- solve(beta_bread)
            beta_inv <- 0.5 * (beta_inv + t(beta_inv))
            se_beta <- sqrt(diag(beta_inv))
            out$se_model <- list(alpha = rep(NA, p_shared),
                                 beta = matrix(se_beta, nrow = p_specific, ncol = num_species, byrow = FALSE))
            out$se_sandwich_alpha <- rep(NA, p_shared)
            out$se_sandwich_alpha_adj <- rep(NA, p_shared)

            if(return_vcov) {
                out$vcov_model <- matrix(NA, nrow = p_total, ncol = p_total)
                out$vcov_model[-(1:p_shared), -(1:p_shared)] <- beta_inv
                out$vcov_sandwich_alpha <- out$vcov_sandwich_alpha_adj <- matrix(NA, nrow = p_shared, ncol = p_shared)
                }
            }

        if(!fix_alpha) {
            bread <- final_gee$M_pen
            bread_inv <- solve(bread)
            bread_inv <- 0.5 * (bread_inv + t(bread_inv))
            se_model <- sqrt(diag(bread_inv))
            out$se_model <- list(alpha = se_model[1:p_shared],
                                 beta = matrix(se_model[-(1:p_shared)], nrow = p_specific, ncol = num_species, byrow = FALSE))

            #' The sandwich covariance is used only to compute the standard errors of the common regression coefficients. It does not make much sense IMO to compute the sandwich covariance for species-specific regression coefficients, since they are informed by only a single species and so "meat-ing" them is not meaningful.
            vcov_sandwich <- bread_inv %*% tcrossprod(final_gee$species_U) %*% bread_inv
            vcov_sandwich <- 0.5 * (vcov_sandwich + t(vcov_sandwich))
            se_sandwich <- sqrt(diag(vcov_sandwich))
            out$se_sandwich_alpha <- se_sandwich[1:p_shared]
            #beta = matrix(se_sandwich[-(1:p_shared)], nrow = p_specific, ncol = num_species, byrow = FALSE))
            #out$species_U <- final_gee$species_U

            #' Makes the degrees of freedom adjustment to the sandwich covariance matrix, as suggested by Mackinnon (1985).
            # coefficient_influence <- solve(final_gee$M_pen, final_gee$M)
            # edf_alpha <- sum(diag(coefficient_influence)[1:p_shared])
            df_adjustment <- num_species/(num_species - p_shared)
            if((num_species - p_shared) <= 0) {
                warning("Degrees of freedom adjustment for sandwich covariance is non-positive. Setting to one.")
                df_adjustment <- 1
            }
            vcov_sandwich_adj <- df_adjustment * vcov_sandwich
            se_sandwich_adj <- sqrt(diag(vcov_sandwich_adj))
            out$se_sandwich_alpha_adj <- se_sandwich_adj[1:p_shared]

            if(return_vcov) {
                out$vcov_model <- bread_inv
                out$vcov_sandwich_alpha <- vcov_sandwich[1:p_shared, 1:p_shared, drop = FALSE]
                out$vcov_sandwich_alpha_adj <- vcov_sandwich_adj[1:p_shared, 1:p_shared, drop = FALSE]
                }
            }

        # kcH_U <- kc_scores(theta = c(new_alpha, as.vector(new_beta)),
        #                    disp_param = new_disp_param,
        #                    Rinv = Rinv,
        #                    bread_inv = bread_inv)
        #
        # vcov_sandwich_kc <- bread_inv %*% tcrossprod(kcH_U) %*% bread_inv
        # se_sandwich_kc <- sqrt(diag(vcov_sandwich_kc))
        # out$se_sandwich_kc_alpha <- se_sandwich_kc[1:p_shared]
        }

    return(out)
    }



