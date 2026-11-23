#' @title Fit a Species Archetype Generalized Estimating Equation Model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Fits a species archetype generalized estimating equation (SAGEE) model to a matrix of multivariate responses. Species are hard-classified into a fixed number of archetypes. Species assigned to the same archetype share a vector of archetypal regression coefficients, while each species also has one or more species-specific regression/dispersion coefficients.
#'
#' For species \eqn{i} assigned to archetype \eqn{g_i}, the marginal linear predictor is
#'
#' \deqn{
#'   \eta_i = X_{\alpha}\alpha_{g_i} + X_{\beta}\beta_i,
#' }
#'
#' where \eqn{\alpha_{g_i}} denotes the vector of coefficients for the assigned species archetype, and \eqn{\beta_i} is the coefficient vector specific to species \eqn{i}. The species-specific component should generally contain the intercept and any covariates whose effects are not intended to characterize the species archetypes.
#'
#' @param y A matrix of responses. Rows represent sites or observational units and columns represent species or response variables.
#' @param Xalpha_formula A one-sided formula specifying the covariates whose coefficients vary among archetypes. The corresponding model matrix is constructed from `data` using [stats::model.matrix()]. This formula *must not* contain an intercept term.
#' @param Xbeta_formula A one-sided or two-sided formula specifying the covariates whose coefficients vary among species. The corresponding model matrix is constructed from `data` using [stats::model.matrix()]. The default is `~ 1`, which gives each species its own intercept.
#' @param data A data frame containing the variables referenced by `Xalpha_formula` and `Xbeta_formula`.
#' @param family A `family` object describing the assumed link function and marginal mean-variance relationship. Currently supported families include Gaussian, binary, Poisson, gamma, inverse Gaussian, along with the negative binomial distribution with a log link function courtesy of the custom [nb2()] function, and the Tweedie distribution with a log link function courtesy of the custom [tweedielogfam()] function.
#' @param num_archetypes A positive integer specifying the number of species archetypes to fit.
#' @param start_alpha An optional matrix of starting values for the archetypal coefficients, with `num_archetypes` rows and `ncol(model.matrix(Xalpha_formula, data))` columns.
#' @param start_beta An optional matrix of starting values for the species-specific coefficients, with `ncol(model.matrix(Xbeta_formula, data))` rows and `ncol(y)` columns.
#' @param working_correlation_control Working correlation structure among the rows of `y` within each species. A named list with the following components:
#' \itemize{
#'  \item \code{cor_type} A character string specifying the type of working correlation structure. Currently available options are `"independence"`, `"matern"`, `"NNGP"`, and `"supplied"`.
#'  \item \code{supply_correlation} A numeric matrix with `nrow(y)` rows and columns, or a list of length `num_archetypes` where each element is a numeric matrix with `nrow(y)` rows and columns. This is used when `cor_type = "supplied"`. Note that the supplied working correlation matrix/set of matrices are assumed to be positive definite and symmetric, with the same dimensions as `nrow(y)`. But this is currently not checked by the SAGEE function!
#'  \item \code{matern_smoothness} A positive numeric value specifying the smoothness parameter of the Matern correlation function. This is used when `cor_type == "matern"` or `cor_type == "NNGP"`. Note for `cor_type = "NNGP"`, setting `matern_smoothness` to be 0.5, 1.5, or 2.5 generally gives faster computation.
#'  \item \code{num_neighbors} A positive integer specifying the number of nearest neighbors to use when constructing the NNGP working correlation structure. This is used when `cor_type = "NNGP"`.
#'  \item \code{spatial_coords} A numeric matrix with `nrow(y)` rows and two columns containing the spatial coordinates of each site. This is used when `cor_type == "matern"` or `cor_type == "NNGP"`.
#'  \item \code{common} When set to `TRUE`, the same working correlation structure is assumed for all species. When set to `FALSE`, a separate working correlation structure is estimated for each archetype.
#'   }
#' @param num_cores `NULL` for sequential computation, or a positive integer specifying the number of cores to use for parallel computation.
#' @param control A named list controlling the outer SAGEE iterations:
#'  \itemize{
#'     \item \code{maxit} Maximum number of outer iterations. Note this is deliberately not set to a very large number, as the SAGEE algorithm typically converges quite quickly unless it "bounces" between local solutions.
#'     \item \code{tol} Convergence tolerance for the relative change in the regression coefficients.
#'     \item \code{grouping_tol} Non-negative threshold for the relative improvement in the pseudo-log-likelihood required before a species is reassigned to a different archetype. Larger values reduce switching between nearly equivalent archetypes, but risk a locally suboptimal grouping.
#'     \item \code{trace} Logical value indicating whether progress from the outer iterations should be printed.
#'     }
#' @param GEE_control A list of control parameters for the inner iterative fitting algorithm passed to [GEEfn()]:
#'   \itemize{
#'     \item \code{maxit} Maximum number of GEE iterations.
#'     \item \code{tol} Convergence tolerance for the relative change in the regression coefficients.
#'     \item \code{step_halving_maxit} Maximum number of step-halving attempts made during each coefficient update.
#'     \item \code{trace} Logical value indicating whether iteration progress should be printed.
#'     }
#'
#'
#' @details
#' SAGEE is similar in spirit to the species archetype models of Dunstan et al. (2011, 2013), but uses a generalized estimating equation framework. That is, instead of explicitly assuming distributions for the species responses and building what is effectively a finite mixture of regression models, SAGEE (only) requires specification of a marginal mean model, a marginal variance function, along with a working correlation matrix describing spatial dependence among sites within each species. The approach can either assume a common working correlation structure for all species, or allow each archetype to have its own working correlation structure.
#'
#' See Ito and Sugasawa (2022) for the original grouped GEE methodology on which the SAGEE algorithm is based. Additional methodological details are provided in the accompanying manuscript, currently in preparation.
#'
#' @return An object of class `"sagee"` containing the following components (not necessarily in order):
#' \itemize{
#'   \item \code{call} The matched function call.
#'   \item \code{converged} Logical value indicating whether the algorithm converged.
#'   \item \code{Xalpha_formula/Xbeta_formula/data/family/working_correlation_control} As per the input arguments.
#'   \item \code{archetypal_coefficients} Estimated matrix of archetypal coefficients, with `num_archetypes` rows and `ncol(model.matrix(Xalpha_formula, data))` columns. Row \eqn{g} contains the estimated coefficient vector for archetype \eqn{g}.
#'   \item \code{species_coefficients} Estimated matrix of species-specific coefficients, with `ncol(model.matrix(Xbeta_formula, data))` rows and `ncol(y)` columns. Column \eqn{i} contains the estimated species-specific coefficient vector for species \eqn{i}.
#'   \item \code{dispersion_parameters} Estimated vector of species-specific dispersion parameters (as appropriate).
#'   \item \code{power_parameters} Estimated vector of species-specific power parameters (as appropriate).
#'   \item \code{species_groupings} Estimated species allocations to archetypes.
#'   \item \code{archetype_sizes} A vector of the number of species assigned to each archetype.
#'   \item \code{pseudo_logLik} The Pearson residual pseudo-log-likelihood of the fitted SAGEE, evaluated at the estimated parameters and species groupings.
#'   \item \code{working_correlation} If `working_correlation_control$common = TRUE`, then a single estimated working correlation matrix with `nrow(y)` rows and columns. If `working_correlation_control$common = FALSE`, then a list of length `num_archetypes`, where each element is an estimated working correlation matrix with `nrow(y)` rows and columns.
#'   \item \code{working_correlation_params} The estimated parameters defining the working correlation structure in the case where `cor_type = "matern"` or `cor_type = "NNGP"`. If `working_correlation_control$common = TRUE`, then it is a vector where the first element is the estimated spatial range parameter, and the second element is the estimated proportion of the marginal variance due to spatial correlation. A value closer to one indicates a stronger spatial correlation, while the limiting case of zero corresponds to a pure nugget effect. If `working_correlation_control$common = FALSE`, then it is a list of length `num_archetypes`, where each element is a 2-vector as described.
#'   \item \code{se_species_coefficients} A matrix with `ncol(model.matrix(Xbeta_formula, data))` rows and `ncol(y)` columns, containing the model-based standard errors of the species-specific coefficients.
#'   \item \code{se_archetypal_coefficients} An array with dimensions `num_archetypes` by `ncol(model.matrix(Xalpha_formula, data))` by 3. The third dimension corresponds to the type of standard error, with the package currently offering model-based, sandwich, and sandwich-adjusted (MacKinnon, 1985) standard errors of the archetypal coefficients.
#'   \item \code{vcov_by_archetype} A list with `num_archetypes` elements, where each element contains covariance matrices associated with the standard errors for the archetypal and species-specific coefficients for species assigned to that archetype.
#'   \item \code{archetype_linear_predictor} A matrix with `nrow(y)` rows and `num_archetypes` columns, where column \eqn{g} contains the estimated linear predictor of archetype \eqn{g} i.e., \eqn{X_{\alpha}\widehat{\alpha}_g}.
#'   \item \code{species_linear_predictor} A matrix with the same dimensions as `y`, where column \eqn{i} contains the full linear predictor for species \eqn{i} i.e., \eqn{X_{\alpha}\widehat{\alpha}_{g_i} + X_{\beta}\widehat{\beta}_i}.
#'   \item \code{species_fitted_values} A matrix with the same dimensions as `y`, containing the fitted marginal means for each species.
#'   \item \code{residuals} A matrix with the same dimensions as `y`, containing raw residuals for each species.
#'   \item \code{pearson_residuals} A numeric matrix with the same dimensions as `y`, containing the Pearson residuals for each species.
#'   \item \code{control} As per the input argument.
#'   \item \code{GEE_control} As per the input argument.
#'   }
#'
#' @note
#' Note sandwich standard errors are not computed for the species-specific regression coefficients, as they are informed by only a single species (response) and hence lack replication for conventional sandwich inference. Standard errors for species-specific coefficients are therefore model-based and rely on the specified marginal variance and working correlation model.
#'
#'
#' @references
#' Dunstan, P. K., Foster, S. D., and Darnell, R. (2011).
#' Model-based grouping of species across environmental gradients.
#' \emph{Ecological Modelling}.
#' \doi{10.1016/j.ecolmodel.2010.11.030}
#'
#' Dunstan, P. K., Foster, S. D., Hui, F. K. C., and Warton, D. I. (2013).
#' Finite mixture of regression modeling for high-dimensional count and
#' biomass data in ecology.
#' \emph{Journal of Agricultural, Biological, and Environmental Statistics},
#' 18, 357--375.
#' \doi{10.1007/s13253-013-0146-x}
#'
#' Ito, T. and Sugasawa, S. (2022).
#' Grouped generalized estimating equations for longitudinal data analysis.
#' \emph{Biometrics}.
#' \doi{10.1111/biom.13718}
#'
#' MacKinnon, J. (1985).
#' Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties.
#' \emph{Journal of Econometrics}.
#' \doi{10.1016/0304-4076(85)90158-7}
#'
#' @seealso
#' [fitted.sagee()], [logLik.sagee()], [model.matrix.sagee()], [plot.sagee()], [predict.sagee()], [residuals.sagee()], [summary.sagee()], [vcov.sagee()], [classify.sagee()], [select_curvature()], [compute_working_correlation()]
#'
#' @author
#' Francis K. C. Hui \email{fhui28@gmail.com}
#'
#' @importFrom fields stationary.cov spatialProcess
#' @importFrom GpGp order_maxmin find_ordered_nn vecchia_Linv
#' @importFrom mgcv gam nb tw
#' @importFrom Matrix sparseMatrix
#' @importFrom stats coef glm kmeans model.matrix optim median dist cov2cor nlminb
#' @importFrom utils modifyList
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#'
#' @export
#' @md
SAGEE <- function(y,
                  Xalpha_formula,
                  Xbeta_formula = ~ 1,
                  data,
                  family,
                  num_archetypes,
                  start_alpha = NULL,
                  start_beta = NULL,
                  working_correlation_control = list(cor_type = "independence",
                                                     supply_correlation = NULL,
                                                     matern_smoothness = 0.5,
                                                     num_neighbors = 30,
                                                     spatial_coords = NULL,
                                                     common = TRUE),
                  num_cores = NULL,
                  control = list(maxit = 30, tol = 1e-5, grouping_tol = 1e-2, trace = FALSE),
                  GEE_control = list(maxit = 30, tol = 1e-5, step_halving_maxit = 30, trace = FALSE)) {

    ##---------------------------------
    #' # Initial set up
    ##---------------------------------
    weights = rep(1, ncol(y))

    if(!is.null(num_cores)) {
        if(!is.numeric(num_cores) || length(num_cores) != 1 || num_cores < 1 || num_cores != as.integer(num_cores)) {
            stop("num_cores must be a positive integer.")
            }
        message("Using ", num_cores, " cores for parallel computing...")

        doParallel::registerDoParallel(cores = num_cores)
        on.exit(doParallel::stopImplicitCluster(), add = TRUE)
        }

    Xmat_alpha <- stats::model.matrix(Xalpha_formula, data = data)
    if("(Intercept)" %in% colnames(Xmat_alpha))
        stop("Xalpha_formula must not contain an intercept.")
    Xmat_beta <- stats::model.matrix(Xbeta_formula, data = data)
    p_shared <- ncol(Xmat_alpha)
    p_specific <- ncol(Xmat_beta)
    if(p_shared < 1)
        stop("Xalpha_formula must define at least one non-intercept column.")
    if(p_specific < 1)
        stop("Xbeta_formula must define at least one species-specific column.")

    Xmat_initial <- cbind(Xmat_alpha, Xmat_beta)
    if(qr(Xmat_initial)$rank < ncol(Xmat_initial)) {
        stop("The combined alpha and beta design matrix is rank deficient. Please check for overlapping or linearly dependent terms.")
        }
    rm(Xmat_initial)

    if(!is.matrix(y))
        y <- as.matrix(y)
    if(!is.numeric(y))
        stop("y must be a numeric matrix.")
    if(any(!is.finite(y)))
        stop("y must not contain missing or non-finite values.")
    if(is.null(rownames(y)))
        rownames(y) <- paste0("site", 1:nrow(y))
    if(is.null(colnames(y)))
        colnames(y) <- paste0("species", 1:ncol(y))
    if(nrow(Xmat_alpha) != nrow(y))
        stop("Xmat_alpha and y must have the same number of rows.")
    if(nrow(Xmat_beta) != nrow(y))
        stop("Xmat_beta and y must have the same number of rows.")
    if(!(family$family %in% c("binomial", "poisson", "gaussian", "Gamma", "inverse.gaussian", "negative_binomial", "tweedie")))
        stop("Proposed family are not supported...sorry!")

    num_sites <- nrow(y)
    num_species <- ncol(y)
    if(length(num_archetypes) != 1 || !is.finite(num_archetypes) || num_archetypes < 1 || num_archetypes != as.integer(num_archetypes)) {
        stop("num_archetypes must be a positive integer.")
        }
    num_archetypes <- as.integer(num_archetypes)
    if(num_archetypes > num_species)
        stop("num_archetypes must not exceed the number of columns in y.")

    # if(!is.numeric(weights) || length(weights) != num_species || any(!is.finite(weights)) || any(weights <= 0))
    #     stop("weights must be a numeric vector of length ncol(y) and contain only positive values.")
    # if(working_correlation_control$cor_type == "matern" && length(unique(weights)) > 1)
    #     warning("Weights are currently not accounted for when estimating the spatial correlation parameters. This may be addressed in a future version of SAGEE...")

    working_correlation_control_defaults <- list(cor_type = "independence",
                                                 supply_correlation = NULL,
                                                 matern_smoothness = 0.5,
                                                 num_neighbors = 30,
                                                 spatial_coords = NULL,
                                                 common = TRUE)
    working_correlation_control <- utils::modifyList(working_correlation_control_defaults, working_correlation_control)

    working_correlation_control$cor_type <- match.arg(working_correlation_control$cor_type, c("independence", "matern", "NNGP", "supplied"))
    if(working_correlation_control$cor_type %in% c("matern", "NNGP") && is.null(working_correlation_control$spatial_coords))
        stop("working_correlation_control$spatial_coords must be provided when working_correlation_control$cor_type is not set to independence.")
    if(working_correlation_control$cor_type %in% c("matern", "NNGP") && (!is.numeric(working_correlation_control$matern_smoothness) || length(working_correlation_control$matern_smoothness) != 1 || working_correlation_control$matern_smoothness <= 0))
        stop("working_correlation_control$matern_smoothness must be a positive numeric value.")
    if(working_correlation_control$cor_type == "NNGP" && (!is.numeric(working_correlation_control$num_neighbors) || length(working_correlation_control$num_neighbors) != 1 || working_correlation_control$num_neighbors <= 0))
        stop("working_correlation_control$num_neighbors must be a positive numeric value.")
    if(!is.null(working_correlation_control$spatial_coords)) {
        working_correlation_control$spatial_coords <- as.matrix(working_correlation_control$spatial_coords)
        if(!is.matrix(working_correlation_control$spatial_coords) || ncol(working_correlation_control$spatial_coords) != 2 || nrow(working_correlation_control$spatial_coords) != nrow(y))
            stop("working_correlation_control$spatial_coords must be a numeric matrix with two columns and the same number of rows as y.")
        }

    if(working_correlation_control$cor_type != "supplied" && !is.null(working_correlation_control$supply_correlation))
        warning("working_correlation_control$supply_correlation is ignored when working_correlation_control$cor_type is not set to supplied.")
    if(working_correlation_control$cor_type == "supplied") {
        message("Using a supplied working correlation matrix/set of matrices. Note that the supplied working correlation matrix/set of matrices are assumed to be positive definite and symmetric, with the same dimensions as nrow(y), but this is not checked by the SAGEE function!")

        if(working_correlation_control$common && !is.matrix(working_correlation_control$supply_correlation))
            stop("working_correlation_control$supply_correlation must be a matrix when working_correlation_control$common = TRUE.")
        if(!working_correlation_control$common && !is.list(working_correlation_control$supply_correlation))
            stop("working_correlation_control$supply_correlation must be a list of length num_archetypes when working_correlation_control$common = FALSE.")
        if(!working_correlation_control$common && length(working_correlation_control$supply_correlation) != num_archetypes)
            stop("working_correlation_control$supply_correlation must be a list of length num_archetypes when working_correlation_control$common = FALSE.")
        }


    control_defaults <- list(maxit = 30, tol = 1e-5, grouping_tol = 1e-2, trace = FALSE)
    control <- utils::modifyList(control_defaults, control)

    ##---------------------------------
    #' ## Initial fit and starting values
    ##---------------------------------
    if(working_correlation_control$cor_type == "independence") {
        if(working_correlation_control$common)
            working_cor <- diag(x = 1, nrow = num_sites)
        if(!working_correlation_control$common) {
            message("Setting working_correlation_control$common = FALSE is rather superfluous here since the working correlation is independence for all archetypes. But proceeding nevertheless...")

            working_cor <- replicate(num_archetypes, diag(x = 1, nrow = num_sites), simplify = FALSE)
            names(working_cor) <- paste0("archetype", 1:num_archetypes)
        }

        new_working_cor_params <- NULL
        }
    if(working_correlation_control$cor_type %in% c("matern", "NNGP")) {
        initial_cor_params <- c(stats::median(stats::dist(working_correlation_control$spatial_coords[sample(1:num_sites, size = min(100, num_sites)), , drop = FALSE])) / 4,
                                0.1)
        names(initial_cor_params) <- c("spatial_range", "spatial_proportion")

        if(working_correlation_control$common) {
            new_working_cor_params <- initial_cor_params
            working_cor <- diag(x = 1, nrow = num_sites)
        }
        if(!working_correlation_control$common) {
            new_working_cor_params <- replicate(num_archetypes, initial_cor_params, simplify = FALSE)
            working_cor <- replicate(num_archetypes, diag(x = 1, nrow = num_sites), simplify = FALSE)
            names(new_working_cor_params) <- names(working_cor) <- paste0("archetype", 1:num_archetypes)
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
            }
    }
    if(working_correlation_control$cor_type == "supplied") {
        if(working_correlation_control$common) {
            working_cor <- working_correlation_control$supply_correlation
            }
        if(!working_correlation_control$common) {
            working_cor <- working_correlation_control$supply_correlation
            names(working_cor) <- paste0("archetype", 1:num_archetypes)
            }
        new_working_cor_params <- NULL
        }

    Xmat_initial <- cbind(Xmat_alpha, Xmat_beta)
    all_alpha <- matrix(NA, num_species, p_shared)
    new_beta <- matrix(NA, p_specific, num_species)
    new_disp_param <- rep(1, num_species)
    new_power_param <- rep(1.6, num_species)
    for(i in 1:num_species) {
        if(!(family$family[1] %in% c("negative_binomial", "tweedie")))
            dostartfit <- stats::glm(y[,i] ~ Xmat_initial - 1, family = family)
        if(family$family[1] == "negative_binomial") {
            dostartfit <- mgcv::gam(y[,i] ~ Xmat_initial - 1, family = mgcv::nb())
            new_disp_param[i] <- 1/dostartfit$family$getTheta(trans = TRUE)
            }
        if(family$family[1] == "tweedie") {
            dostartfit <- mgcv::gam(y[,i] ~ Xmat_initial - 1, family = mgcv::tw())
            new_disp_param[i] <- dostartfit$sig2
            new_power_param[i] <- dostartfit$family$getTheta(trans = TRUE)
            }
        all_alpha[i,] <- coef(dostartfit)[1:p_shared]
        new_beta[,i] <- coef(dostartfit)[p_shared + 1:p_specific]

        if(family$family[1] %in% c("gaussian", "Gamma", "inverse.gaussian")) {
            new_disp_param[i] <- summary(dostartfit)$dispersion
            }
        }

    if(is.null(start_alpha)) {
        start_fit <- stats::kmeans(all_alpha, centers = num_archetypes, nstart = 50)
        new_species_groupings <- start_fit$cluster
        new_alpha <- start_fit$centers
        }
     if(!is.null(start_alpha)) {
         start_alpha <- as.matrix(start_alpha)
         if(nrow(start_alpha) != num_archetypes || ncol(start_alpha) != p_shared)
             stop("start_alpha must be a matrix with num_archetypes rows and ncol(model.matrix(Xalpha_formula, data)) columns.")
         new_alpha <- start_alpha
         new_species_groupings <- apply(all_alpha, 1, function(x) which.min(colSums((t(new_alpha) - x)^2)))
        }
    if(!is.null(start_beta)) {
        start_beta <- as.matrix(start_beta)
        if(nrow(start_beta) != p_specific || ncol(start_beta) != num_species)
            stop("start_beta must be a matrix with ncol(model.matrix(Xbeta_formula, data)) rows and ncol(y) columns.")
        new_beta <- start_beta
        }
    rm(dostartfit, all_alpha, Xmat_initial, i)


    .get_working_cor <- function(working_cor,
                                 archetype,
                                 common) {
        if(common)
            return(working_cor)
        if(!common)
            return(working_cor[[archetype]])
        }

    ##---------------------------------
    #' Go go algorithm!
    ##---------------------------------
    converged <- FALSE
    for(k in 1:control$maxit) {
        cw_alpha <- new_alpha
        cw_beta <- new_beta
        cw_theta <- c(as.vector(cw_alpha), as.vector(cw_beta))
        cw_species_groupings <- new_species_groupings

        ##---------------------------------
        #' ## Update beta and species-specific dispersion parameters
        ##---------------------------------
        update_beta_dispersion_fn <- function(j) {
            Ind <- (1:num_species)[cw_species_groupings == j]
            if(length(Ind) >= 1) {
                update_fit <- GEEfn(y = y[, Ind, drop = FALSE],
                                    Xmat_alpha = Xmat_alpha,
                                    Xmat_beta = Xmat_beta,
                                    start_alpha = cw_alpha[j, ],
                                    start_beta = cw_beta[, Ind, drop = FALSE],
                                    family = family,
                                    weights = weights[Ind],
                                    workingcorrelation = .get_working_cor(working_cor = working_cor,
                                                                          archetype = j,
                                                                          common = working_correlation_control$common),
                                    se = FALSE,

                                    control = GEE_control)
                return(update_fit)
                }
            if(length(Ind) < 1)
                return(NULL)
            }

        if(is.null(num_cores))
            update_results <- lapply(1:num_archetypes, update_beta_dispersion_fn)
        if(!is.null(num_cores))
            update_results <- foreach(j = 1:num_archetypes,
                                      .export = c(".get_working_cor",
                                                  "GEEfn",
                                                  "cw_species_groupings",
                                                  "y",
                                                  "Xmat_alpha",
                                                  "Xmat_beta",
                                                  "cw_alpha",
                                                  "cw_beta",
                                                  "family",
                                                  "weights",
                                                  "working_cor",
                                                  "working_correlation_control",
                                                  "GEE_control"),
                                      .errorhandling = "stop") %dopar% update_beta_dispersion_fn(j)

        for(j in 1:num_archetypes) {
            Ind <- (1:num_species)[cw_species_groupings == j]
            if(length(Ind) < 1)
                next;

            new_alpha[j, ] <- update_results[[j]]$alpha
            new_beta[, Ind] <- update_results[[j]]$beta
            new_disp_param[Ind] <- update_results[[j]]$disp_param
            new_power_param[Ind] <- update_results[[j]]$pow_param
            rm(Ind)
            }
        new_theta <- c(as.vector(new_alpha), as.vector(new_beta))
        rm(update_results, update_beta_dispersion_fn)

        ##---------------------------------
        #' ## Update correlation parameter
        #' Might be worth considering using a dampened update, i.e., new correlation parameters or matrix is a weighted combination of the existing and new estimates. Alternatively, update the working correlation less frequently?
        ##---------------------------------
        if(working_correlation_control$cor_type == "independence") { }
        if(working_correlation_control$cor_type == "supplied") { }
        if(working_correlation_control$cor_type == "matern") {
            marginal_mean <- tcrossprod(Xmat_alpha, new_alpha)[, cw_species_groupings, drop = FALSE] + Xmat_beta %*% new_beta
            marginal_mean <- family$linkinv(marginal_mean)
            if(family$family[1] %in% c("poisson", "binomial", "gaussian", "Gamma", "inverse.gaussian"))
                get_variance <- family$variance(marginal_mean) * matrix(new_disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE)
            if(family$family[1] == "negative_binomial")
                get_variance <- family$variance(marginal_mean, phi = matrix(new_disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE))
            if(family$family[1] == "tweedie")
                get_variance <- family$variance(marginal_mean, phi = matrix(new_disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE), power = matrix(new_power_param, nrow = num_sites, ncol = num_species, byrow = TRUE))
            get_variance  <- pmax(get_variance, 1e-6) ## protect against variances too close to zero
            pearson_residuals <- (y - marginal_mean) / sqrt(get_variance)
            pearson_residuals <- scale(pearson_residuals, center = FALSE, scale = TRUE)
            rm(marginal_mean, get_variance)

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
                update_correlation_fn <- function(j) {
                    Ind <- which(cw_species_groupings == j)
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
                                              .export = c("cw_species_groupings",
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
            marginal_mean <- tcrossprod(Xmat_alpha, new_alpha)[, cw_species_groupings, drop = FALSE] + Xmat_beta %*% new_beta
            marginal_mean <- family$linkinv(marginal_mean)
            if(family$family[1] %in% c("poisson", "binomial", "gaussian", "Gamma", "inverse.gaussian"))
                get_variance <- family$variance(marginal_mean) * matrix(new_disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE)
            if(family$family[1] == "negative_binomial")
                get_variance <- family$variance(marginal_mean, phi = matrix(new_disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE))
            if(family$family[1] == "tweedie")
                get_variance <- family$variance(marginal_mean, phi = matrix(new_disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE), power = matrix(new_power_param, nrow = num_sites, ncol = num_species, byrow = TRUE))
            get_variance  <- pmax(get_variance, 1e-6) ## protect against variances too close to zero
            pearson_residuals <- (y - marginal_mean) / sqrt(get_variance)
            pearson_residuals <- scale(pearson_residuals, center = FALSE, scale = TRUE)
            rm(marginal_mean, get_variance)

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
                #     ll <- lapply(1:ncol(pearson_residuals), function(k)
                #         GpGp::vecchia_meanzero_loglik(covparms = exp(th),
                #                                       covfun_name = "exponential_isotropic",
                #                                       y = pearson_residuals_ord[, k],
                #                                       locs = sp_coords_ord,
                #                                       NNarray = NNarray))
                #     return(-sum(unlist(ll)))

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

            if(working_correlation_control$common) {
                # spatial_fit <- autoFRK::autoFRK(Data = pearson_residuals,
                #                                 loc = working_correlation_control$spatial_coords,
                #                                 finescale = FALSE,
                #                                 maxit = 1000)
                # working_cor <- spatial_fit$G %*% tcrossprod(spatial_fit$M, spatial_fit$G) + diag(x = spatial_fit$s, nrow = num_sites)
                # working_cor <- 0.5 * (working_cor + t(working_cor))
                # diag(working_cor) %>% summary
                # working_cor <- cov2cor(working_cor)
                # rm(spatial_fit)

                update_nngpfit <- .nngp_fit_fn(covfun = nngp_config$covfun,
                                               resp_ord = pearson_residuals[nngp_config$ord, , drop = FALSE],
                                               sp_coords_ord = nngp_config$sp_coords_ord,
                                               mat_smooth = working_correlation_control$matern_smoothness,
                                               B = nngp_config$B,
                                               NNarray = nngp_config$NNarray,
                                               NNarray_lowertri = nngp_config$NNarray_lowertri,
                                               perm = nngp_config$perm,
                                               starting_params = c(new_working_cor_params[1], 1/new_working_cor_params[2] - 1))

                new_working_cor_params <- c(exp(update_nngpfit$par[1]), 1/(1+exp(update_nngpfit$par[2])))
                names(new_working_cor_params) <- c("spatial_range", "spatial_proportion")

                working_cor <- .working_cor_from_par(c(new_working_cor_params[1], new_working_cor_params[2]),
                                                     covfun = nngp_config$covfun,
                                                     mat_smooth = working_correlation_control$matern_smoothness,
                                                     locs = working_correlation_control$spatial_coords) # The estimated marginal variance is ignored since we want a working correlation matrix, not a covariance matrix
                rm(update_nngpfit)
                }
            if(!working_correlation_control$common) {
                update_correlation_fn <- function(j) {
                    Ind <- which(cw_species_groupings == j)
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
                                                   starting_params = c(new_working_cor_params[[j]][1], 1/new_working_cor_params[[j]][2] - 1))

                    cw_working_cor_params <- c(exp(update_nngpfit$par[1]), 1/(1+exp(update_nngpfit$par[2])))
                    names(cw_working_cor_params) <- c("spatial_range", "spatial_proportion")

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
                                                          "new_working_cor_params",
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
        #' ## Update grouping assignment -- using max of pseudo-likelihood of the Pearson residuals, and moves species only if it exceeds a threshold. This goes beyond the standard quadratic form of Ito and Sugasawa (2022) as it also accounts for differences in the marginal variance of each species.
        ##---------------------------------
        update_group_assignment_fn <- function(j) {
            updated_fitted_values <- family$linkinv(as.vector(Xmat_alpha %*% new_alpha[j,]) + Xmat_beta %*% new_beta)
            res <- y - updated_fitted_values
            if(family$family[1] %in% c("poisson", "binomial", "gaussian", "Gamma", "inverse.gaussian"))
                marginal_variance <- family$variance(updated_fitted_values) * matrix(new_disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE)
            if(family$family[1] == "negative_binomial")
                marginal_variance <- family$variance(updated_fitted_values, phi = matrix(new_disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE))
            if(family$family[1] == "tweedie")
                marginal_variance <- family$variance(updated_fitted_values, phi = matrix(new_disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE), power = matrix(new_power_param, nrow = num_sites, ncol = num_species, byrow = TRUE))
            marginal_variance <- pmax(marginal_variance, 1e-6) ## protect against variances too close to zero
            pearsonres <- res / sqrt(marginal_variance)

            cw_working_cor <- .get_working_cor(working_cor = working_cor,
                                               archetype = j,
                                               common = working_correlation_control$common)
            Rchol <- t(chol(cw_working_cor))
            whitened_pearsonres <- forwardsolve(Rchol, pearsonres)
            out <- colSums(whitened_pearsonres^2) + colSums(log(marginal_variance)) + 2 * sum(log(diag(Rchol)))
            return(out)
            }

        if(is.null(num_cores))
            QQ <- sapply(1:num_archetypes, update_group_assignment_fn)
        if(!is.null(num_cores))
            QQ <- foreach(j = 1:num_archetypes,
                          .combine = "cbind",
                          .export = c(".get_working_cor",
                                      "y",
                                      "Xmat_alpha",
                                      "Xmat_beta",
                                      "new_alpha",
                                      "new_beta",
                                      "new_disp_param",
                                      "new_power_param",
                                      "family",
                                      "num_sites",
                                      "num_species",
                                      "working_cor",
                                      "working_correlation_control"),
                          .errorhandling = "stop") %dopar% update_group_assignment_fn(j)

        QQ <- as.matrix(QQ)
        candidate_species_groupings <- max.col(-QQ, ties.method = "first")

        #' Compute relative change in negative pseudo logL (smaller QQ values indiciate a more preferred archetype) between the new and existing groups, and move species only if that change exceeds a certain threshold
        cw_QQ <- QQ[cbind(1:num_species, cw_species_groupings)]
        candidate_QQ <- QQ[cbind(1:num_species, candidate_species_groupings)]
        relative_improvement <- (cw_QQ - candidate_QQ) / pmax(abs(cw_QQ), 1)
        species_to_switch <- which(candidate_species_groupings != cw_species_groupings & relative_improvement > control$grouping_tol)
        new_species_groupings <- cw_species_groupings
        new_species_groupings[species_to_switch] <- candidate_species_groupings[species_to_switch]

        #' Check for empty archetypes after species reassignment
        new_archetype_sizes <- tabulate(new_species_groupings, nbins = num_archetypes)
        if(any(new_archetype_sizes == 0))
            stop("One or more archetypes became empty during species reassignment. Please try fewer archetypes, different starting values, a larger control$grouping_tol, etc...")

        # message("Relative_improvement: ", paste(signif(relative_improvement, 4), collapse = ", "))
        # message("Reassigned species: ", paste(species_to_switch, collapse = ", "))
        # message("Relative improvement gaps: ", paste(signif(relative_improvement[species_to_switch], 4), collapse = ", "))
        pseudologLik <- -0.5 * sum(QQ[cbind(1:num_species, new_species_groupings)])
        rm(QQ, cw_QQ, candidate_QQ, relative_improvement, species_to_switch, update_group_assignment_fn)


        ##---------------------------------
        #' ## Convergence check
        ##---------------------------------
        err <- sum(abs(new_theta - cw_theta)) / sum(abs(cw_theta) + 1e-5) ## Relative change in regression coefficients
        # print(new_working_cor_params)

        if(control$trace)
            message("Iteration: ", k, "\t Relative change in regression coefficients: ", round(err, 4))
        if(err < control$tol) {
            converged <- TRUE
            break;
            }
        }

    if(!converged)
        warning("SAGEE did not converge after ", control$maxit, " iterations.")


    ##---------------------------------
    #' # Wrap up
    ##---------------------------------
    #' ## Final fits
    update_beta_dispersion_fn <- function(j) {
        Ind <- (1:num_species)[new_species_groupings == j]
        if(length(Ind) >= 1) {
            update_fit <- GEEfn(y = y[, Ind, drop = FALSE],
                                Xmat_alpha = Xmat_alpha,
                                Xmat_beta = Xmat_beta,
                                family = family,
                                start_alpha = new_alpha[j, ],
                                start_beta = new_beta[, Ind, drop = FALSE],
                                weights = weights[Ind],
                                workingcorrelation = .get_working_cor(working_cor = working_cor,
                                                                      archetype = j,
                                                                      common = working_correlation_control$common),
                                se = TRUE,
                                return_vcov = TRUE,
                                control = GEE_control)
            return(update_fit)
            }
        if(length(Ind) < 1)
            return(NULL)
        }

    if(is.null(num_cores))
        final_results <- lapply(1:num_archetypes, update_beta_dispersion_fn)
    if(!is.null(num_cores))
        final_results <- foreach(j = 1:num_archetypes,
                                 .export = c(".get_working_cor",
                                             "GEEfn",
                                             "new_species_groupings",
                                             "y",
                                             "Xmat_alpha",
                                             "Xmat_beta",
                                             "new_alpha",
                                             "new_beta",
                                             "family",
                                             "weights",
                                             "working_cor",
                                             "working_correlation_control",
                                             "GEE_control"),
                                 .errorhandling = "stop") %dopar% update_beta_dispersion_fn(j)

    save_ses <- vector("list", num_archetypes)
    for(j in 1:num_archetypes) {
        Ind <- (1:num_species)[new_species_groupings == j]
        if(length(Ind) < 1)
            next;

        new_alpha[j, ] <- final_results[[j]]$alpha
        new_beta[, Ind] <- final_results[[j]]$beta
        new_disp_param[Ind] <- final_results[[j]]$disp_param
        new_power_param[Ind] <- final_results[[j]]$pow_param

        save_ses[[j]] <- list(se_model = final_results[[j]]$se_model,
                              se_sandwich_alpha = final_results[[j]]$se_sandwich_alpha,
                              se_sandwich_alpha_adj = final_results[[j]]$se_sandwich_alpha_adj,
                              vcov_model = final_results[[j]]$vcov_model,
                              vcov_sandwich_alpha = final_results[[j]]$vcov_sandwich_alpha,
                              vcov_sandwich_alpha_adj = final_results[[j]]$vcov_sandwich_alpha_adj)
        }
    rm(final_results, update_beta_dispersion_fn, Ind, j)

    .compute_final_QQ_fn <- function(j) {
        updated_fitted_values <- family$linkinv(as.vector(Xmat_alpha %*% new_alpha[j,]) + Xmat_beta %*% new_beta)
        res <- y - updated_fitted_values
        if(family$family[1] %in% c("poisson", "binomial", "gaussian", "Gamma", "inverse.gaussian"))
            marginal_variance <- family$variance(updated_fitted_values) * matrix(new_disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE)
        if(family$family[1] == "negative_binomial")
            marginal_variance <- family$variance(updated_fitted_values, phi = matrix(new_disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE))
        if(family$family[1] == "tweedie")
            marginal_variance <- family$variance(updated_fitted_values, phi = matrix(new_disp_param, nrow = num_sites, ncol = num_species, byrow = TRUE), power = matrix(new_power_param, nrow = num_sites, ncol = num_species, byrow = TRUE))
        marginal_variance <- pmax(marginal_variance, 1e-6) ## protect against variances too close to zero
        pearsonres <- res / sqrt(marginal_variance)

        cw_working_cor <- .get_working_cor(working_cor = working_cor,
                                           archetype = j,
                                           common = working_correlation_control$common)
        Rchol <- t(chol(cw_working_cor))
        whitened_pearsonres <- forwardsolve(Rchol, pearsonres)
        out <- colSums(whitened_pearsonres^2) + colSums(log(marginal_variance)) + 2 * sum(log(diag(Rchol)))
        return(out)
        }

    if(is.null(num_cores))
        QQ <- sapply(1:num_archetypes, .compute_final_QQ_fn)
    if(!is.null(num_cores))
        QQ <- foreach(j = 1:num_archetypes,
                      .combine = "cbind",
                      .export = c(".get_working_cor",
                                  "y",
                                  "Xmat_alpha",
                                  "Xmat_beta",
                                  "new_alpha",
                                  "new_beta",
                                  "new_disp_param",
                                  "new_power_param",
                                  "family",
                                  "num_sites",
                                  "num_species",
                                  "working_cor",
                                  "working_correlation_control"),
                      .errorhandling = "stop") %dopar% .compute_final_QQ_fn(j)

    QQ <- as.matrix(QQ)
    pseudologLik <- -0.5 * sum(QQ[cbind(1:num_species, new_species_groupings)])


    #' ## Formatting
    rownames(new_alpha) <- paste0("archetype", 1:num_archetypes)
    colnames(new_alpha) <- colnames(Xmat_alpha)
    colnames(new_beta) <- colnames(y)
    rownames(new_beta) <- colnames(Xmat_beta)
    names(new_disp_param) <- colnames(y)
    names(new_power_param) <- colnames(y)
    names(new_species_groupings) <- colnames(y)
    if(working_correlation_control$common)
        dimnames(working_cor) <- list(rownames(y), rownames(y))
    if(!working_correlation_control$common) {
        for(j in 1:num_archetypes)
            dimnames(working_cor[[j]]) <- list(rownames(y), rownames(y))
        }

    out <- list(call = match.call(),
                converged = converged,
                Xalpha_formula = Xalpha_formula,
                Xbeta_formula = Xbeta_formula,
                data = data,
                family = family,
                weights = weights,
                working_correlation_control = working_correlation_control,
                archetypal_coefficients = new_alpha,
                species_coefficients = new_beta,
                dispersion_parameters = new_disp_param,
                power_parameters = new_power_param,
                species_groupings = new_species_groupings,
                archetype_sizes = stats::setNames(tabulate(new_species_groupings, nbins = num_archetypes), paste0("archetype", 1:num_archetypes)),
                pseudo_logLik = pseudologLik,
                working_correlation = working_cor,
                working_correlation_params = new_working_cor_params)

    #' ## Standard errors
    out$se_species_coefficients <- matrix(NA, nrow = p_specific, ncol = num_species)
    out$se_archetypal_coefficients <- array(NA,
                                            dim = c(num_archetypes, p_shared, 3),
                                            dimnames = list(paste0("archetype", 1:num_archetypes),
                                                            colnames(Xmat_alpha),
                                                            c("model-based", "sandwich", "sandwich_adjusted")))
    out$vcov_by_archetype <- vector("list", num_archetypes)

    for(k0 in 1:num_archetypes) {
        Ind <- (1:num_species)[out$species_groupings == k0]
        if(length(Ind) < 1)
            next;
        out$se_species_coefficients[, Ind] <- save_ses[[k0]]$se_model$beta
        out$se_archetypal_coefficients[k0,,1] <- save_ses[[k0]]$se_model$alpha
        out$se_archetypal_coefficients[k0,,2] <- save_ses[[k0]]$se_sandwich_alpha
        out$se_archetypal_coefficients[k0,,3] <- save_ses[[k0]]$se_sandwich_alpha_adj

        out$vcov_by_archetype[[k0]] <- list(vcov_model = save_ses[[k0]]$vcov_model,
                                            vcov_sandwich_alpha = save_ses[[k0]]$vcov_sandwich_alpha,
                                            vcov_sandwich_alpha_adj = save_ses[[k0]]$vcov_sandwich_alpha_adj)
        }
    names(out$vcov_by_archetype) <- paste0("archetype", 1:num_archetypes)
    colnames(out$se_species_coefficients) <- colnames(y)
    rownames(out$se_species_coefficients) <- colnames(Xmat_beta)
    rm(save_ses, Ind, k0)

    #' ## Archetype-specific common linear predictors
    out$archetype_linear_predictor <- tcrossprod(Xmat_alpha, new_alpha)
    colnames(out$archetype_linear_predictor) <- paste0("archetype", 1:num_archetypes)

    #' ## Form the full linear predictor corresponding to each species' final archetype.
    out$species_linear_predictor <- out$archetype_linear_predictor[, new_species_groupings, drop = FALSE] + Xmat_beta %*% new_beta
    out$species_fitted_values <- family$linkinv(out$species_linear_predictor)

    out$residuals <- y - out$species_fitted_values
    if(family$family[1] %in% c("poisson", "binomial", "gaussian", "Gamma", "inverse.gaussian"))
        get_variance <- family$variance(out$species_fitted_values) * matrix(out$dispersion_parameters, nrow = num_sites, ncol = num_species, byrow = TRUE)
    if(family$family[1] == "negative_binomial")
        get_variance <- family$variance(out$species_fitted_values, phi = matrix(out$dispersion_parameters, nrow = num_sites, ncol = num_species, byrow = TRUE))
    if(family$family[1] == "tweedie")
        get_variance <- family$variance(out$species_fitted_values, phi = matrix(out$dispersion_parameters, nrow = num_sites, ncol = num_species, byrow = TRUE), power = matrix(out$power_parameters, nrow = num_sites, ncol = num_species, byrow = TRUE))

    out$pearson_residuals <- out$residuals / sqrt(get_variance)
    rm(get_variance)
    dimnames(out$species_linear_predictor) <- dimnames(y)
    dimnames(out$species_fitted_values) <- dimnames(y)
    dimnames(out$residuals) <- dimnames(y)
    dimnames(out$pearson_residuals) <- dimnames(y)

    out$control <- control
    out$GEE_control <- GEE_control
    class(out) <- "sagee"
    return(out)
    }

