#' @title Simulate Data from a Marginal Species Archetype Model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Simulates multivariate responses under a marginal species archetype regression model. In particular, species assigned to the same archetype share a vector of archetypal regression coefficients, while each species also has one or more species-specific regression/dispersion coefficients. Spatial dependence, and subsequently the responses for each species, are generated using a Gaussian copula.
#'
#'
#' For species \eqn{i} assigned to archetype \eqn{g_i}, the marginal linear predictor is
#'
#' \deqn{
#'   \eta_i = X_{\alpha}\alpha_{g_i} + X_{\beta}\beta_i,
#' }
#'
#' where \eqn{\alpha_{g_i}} denotes the vector of coefficients for the assigned species archetype, and \eqn{\beta_i} is the coefficient vector specific to species \eqn{i}. The species-specific component should generally contain the intercept and any covariates whose effects are not intended to characterize the species archetypes.
#'
#' @param Xalpha_formula A one-sided formula specifying the covariates whose coefficients vary among archetypes. The corresponding model matrix is constructed from `data` using [stats::model.matrix()]. This formula *must not* contain an intercept term.
#' @param Xbeta_formula A one-sided or two-sided formula specifying the covariates whose coefficients vary among species. The corresponding model matrix is constructed from `data` using [stats::model.matrix()]. The default is `~ 1`, which gives each species its own intercept.
#' @param data A data frame containing the variables referenced by `Xalpha_formula` and `Xbeta_formula`.
#' @param family A `family` object describing the assumed link function and marginal mean-variance relationship. Currently supported families include those available in [stats::family()], along with the negative binomial distribution with a log link function courtesy of the custom [nb2()] function.
#' @param species_groupings An optional vector of integers specifying the archetype assignment for each species. If not provided, species are assigned based on the `mixing_probabilities` argument.
#' @param mixing_probabilities A vector of probabilities specifying the probability of each species being assigned to each archetype. This must be provided if `species_groupings` is not provided.
#' @param archetypal_coefficients A matrix of archetypal coefficients, with `num_archetypes` rows and `ncol(model.matrix(Xalpha_formula, data))` columns. Row \eqn{g} contains the estimated coefficient vector for archetype \eqn{g}.
#' @param species_coefficients A matrix of species-specific coefficients, with `ncol(model.matrix(Xbeta_formula, data))` rows and `num_species` columns. Column \eqn{i} contains the estimated species-specific coefficient vector for species \eqn{i}.
#' @param dispersion_parameters A vector of species-specific dispersion parameters, required for families with a dispersion parameter (e.g., Gaussian, Gamma, inverse Gaussian, negative binomial).
#' @param power_parameters A vector of species-specific power parameters, required for the Tweedie family.
#' @param common_correlation A logical indicating whether the spatial correlation is common to all species. When set to `TRUE`, the same marginal correlation structure is assumed for all species, in which case `marginal_correlation` must be a single square matrix with dimensions equal to the number of sites. When set to `FALSE`, a separate marginal correlation structure is estimated for each archetype, in which case `marginal_correlation` must be a list of square matrices, one for each archetype.
#' @param marginal_correlation Either a single marginal correlation matrix with `nrow(y)` rows and columns, or a list of length `num_archetypes` where each element is a marginal correlation matrix with `nrow(y)` rows and columns; see the `common_correlation` argument for details.
#' @param seed An optional integer seed for reproducibility.
#' @param control A list of control parameters for the simulation process. Currently supported parameters include
#'  \itemize{
#'     \item \code{max_resp} The maximum value that any response can take. This is used to prevent extreme values in the simulated responses.
#'     \item \code{maxresp_attempts} The maximum number of attempts to generate a response matrix that satisfies the `max_resp` constraint. If this limit is reached, the function will stop and take the last simulated response matrix, even if it exceeds `max_resp`.
#'  }
#'
#' @details
#' The marginal species archetype model is similar in spirit to the species archetype models of Dunstan et al. (2011, 2013), but uses a marginal, Gaussian copula framework to simulate responses so as to make it compatible with the SAGEE approach.
#'
#' See Momozaki et al., (2026) for a similar marginal Gaussian copula approach.
#'
#' @return A named list containing:
#' \itemize{
#'  \item \code{call} The matched call to the function.
#'  \item \code{Y} A data frame of simulated responses. Rows represent sites or observational units and columns represent species or response variables.
#'  \item \code{species_linear_predictor} A matrix with the same dimensions as `Y`, where column \eqn{i} contains the full linear predictor for species \eqn{i} i.e., \eqn{X_{\alpha}\alpha_{g_i} + X_{\beta}\beta_i}.
#'  \item \code{data/family/archetypal_coefficients/species_coefficients/dispersion_parameters/seed} As per the input argument.
#'  \item \code{species_groupings} A vector of integers indicating species allocations to archetypes.
#' }
#'
#' @references
#' Momozaki, T., Nakagawa T., Sugasawa, S., and Solvang H. K. (2026).
#' Semiparametric Copula Estimation for Spatially Correlated Multivariate Mixed Outcomes: Analyzing Visual Sightings of Fin Whales from a Line Transect Survey.
#' \emph{arXiv}.
#' \doi{10.48550/arXiv.2312.12710}
#'
#' @seealso
#' [SAGEE()]
#'
#' @author
#' Francis K. C. Hui \email{fhui28@gmail.com}
#'
#' @importFrom bindata rmvbin
#' @importFrom mvtnorm rmvnorm
#' @importFrom statmod qinvgauss
#' @importFrom stats qnbinom rnorm qgamma qpois
#' @importFrom tweedie qtweedie
#' @importFrom utils modifyList
#'
#' @export
#' @md
create_sageelife <- function(Xalpha_formula,
                             Xbeta_formula = ~ 1,
                             data,
                             family,
                             species_groupings = NULL,
                             mixing_probabilities = NULL,
                             archetypal_coefficients,
                             species_coefficients,
                             dispersion_parameters = NULL,
                             power_parameters = NULL,
                             common_correlation = FALSE,
                             marginal_correlation,
                             seed = NULL,
                             control = list(max_resp = 10000, maxresp_attempts = 20)) {

    ##---------------------------------
    #' # Checks and balances
    ##---------------------------------
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

    num_sites <- nrow(data)
    num_species <- ncol(species_coefficients)
    num_archetypes <- nrow(archetypal_coefficients)
    if(num_archetypes > num_species)
        stop("num_archetypes must not exceed the number of species, i.e., nrow(species_coefficients).")

    if(!(family$family %in% c("binomial", "poisson", "gaussian", "Gamma", "inverse.gaussian", "negative_binomial", "tweedie")))
        stop("Proposed family are supported...sorry!")
    if(family$family[1] %in% c("Gamma", "inverse.gaussian", "tweedie") & family$link[1] != "log")
        stop("Proposed family currently only supports a log link function.")

    if(is.null(species_groupings) & is.null(mixing_probabilities))
        stop("Either species_groupings or mixing_probabilities must be provided.")
    if(!is.null(mixing_probabilities)) {
        if(length(mixing_probabilities) != num_archetypes)
            stop("mixing_probabilities must have length equal to the number of archetypes, i.e., nrow(archetypal_coefficients).")
        if(any(mixing_probabilities < 0) || sum(mixing_probabilities) != 1 || any(mixing_probabilities > 1))
            stop("mixing_probabilities must be non-negative and sum to 1.")
        }

    if(ncol(archetypal_coefficients) != p_shared)
        stop("archetypal_coefficients must have the number of columns equal to the number of shared covariates, i.e., ncol(Xalpha_formula).")
    if(nrow(species_coefficients) != p_specific)
        stop("species_coefficients must have the number of rows equal to the number of species-specific covariates, i.e., ncol(Xbeta_formula).")
    if(!is.null(dispersion_parameters)) {
        if(length(dispersion_parameters) != num_species)
            stop("dispersion_parameters must have length equal to the number of species, i.e., nrow(species_coefficients).")
        }
    if(is.null(dispersion_parameters))
        dispersion_parameters <- rep(1, num_species)
    if(is.null(power_parameters))
        power_parameters <- rep(1.6, num_species)


    if(common_correlation) {
        if(!is.matrix(marginal_correlation) || nrow(marginal_correlation) != num_sites || ncol(marginal_correlation) != num_sites)
            stop("marginal_correlation must be a square matrix with dimensions equal to the number of sites, i.e., nrow(data).")
        }
    if(!common_correlation) {
        if(!is.list(marginal_correlation) || length(marginal_correlation) != num_archetypes)
            stop("marginal_correlation must be a list of length equal to the number of archetypes, i.e., nrow(archetypal_coefficients).")
        for(i in 1:num_archetypes)
            if(!is.matrix(marginal_correlation[[i]]) || nrow(marginal_correlation[[i]]) != num_sites || ncol(marginal_correlation[[i]]) != num_sites)
                stop(paste0("marginal_correlation[[", i, "]] must be a square matrix with dimensions equal to the number of sites, i.e., nrow(data)."))
        }


    control_defaults <- list(max_resp = 10000, maxresp_attempts = 20)
    control <- utils::modifyList(control_defaults, control)

    set.seed(seed)

    ##---------------------------------
    #' # Generate species archetype labels if not provided
    ##---------------------------------
    if(!is.null(species_groupings)) {
        message("Using provided species_groupings...")
        if(length(species_groupings) != num_species)
            stop("species_groupings must have length equal to the number of species, i.e., nrow(species_coefficients).")
        if(any(!species_groupings %in% 1:num_archetypes) || any(species_groupings != as.integer(species_groupings)))
            stop("species_groupings must contain integers between 1 and num_archetypes, i.e., nrow(archetypal_coefficients).")
        }

    if(is.null(species_groupings))
        species_groupings <- sample(1:num_archetypes, num_species, replace = TRUE, prob = mixing_probabilities)


    ##---------------------------------
    #' # Generate species responses
    ##---------------------------------
    true_linpred <- sapply(1:num_species, function(x) {
        as.vector(Xmat_alpha %*% archetypal_coefficients[species_groupings[x], ]) + Xmat_beta %*% species_coefficients[, x]
        })


    .get_species_correlation <- function(i) {
        if (common_correlation)
            return(marginal_correlation)

        return(marginal_correlation[[species_groupings[i]]])
        }

    simY <- matrix(NA, num_sites, num_species)
    for(i in 1:num_species) {
        cw_marginal_correlation <- .get_species_correlation(i)

        if(family$family[1] == "binomial")
            simY[,i] <- as.vector(bindata::rmvbin(1, margprob = family$linkinv(true_linpred[,i]), sigma = cw_marginal_correlation))
        if(family$family[1] == "gaussian")
            simY[,i] <- as.vector(mvtnorm::rmvnorm(1, mean = family$linkinv(true_linpred[,i]), sigma = dispersion_parameters[i] * cw_marginal_correlation))
        if(family$family[1] == "poisson") {
            u <- mvtnorm::rmvnorm(n = 1, mean = rep(0, num_sites), sigma = cw_marginal_correlation)
            simY[,i] <- qpois(pnorm(u), lambda = family$linkinv(true_linpred[,i]))
            try_counter <- 0
            while(try_counter < control$maxresp_attempts && any(simY[,i] > control$max_resp)) {
                u <- mvtnorm::rmvnorm(n = 1, mean = rep(0, num_sites), sigma = cw_marginal_correlation)
                simY[,i] <- qpois(pnorm(u), lambda = family$linkinv(true_linpred[,i]))
                try_counter <- try_counter + 1
                }
            }
        if(family$family[1] == "negative_binomial") {
            u <- mvtnorm::rmvnorm(n = 1, mean = rep(0, num_sites), sigma = cw_marginal_correlation)
            simY[, i] <- qnbinom(pnorm(u), mu = family$linkinv(true_linpred[, i]), size = 1/dispersion_parameters[i])
            try_counter <- 0
            while(try_counter < control$maxresp_attempts && any(simY[,i] > control$max_resp)) {
                u <- mvtnorm::rmvnorm(n = 1, mean = rep(0, num_sites), sigma = cw_marginal_correlation)
                simY[, i] <- qnbinom(pnorm(u), mu = family$linkinv(true_linpred[, i]), size = 1/dispersion_parameters[i])
                try_counter <- try_counter + 1
                }
            }
        if(family$family[1] == "Gamma") {
            u <- mvtnorm::rmvnorm(n = 1, mean = rep(0, num_sites), sigma = cw_marginal_correlation)
            simY[,i] <- qgamma(pnorm(u), shape = 1/dispersion_parameters[i], scale = exp(true_linpred[,i])*dispersion_parameters[i])
            try_counter <- 0
            while(try_counter < control$maxresp_attempts && any(simY[,i] > control$max_resp)) {
                u <- mvtnorm::rmvnorm(n = 1, mean = rep(0, num_sites), sigma = cw_marginal_correlation)
                simY[,i] <- qgamma(pnorm(u), shape = 1/dispersion_parameters[i], scale = exp(true_linpred[,i])*dispersion_parameters[i])
                try_counter <- try_counter + 1
                }
            }
        if(family$family[1] == "inverse.gaussian") {
            u <- mvtnorm::rmvnorm(n = 1, mean = rep(0, num_sites), sigma = cw_marginal_correlation)
            simY[,i] <- statmod::qinvgauss(pnorm(u), dispersion = dispersion_parameters[i], mean = exp(true_linpred[,i]))
            try_counter <- 0
            while(try_counter < control$maxresp_attempts && any(simY[,i] > control$max_resp)) {
                u <- mvtnorm::rmvnorm(n = 1, mean = rep(0, num_sites), sigma = cw_marginal_correlation)
                simY[,i] <- statmod::qinvgauss(pnorm(u), dispersion = dispersion_parameters[i], mean = exp(true_linpred[,i]))
                try_counter <- try_counter + 1
                }
            }
        if(family$family[1] == "tweedie") {
            u <- mvtnorm::rmvnorm(n = 1, mean = rep(0, num_sites), sigma = cw_marginal_correlation)
            simY[,i] <- tweedie::qtweedie(pnorm(u), mu = exp(true_linpred[,i]), phi = dispersion_parameters[i], xi = power_parameters[i])
            try_counter <- 0
            while(try_counter < control$maxresp_attempts && any(simY[,i] > control$max_resp)) {
                u <- mvtnorm::rmvnorm(n = 1, mean = rep(0, num_sites), sigma = cw_marginal_correlation)
                simY[,i] <- tweedie::qtweedie(pnorm(u), mu = exp(true_linpred[,i]), phi = dispersion_parameters[i], xi = power_parameters[i])
                try_counter <- try_counter + 1
                }
            }
    }


    ##---------------------------------
    #' # Bling and output
    ##---------------------------------
    set.seed(NULL)

    rownames(simY) <- rownames(data)
    colnames(simY) <- paste0("species_", 1:num_species)
    rownames(true_linpred) <- rownames(data)
    colnames(true_linpred) <- paste0("species_", 1:num_species)

    out <- list(call = match.call(),
                Y = as.data.frame(simY),
                species_linear_predictor = true_linpred,
                data = as.data.frame(data),
                family = family,
                archetypal_coefficients = archetypal_coefficients,
                species_groupings = species_groupings,
                species_coefficients = species_coefficients,
                dispersion_parameters = dispersion_parameters,
                power_parameters = power_parameters,
                seed = seed)

    return(out)
    }
