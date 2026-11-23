#' @title Construct standard and adjusted predictions from a fitted SAGEE model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Extends the adjusted GEE prediction of Hui, Muller and Welsh (2025) from independent clustered longitudinal data to the case of a species archetype GEE model fitted by [SAGEE()], potentially with a spatial working correlation matrix/set of matrices.
#'
#' @param object An object of class `"sagee"`.
#' @param newdata A data frame containing the variables referenced by the chosen formula in `object` for the new sites.
#' @param new_spatial_coords A numeric matrix with `nrow(newdata)` rows and two columns, giving the spatial coordinates of the new sites on the same coordinate system as `object$working_correlation_control$spatial_coords`.
#' @param intervals Should uncertainty intervals also be produced for the point predictions? Defaults to `FALSE`.
#' @param intervals_quantiles Quantiles to use for constructing uncertainty intervals.
#' @param vcov_type The type of covariance matrix to use for constructing uncertainty intervals. Either `"model-based"` (the default), `"sandwich"`, or `"sandwich_adjusted"`. The latter two options replace the "archetypal block" of the covariance matrix with the corresponding sandwich estimate, while keeping the model-based blocks for the remainder of the matrix; see [vcov.sagee()] for details.
#' @param ... Unused.
#'
#' @details
#' Briefly, the standard GEE predictor is based on the marginal mean model alone, i.e.,
#'
#' \deqn{
#'  \mu_{i} = h^{-1}(X_{\alpha} \alpha_{g_i} + X_{\beta} \beta_i),
#' }
#'
#' where \eqn{\alpha_{g_i}} denotes the vector of coefficients for the assigned species archetype, and \eqn{\beta_i} is the coefficient vector specific to species \eqn{i}. The species-specific component should generally contain the intercept and any covariates whose effects are not intended to characterize the species archetypes. Finally, \eqn{h^{-1}} denotes the inverse link function of the specified family. The adjusted GEE predictor views the GEE as solving an iterative working linear model, and borrows ideas from universal kriging to exploit the *working cross-correlations* between the observed sites and the new sites. That is, for species i, we have
#'
#'\deqn{
#' \mu_{adj,i} = \mu_i + A_{i,new}^{1/2} R_{curr,new}^\top R_{curr,curr}^{-1} r_{i,curr},
#' }
#' where \eqn{r_{i,curr}} denote the Pearson residuals of species i, \eqn{R_{curr,curr}} is the working correlation matrix for current (training) data only, \eqn{R_{curr,new}} is the working cross-correlation matrix between current and new sites, and \eqn{A_{i,new}} is the diagonal matrix of marginal variances for species i at the new sites. Adjusted GEE predictions are available for all cases except `object$working_correlation_control$cor_type = "supplied"`. Of course if `cor_type = "independence"` then the cross-correlations are zero and the adjusted and standard GEE predictions coincide.
#'
#' Uncertainty (Wald-type) intervals for standard and adjusted GEE predictions are also provided, and we refer to Hui, Muller and Welsh (2025) for details on the construction of the adjusted GEE intervals. If an adjusted GEE prediction falls outside the valid range of the mean, then it is reverted to the standard GEE prediction. Moreover, for all uncertainty intervals, truncation is applied to the limits as appropriate to ensure they stay within the range of the response.
#'
#' Note the archetype assignments are treated as fixed, so uncertainty in the species groupings is not accounted for.
#'
#' @return
#' If `intervals = FALSE`, then a list containing the standard and adjusted GEE predictions, each of which is a matrix with `nrow(newdata)` rows and one column per species. If `intervals = TRUE`, then a list containing the standard and adjusted GEE predictions, along with their uncertainty intervals, each of which is again a matrix with `nrow(newdata)` rows and one column per species.
#'
#' @author
#' Francis K. C. Hui \email{fhui28@gmail.com}
#'
#' @seealso
#' [SAGEE()]
#'
#' @export
#' @md
predict.sagee <- function(object,
                          newdata,
                          new_spatial_coords,
                          intervals = FALSE,
                          intervals_quantiles = c(0.025, 0.975),
                          vcov_type = "model-based",
                          ...) {

    ##---------------------------------
    #' # Checks and set up
    ##---------------------------------
    if(!inherits(object, "sagee"))
        stop("object must be of class 'sagee'.")

    vcov_type <- match.arg(vcov_type, choices = c("model-based", "sandwich", "sandwich_adjusted"))
    cor_type <- object$working_correlation_control$cor_type
    if(cor_type == "supplied")
        stop("Adjusted predictions require the working cross-correlations between observed and new sites. This can not be computed when working_correlation_control$cor_type = \"supplied\"!")
    if(cor_type == "independence")
        message("The working correlation is independence, so the adjusted and standard GEE predictions coincide.")

    if(cor_type != "independence") {
        new_spatial_coords <- as.matrix(new_spatial_coords)
        if(!is.numeric(new_spatial_coords) || ncol(new_spatial_coords) != 2 || nrow(new_spatial_coords) != nrow(newdata))
            stop("new_spatial_coords must be a numeric matrix with two columns and the same number of rows as newdata.")
        spatial_coords_train <- NULL
        if(cor_type %in% c("matern", "NNGP"))
            spatial_coords_train <- as.matrix(object$working_correlation_control$spatial_coords)
        if(length(intervals_quantiles) != 2)
            stop("intervals_quantiles must be a vector of length 2.")
        }

    num_species <- ncol(object$species_coefficients)
    num_sites_new <- nrow(newdata)
    new_site_names <- rownames(newdata)
    p_alpha <- ncol(object$archetypal_coefficients)
    p_beta <- nrow(object$species_coefficients)
    species_names <- colnames(object$species_fitted_values)


    #' ## Construct model matrices for observed and new sites, along with marginal means and variances for new sites
    Xalpha_train <- model.matrix.sagee(object = object, which = "alpha")
    Xbeta_train <- model.matrix.sagee(object = object, which = "beta")
    Xalpha_new <- model.matrix.sagee(object = object, which = "alpha", newdata = newdata)
    Xbeta_new <- model.matrix.sagee(object = object, which = "beta", newdata = newdata)
    if(ncol(Xalpha_new) != p_alpha || ncol(Xbeta_new) != p_beta)
        stop("The model matrices constructed from newdata do not match the dimensions of the estimated coefficients.")

    eta_new <- tcrossprod(Xalpha_new, object$archetypal_coefficients)[, object$species_groupings, drop = FALSE] + Xbeta_new %*% object$species_coefficients
    mu_new <- object$family$linkinv(eta_new)

    if(object$family$family[1] %in% c("poisson", "binomial", "gaussian", "Gamma", "inverse.gaussian"))
        var_new <- object$family$variance(mu_new) * matrix(object$dispersion_parameters, nrow = num_sites_new, ncol = num_species, byrow = TRUE)
    if(object$family$family[1] == "negative_binomial")
        var_new <- object$family$variance(mu_new,
                                          phi = matrix(object$dispersion_parameters, nrow = num_sites_new, ncol = num_species, byrow = TRUE))
    if(object$family$family[1] == "tweedie")
        var_new <- object$family$variance(mu_new,
                                          phi = matrix(object$dispersion_parameters, nrow = num_sites_new, ncol = num_species, byrow = TRUE),
                                          power = matrix(object$power_parameters, nrow = num_sites_new, ncol = num_species, byrow = TRUE))

    get_alpha <- qnorm(intervals_quantiles)


    ##--------------------------------
    #' # Standard GEE predictions and their uncertainty intervals
    #' Wald interval on linear predictor and then apply inverse link to the limits
    ##--------------------------------
    if(intervals) {
        #' Function to compute the one covariance matrix of the coefficients for each species, i.e., (alpha_{g_i}, beta_i) in that order. Note this is indeed for both standard and adjusted predictions.
        .sagee_species_vcov <- function(object,
                                        vcov_type = "model-based") {
            num_species <- ncol(object$species_coefficients)
            alpha_names <- colnames(object$archetypal_coefficients)
            beta_names <- rownames(object$species_coefficients)
            species_names <- colnames(object$species_coefficients)
            species_groupings <- as.integer(object$species_groupings)

            V_model <- vcov.sagee(object, type = "model-based")
            if(vcov_type != "model-based")
                V_alpha <- vcov(object, type = vcov_type)

            out <- vector("list", num_species)
            names(out) <- species_names
            for(i in 1:num_species) {
                V <- V_model[[species_groupings[i]]]

                sel <- c(alpha_names, paste0(species_names[i], ":", beta_names)) # Extract the archetypal coefficients and the coefficients relevant to species i
                S <- as.matrix(V[sel, sel, drop = FALSE])
                #' Replace with sandwich estimates for alpha_g if necessary
                if(vcov_type != "model-based")
                    S[alpha_names, alpha_names] <- as.matrix(V_alpha[[species_groupings[i]]])[alpha_names, alpha_names]

                out[[i]] <- S
                }

            return(out)
            }

        species_vcov <- .sagee_species_vcov(object, vcov_type = vcov_type)

        intervals_standard_lower <- intervals_standard_upper <- matrix(NA, nrow = num_sites_new, ncol = num_species)
        for(i in 1:num_species) {
            se <- sqrt(rowSums((cbind(Xalpha_new, Xbeta_new) %*% species_vcov[[i]]) * cbind(Xalpha_new, Xbeta_new)))
            se <- pmax(se, 0)

            intervals_standard_lower[,i] <- object$family$linkinv(eta_new[,i] + get_alpha[1] * se)
            intervals_standard_upper[,i] <- object$family$linkinv(eta_new[,i] + get_alpha[2] * se)
            }
        rm(se)
        }

    ##--------------------------------
    #' # Adjusted GEE predictions, constructed for each distinct working correlation (one if common = TRUE, otherwise one per archetype)
    ##--------------------------------
    .build_cross_correlation <- function(train_coords,
                                         new_coords,
                                         params,
                                         smoothness) {

        out <- fields::stationary.cov(x1 = train_coords,
                                      x2 = new_coords,
                                      Covariance = "Matern",
                                      smoothness = smoothness,
                                      aRange = params[1])

        out <- params[2] * out + (1 - params[2]) * 0 # Nugget is only on the diagonal!
        return(out)
        }

    #' ## Construct kriging weights K = R_{train,train}^{-1} %*% R_{train, new} for each distinct working correlation: one if common = TRUE, otherwise one per archetype
    #' Also construct one of the bits relevant to computing uncertainty intervals for adjusted predictions, namely the diagonal of R_{curr,new}^T R_{curr,curr}^{-1} R_{curr,new}; see middle component of C_i matrix in the Biometrics paper...
    if(cor_type == "independence" || object$working_correlation_control$common) {
        map_correlation_index <- rep(1, num_species)
        correlation_params_list <- list(object$working_correlation_params)
        } else {
        map_correlation_index <- object$species_groupings
        correlation_params_list <- object$working_correlation_params
        }

    kriging_list <- lapply(1:length(correlation_params_list), function(x) {
        if(cor_type == "independence")
            return(list(krig_weights = matrix(0, nrow = nrow(Xalpha_train), ncol = num_sites_new),
                        interval_component = rep(0, num_sites_new)))

        if(object$working_correlation_control$common)
            R_traintrain <- object$working_correlation
        if(!object$working_correlation_control$common)
            R_traintrain <- object$working_correlation[[x]]

        R_trainnew <- .build_cross_correlation(train_coords = spatial_coords_train,
                                               new_coords = new_spatial_coords,
                                               params = correlation_params_list[[x]],
                                               smoothness = object$working_correlation_control$matern_smoothness)
        krig_weights <- base::solve(R_traintrain, R_trainnew)
        interval_component <- colSums(R_trainnew * krig_weights) # diag(R_cn^T R_cc^{-1} R_cn)

        return(list(krig_weights = krig_weights, # Vector of length num_sites_new
                    interval_component = interval_component)) # Vector of length num_sites_new
        })


    adjusted_mu_new <- matrix(NA, nrow = num_sites_new, ncol = num_species)
    if(intervals) {
        intervals_adjusted_lower <- intervals_adjusted_upper <- matrix(NA, nrow = num_sites_new, ncol = num_species)

        if(object$family$family[1] %in% c("poisson", "binomial", "gaussian", "Gamma", "inverse.gaussian"))
            var_train <- object$family$variance(object$species_fitted_values) * matrix(object$dispersion_parameters, nrow = nrow(Xalpha_train), ncol = num_species, byrow = TRUE)
        if(object$family$family[1] == "negative_binomial")
            var_train <- object$family$variance(object$species_fitted_values,
                                                phi = matrix(object$dispersion_parameters, nrow = nrow(Xalpha_train), ncol = num_species, byrow = TRUE))
        if(object$family$family[1] == "tweedie")
            var_train <- object$family$variance(object$species_fitted_values,
                                                phi = matrix(object$dispersion_parameters, nrow = nrow(Xalpha_train), ncol = num_species, byrow = TRUE),
                                                power = matrix(object$power_parameters, nrow = nrow(Xalpha_train), ncol = num_species, byrow = TRUE))
        }

    for(i in 1:num_species) {
        rootA_new <- sqrt(var_new[, i])
        make_adjusted_pred <- mu_new[, i] + rootA_new * as.vector(crossprod(kriging_list[[map_correlation_index[i]]]$krig_weights, object$pearson_residuals[, i]))
        check_issues <-  sapply(make_adjusted_pred, object$family$validmu)
        make_adjusted_pred[!check_issues] <- mu_new[,i][!check_issues]
        adjusted_mu_new[, i] <- make_adjusted_pred

        if(intervals) {
            T1 <- var_new[, i] # Diagonal of A_new^{1/2} R_{new,new} A_new^{1/2}, noting diag(R_new,new) = 1
            T2 <- var_new[, i] * kriging_list[[map_correlation_index[i]]]$interval_component # Diagonal of A_new^{1/2} R_{train,new}^T R_{train,train}^{-1} R_{train,new} A_new^{1/2}

            # Q = D_new - A_new^{1/2} R_cn^T R_cc^{-1} A_curr^{-1/2} D_curr, with D taken with respect to (alpha_{g_i}, beta_i)
            D_curr_scaled <- (object$family$mu.eta(object$species_linear_predictor[, i]) / sqrt(var_train[, i])) * cbind(Xalpha_train, Xbeta_train)
            Q <- object$family$mu.eta(eta_new[, i]) * cbind(Xalpha_new, Xbeta_new) - rootA_new * crossprod(kriging_list[[map_correlation_index[i]]]$krig_weights, D_curr_scaled)
            T3 <- rowSums((Q %*% species_vcov[[i]]) * Q) # Diagonal of Q S Q^T

            root_MSEP <- sqrt(pmax(T1 - T2 + T3, 0))
            UI <- cbind(as.vector(make_adjusted_pred + get_alpha[1] * root_MSEP),
                        as.vector(make_adjusted_pred + get_alpha[2] * root_MSEP))
            if(object$family$family[1] %in% c("Gamma", "poisson", "binomial", "inverse.gaussian", "negative_binomial", "tweedie"))
                UI[UI[,1] < 0, 1] <- 0
            if(object$family$family[1] == "binomial")
                UI[UI[,2] > 1, 2] <- 1

            intervals_adjusted_lower[, i] <- UI[,1]
            intervals_adjusted_upper[, i] <- UI[,2]

            rm(T1, T2, D_curr_scaled, Q, T3, root_MSEP, UI)
            }

        rm(rootA_new, make_adjusted_pred, check_issues)
        }


    ##--------------------------------
    #' # Format output
    ##--------------------------------
    out <- list(standard = mu_new, adjusted = adjusted_mu_new)
    if(intervals) {
        out <- list(standard = mu_new,
                    standard_lower = intervals_standard_lower,
                    standard_upper = intervals_standard_upper,
                    adjusted = adjusted_mu_new,
                    adjusted_lower = intervals_adjusted_lower,
                    adjusted_upper = intervals_adjusted_upper)
        }

    for(k in 1:length(out))
        dimnames(out[[k]]) <- list(new_site_names, species_names)

    return(out)
    }


