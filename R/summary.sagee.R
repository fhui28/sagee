#' @title Summarize a fitted SAGEE model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Summarizes a fitted species archetype generalized estimating equation (SAGEE) model, including performing Wald-based tests for the archetypal regression coefficients.
#'
#' @param object An object of class `"sagee"`.
#' @param se_type Character string specifying the standard errors used for inference on the archetypal coefficients. One of `"model-based"`, `"sandwich"`, or `"sandwich_adjusted"`, corresponding to model-based, sandwich, and sandwich-adjusted (MacKinnon, 1985) standard errors of the archetypal coefficients, respectively.
#' @param conf_level Numeric value between 0 and 1 specifying the confidence level for the Wald-based tests. Default is 0.95.
#' @param ... Additional arguments. Currently unused.
#'
#' @details
#' The Wald-type hypothesis tests are conditional on the estimated species memberships, and do not account for uncertainty arising from these archetype assignments.
#'
#' @return An object of class `"summary.sagee"` containing model information, archetypal coefficient tables, and working correlation summaries as appropriate.
#'
#' @references
#' MacKinnon, J. (1985).
#' Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties.
#' \emph{Journal of Econometrics}.
#' \doi{10.1016/0304-4076(85)90158-7}
#'
#' @author
#' Francis K. C. Hui \email{fhui28@gmail.com}
#'
#' @importFrom stats pnorm
#'
#' @method summary sagee
#' @export
summary.sagee <- function(object,
                          se_type = c("sandwich_adjusted", "sandwich", "model-based"),
                          conf_level = 0.95,
                          ...) {

    if(!inherits(object, "sagee"))
        stop("object must inherit from class 'sagee'.")

    se_type <- match.arg(se_type)
    if(!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1)
        stop("conf_level must be a scalar between zero and one.")

    coefficients <- object$archetypal_coefficients
    #num_terms <- ncol(coefficients)
    se_names <- dimnames(object$se_archetypal_coefficients)[[3]]

    selected_se <- object$se_archetypal_coefficients[, , match(se_type, se_names), drop = FALSE]
    selected_se <- matrix(selected_se, nrow = nrow(coefficients), ncol = ncol(coefficients), dimnames = dimnames(coefficients))
    #critical_value <- stats::qnorm(1 - (1 - conf_level) / 2)
    statistic <- coefficients / selected_se
    p_value <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)

    coefficient_tables <- vector("list", length = nrow(coefficients))
    names(coefficient_tables) <- rownames(coefficients)
    for(k0 in 1:nrow(coefficients)) {
        coefficient_tables[[k0]] <- data.frame(term = colnames(coefficients),
                                               estimate = unname(coefficients[k0, ]),
                                               std_error = unname(selected_se[k0, ]),
                                               statistic = unname(statistic[k0, ]),
                                               p_value = unname(p_value[k0, ]),
                                               row.names = colnames(coefficients))
        }

    #species_groupings <- object$species_groupings

    result <- list(call = object$call,
                   family = object$family$family[1],
                   link = object$family$link[1],
                   num_sites = nrow(object$data),
                   num_species = length(object$species_groupings),
                   num_archetypes = nrow(coefficients),
                   archetype_sizes = object$archetype_sizes,
                   Xalpha_formula = object$Xalpha_formula,
                   Xbeta_formula = object$Xbeta_formula,
                   correlation_type = object$working_correlation_control$cor_type,
                   working_correlation_params = object$working_correlation_params,
                   quasi_logLik = object$quasi_logLik,
                   se_type = se_type,
                   #conf_level = conf_level,
                   coefficients = coefficients,
                   standard_errors = selected_se,
                   coefficient_tables = coefficient_tables)

    class(result) <- "summary.sagee"
    return(result)
    }


#' @title Print a SAGEE model summary
#'
#' @description
#' Prints some model information and Wald-based tests for the archetypal regression coefficients from a summarized SAGEE model fit.
#'
#' @param x An object of class `"summary.sagee"`.
#' @param digits Number of significant digits used for printing.
#' @param signif_stars Logical value indicating whether significance stars should be printed.
#' @param ... Additional arguments. Currently unused.
#'
#' @return `x`, invisibly.
#'
#' @export
print.summary.sagee <- function(x,
                                digits = max(3L, getOption("digits") - 3L),
                                signif_stars = getOption("show.signif.stars"),
                                ...) {

    if(!inherits(x, "summary.sagee"))
        stop("x must inherit from class 'summary.sagee'.")


    cat("========================================================")
    cat("\nSpecies Archetype Generalized Estimating Equation Model\n")
    cat("========================================================\n\n")

    cat("Call:\n")
    print(x$call)

    cat("\nModel specification:\n")
    cat(" Family: ", x$family, "\n", sep = "")
    cat(" Link function: ", x$link, "\n", sep = "")
    cat(" Number of sites: ", x$num_sites, "\n", sep = "")
    cat(" Number of species: ", x$num_species, "\n", sep = "")
    cat(" Number of archetypes: ", x$num_archetypes, "\n", sep = "")
    cat("Quasi-log-likelihood: ", formatC(x$quasi_logLik, digits = digits, format = "g"), "\n", sep = "")

    cat("\nSpecies allocation:\n")
    allocation_table <- data.frame(Archetype = names(x$archetype_sizes),
                                   `Number of species` = as.integer(x$archetype_sizes))
    print(allocation_table,
          row.names = FALSE)

    cat("\nWorking correlation:\n")
    cat("  Structure: ", x$correlation_type, "\n", sep = "")
    if(!is.null(x$working_correlation_params)) {
        cat("  Estimated parameters:\n")
        print(x$working_correlation_params, digits = digits)
        }

    #confidence_percent <- format(100 * x$conf_level, trim = TRUE, scientific = FALSE)

    cat("\nArchetypal coefficient tests (using ", x$se_type, " standard errors):\n", sep = "")
    coefficient_table <- do.call(rbind,
                                 lapply(names(x$coefficient_tables), function(archetype) {
                                     table_g <- x$coefficient_tables[[archetype]]
                                     data.frame(Archetype = archetype,
                                                Term = table_g$term,
                                                Estimate = table_g$estimate,
                                                `Std. Error` = table_g$std_error,
                                                `z value` = table_g$statistic,
                                                `Pr(>|z|)` = table_g$p_value,
                                                check.names = FALSE,
                                                row.names = NULL)
                                    })
                                 )

    if(isTRUE(signif_stars))
        coefficient_table[[" "]] <- as.character(cut(coefficient_table[["Pr(>|z|)"]],
                                                           breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                                                           labels = c("***", "**", "*", ".", " ")))

    numeric_columns <- sapply(coefficient_table, is.numeric)
    coefficient_table[numeric_columns] <- lapply(coefficient_table[numeric_columns],
                                                 function(column) { signif(column, digits = digits) })

    cat(" ---\n")
    print(coefficient_table,
          row.names = FALSE,
          right = TRUE)
    cat(" ---")

    if(isTRUE(signif_stars))
        cat("\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n", sep = "")

    invisible(x)
    }


