#' Print a Fitted SAGEE Model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' 
#' Prints some basic information about a fitted species archetype generalized estimating equation (SAGEE) model.
#'
#' @param x A fitted SAGEE model object of class `"sagee"`, typically returned by [SAGEE()].
#' @param digits A positive integer specifying the number of significant digits used when printing numeric summaries.
#' @param ... Additional arguments. These are currently unused.
#'
#' @return The fitted model object `x`, returned invisibly.
#'
#' @author
#' Francis K. C. Hui \email{fhui28@gmail.com}
#' 
#' @seealso
#' [SAGEE()]
#'
#' @method print sagee
#' @export
print.sagee <- function(x,
                        digits = max(3L, getOption("digits") - 3L),
                        ...) {
    
    if(!is.numeric(digits) || length(digits) != 1 || !is.finite(digits) || digits < 1 || digits != as.integer(digits)) 
        stop("digits must be a positive integer.")
    digits <- as.integer(digits)
    
    ## Infer dimensions from stable fitted-object components.
    num_species <- length(x$species_groupings)
    num_sites <- nrow(x$data)
    num_archetypes <- nrow(x$archetypal_coefficients)

    family_name <- x$family$family[1L]
    link_name <- x$family$link[1L]
    correlation_type <- x$working_correlation_control$cor_type
    common_correlation <- x$working_correlation_control$common
    
    # cat("========================================================")
    # cat("\nSpecies Archetype Generalized Estimating Equation Model\n")
    # cat("========================================================\n\n")
    
    cat("Call:\n")
    print(x$call)
    
    cat("\nModel specification:\n")
    cat(" Family: ", x$family$family[1], "\n", sep = "")
    cat(" Link function: ", x$family$link[1], "\n", sep = "")
    # cat(" Number of sites: ", x$num_sites, "\n", sep = "")
    # cat(" Number of species: ", x$num_species, "\n", sep = "")
    cat(" Number of archetypes: ", nrow(x$archetypal_coefficients), "\n", sep = "")
    
    cat("\nWorking correlation:\n")
    if(!is.null(correlation_type)) {
        cat("  Structure: ", correlation_type, "\n", sep = "")
        if(correlation_type != "independence")
            cat("  Common: ", common_correlation, "\n", sep = "")
        } 

    invisible(x)
    }