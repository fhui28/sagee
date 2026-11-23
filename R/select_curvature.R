#' @title Select the Number of Archetypes by the Simple Curvature Method
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Fits a sequence of species archetype generalized estimating equation (SAGEE) models and selects the number of archetypes using the simple curvature ratio
#'
#' \deqn{
#' C_g = \frac{J_{g-1} - J_g}{J_g - J_{g+1}},
#' }
#'
#' where \eqn{J_g = -2 \ell_g} and \eqn{\ell_g} denotes the SAGEE pseudo-log-likelihood for a model with \eqn{g} archetypes. The selected number of archetypes is the value of \eqn{g} with the largest curvature ratio.
#'
#' @param candidate_archetypes An increasing integer vector containing the candidate numbers of archetypes. At least three consecutive candidate values are required.
#' @param trace Logical value indicating whether progress messages should be printed.
#' @param ... Arguments passed to [SAGEE()]. The `num_archetypes` argument must *not* be supplied because it is set internally from `candidate_archetypes`.
#'
#' @details
#' See Zhang et al., (2017) for the (very simple!) idea and mathetical principles underlying the curvature method for selecting the number of clusters.
#'
#' @return A list containing the following elements (not necessarily in order):
#' \itemize{
#'   \item \code{call} The matched function call.
#'   \item \code{selected_num_archetypes} The selected number of archetypes.
#'   \item \code{selected_fit} The fitted SAGEE model corresponding to the selected number of archetypes.
#'   \item \code{selection_table} A data frame summarizing the results, including the number of archetypes, pseudo-log-likelihood, criterion, improvement, and curvature ratio.
#' }
#'
#'
#' @references
#' Zhang, Y., Mandziuk J, Quek C. H., and Goh, B.W. (2017).
#' Curvature-based method for determining the number of clusters
#' \emph{Information Sciences}.
#' \doi{https://doi.org/10.1016/j.ins.2017.05.024}
#'
#' @seealso
#' [SAGEE()], [summary.sagee()]
#'
#' @author
#' Francis K. C. Hui \email{fhui28@gmail.com}
#'
#'
#' @export
#' @md
select_curvature <- function(...,
                             candidate_archetypes,
                             trace = TRUE) {

    ##---------------------------------
    #' ## Validate candidate values
    ##---------------------------------
    if(!is.numeric(candidate_archetypes) || length(candidate_archetypes) < 3 || any(candidate_archetypes != as.integer(candidate_archetypes)))
        stop("candidate_archetypes must contain at least three ", "positive integers in consecutive order.")

    candidate_archetypes <- as.integer(candidate_archetypes)
    if(any(diff(candidate_archetypes) != 1))
        stop("candidate_archetypes must contain consecutive integers.")

    sagee_args <- list(...)
    if("num_archetypes" %in% names(sagee_args))
        stop("Do not supply num_archetypes through ...")


    ##---------------------------------
    #' ## Fit candidate models
    ##---------------------------------
    fits <- lapply(candidate_archetypes, function(g) {
        if(trace)
            message("Fitting SAGEE with ", g, " archetypes...")

        do.call(SAGEE, c(sagee_args, list(num_archetypes = g)))
        })

    names(fits) <- paste0("archetype", candidate_archetypes)

    #' ## Calculate the criterion
    criterion <- sapply(fits, function(x) { -2 * x$pseudo_logLik })
    improvement <- -diff(criterion)
    names(improvement) <- paste0(candidate_archetypes[-length(candidate_archetypes)], "_to_", candidate_archetypes[-1])
    curvature <- improvement[-length(improvement)] / improvement[-1]

    curvature_archetypes <- candidate_archetypes[2:(length(candidate_archetypes) - 1)]
    selected_num_archetypes <- curvature_archetypes[which.max(curvature)]
    selected_fit_index <- match(selected_num_archetypes, candidate_archetypes)
    selected_fit <- fits[[selected_fit_index]]

    ##---------------------------------
    #' # Wrap up and create summary table
    ##---------------------------------
    selection_table <- data.frame(num_archetypes = candidate_archetypes,
                                  criterion = criterion,
                                  improvement = c(NA, improvement),
                                  curvature = c(NA, curvature, NA),
                                  row.names = NULL)

    out <- list(call = match.call(),
                selected_num_archetypes = selected_num_archetypes,
                selected_fit = selected_fit,
                selection_table = selection_table)
    return(out)
    }
