#' @title Construct model matrices from a fitted SAGEE model
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Constructs the model matrix for either the archetypal coefficients (from `Xalpha_formula`) or the species-specific coefficients (from `Xbeta_formula`) of a fitted `"sagee"` object.
#'
#' @param object An object of class `"sagee"`.
#' @param which Either `"alpha"` (the default) for the model matrix built from `Xalpha_formula` corresponding to the archetypal coefficients, or `"beta"` for the model matrix built from `Xbeta_formula` corresponding to the species-specific coefficients.
#' @param newdata An optional data frame containing the variables referenced by the chosen formula. If `NULL` (the default), the data used to fit the model are used.
#' @param ... Not used.
#'
#' @details
#' When `newdata` is supplied, the model matrix is constructed in the same way as for the fitted model, in the style of [stats::predict.lm()]. Any response on the left-hand side of `Xbeta_formula` is ignored.
#'
#' @return A numeric model matrix. When `which = "alpha"`, the number of columns correspond to the columns of `object$archetypal_coefficients`, and `which = "beta"` then the number of columns correspond to the rows of `object$species_coefficients`.
#'
#' @seealso
#' [SAGEE()]
#'
#' @author
#' Francis K. C. Hui \email{fhui28@gmail.com}
#'
#' @importFrom stats model.matrix model.frame terms delete.response .getXlevels na.pass
#'
#' @export
#' @md
model.matrix.sagee <- function(object,
                               which = c("alpha", "beta"),
                               newdata = NULL,
                               ...) {

    if(!inherits(object, "sagee"))
        stop("object must be of class 'sagee'.")
    which <- match.arg(which, choices = c("alpha", "beta"))

    form_formula <- if(which == "alpha") object$Xalpha_formula else object$Xbeta_formula
    expected_names <- if(which == "alpha") colnames(object$archetypal_coefficients) else rownames(object$species_coefficients)

    #' Model frame from the original data
    getterms <- stats::delete.response(stats::terms(form_formula, data = object$data))
    modframe_orig <- stats::model.frame(getterms, data = object$data)
    getterms_orig <- attr(modframe_orig, "terms")
    X_orig <- stats::model.matrix(getterms_orig, data = modframe_orig)

    if(is.null(newdata)) {
        X <- X_orig
        }
    if(!is.null(newdata)) {
        modframe_new <- stats::model.frame(getterms_orig,
                                           data = as.data.frame(newdata),
                                           xlev = stats::.getXlevels(getterms_orig, modframe_orig),
                                           na.action = stats::na.pass)
        X <- stats::model.matrix(getterms_orig,
                                 data = modframe_new,
                                 contrasts.arg = attr(X_orig, "contrasts"))
        }

    return(X)
    }
