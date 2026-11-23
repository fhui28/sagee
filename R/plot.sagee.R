#' @title Residual diagnostics for a fitted SAGEE model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Produces diagnostic plots based on either Pearson or Dunn-Smyth residuals of a fitted `"sagee"` object.
#'
#' @param x An object of class `"sagee"`.
#' @param type Either `"pearson"` (the default) or `"dunnsmyth"` residuals to use as the basis for plotting; see [residuals.sagee()] for details.
#' @param which_plot A vector of integers between 1 and 4. See details below for the plots available.
#' @param color_by Whether points should be colored by `"archetype"` (the default) or by `"species"`.
#' @param titles A character vector giving the titles of the four plots.
#' @param smooth If `TRUE`, then a smooth along with an approximate 95% confidence band is added to plots 1 and 3, courtesy of [mgcv::gam()].
#' @param envelope Should a simulation envelope be added to the normal probability plot (plot 2)?
#' @param envelope_color A character vector of length 2, giving the colors of the envelope boundaries and the envelope fill, respectively.
#' @param seed An optional seed for the simulation envelope, as well as for the Dunn-Smyth residuals if `type = "dunnsmyth"`.
#' @param ask If `TRUE`, the user is prompted before each new plot.
#' @param ... Additional graphical arguments passed to [graphics::plot()].
#'
#' @details
#' The four plots available are as follows:
#' \itemize{
#'   \item Residuals versus the linear predictors.
#'   \item A normal probability plot of the residuals.
#'   \item A scale-location plot, i.e., the square root of the absolute residuals versus the linear predictors.
#'   \item Boxplots of the residuals for each species, grouped by archetype.
#'   }
#'
#' Note that for discrete responses, Pearson residuals are far from normally distributed and form visible bands in plots 1 and 3, even when the assumed SAGEE model is correct. The normal probability plot should also be interpreted with a grain of salt for such response types.
#'
#' In plot 4, species whose boxplots are clearly shifted from zero and/or are much wider than others may be poorly described by their assigned archetype.
#'
#' @return The object `x`, invisibly.
#'
#' @seealso
#' [SAGEE()]
#'
#' @author
#' Francis K. C. Hui \email{fhui28@gmail.com}
#'
#' @importFrom graphics plot abline lines points polygon boxplot
#' @importFrom grDevices hcl.colors adjustcolor devAskNewPage
#' @importFrom stats qqnorm rnorm predict
#' @importFrom mgcv gam
#'
#' @export
#' @md
plot.sagee <- function(x,
                       type = "pearson",
                       which_plot = 1:4,
                       color_by = "archetype",
                       titles = c("Residuals vs. linear predictors",
                                  "Normal probability plot",
                                  "Scale-location plot",
                                  "Residuals vs. species"),
                       smooth = TRUE,
                       envelope = TRUE,
                       envelope_color = c("grey3", "lightgrey"),
                       seed = NULL,
                       ask = FALSE,
                       ...) {

    ##---------------------------------
    #' # Checks and set up
    ##---------------------------------
    if(!inherits(x, "sagee"))
        stop("x must be of class 'sagee'.")
    if(!is.numeric(which_plot) || any(!which_plot %in% 1:4))
        stop("which_plot must contain integers between 1 and 4.")
    if(length(titles) != 4)
        stop("titles must be a character vector of length 4.")
    color_by <- match.arg(color_by, choices = c("archetype", "species"))
    if(!is.null(seed))
        set.seed(seed)

    type <- match.arg(type, choices = c("pearson", "dunnsmyth"))
    getres <- residuals.sagee(object = x, type = type)
    linpred <- x$species_linear_predictor
    num_sites <- nrow(getres)
    num_species <- ncol(getres)
    species_groupings <- as.integer(x$species_groupings)
    num_archetypes <- length(x$archetype_sizes)
    residual_label <- if(type == "pearson") "Pearson residuals" else "Dunn-Smyth residuals"

    ok <- is.finite(getres)
    getres[!ok] <- NA

    archetype_cols <- grDevices::hcl.colors(num_archetypes, palette = "Dark3")
    if(color_by == "archetype") {
        species_cols <- archetype_cols[species_groupings]
        }
    if(color_by == "species") {
        species_cols <- grDevices::hcl.colors(num_species, palette = "Viridis")
        }
    point_cols <- matrix(grDevices::adjustcolor(species_cols, alpha.f = 0.5), nrow = num_sites, ncol = num_species, byrow = TRUE)

    # if(is.null(x$working_correlation_control$spatial_coords) && any(which_plot %in% 5)) {
    #     message("Plot 5 requires spatial coordinates stored in the fitted sagee object. It will be omitted from the set of constructed plots.")
    #     which_plot <- setdiff(which_plot, 5:6)
    #     }
    # spatial_coords <- x$working_correlation_control$spatial_coords
    # if(!is.null(spatial_coords))
    #     spatial_coords <- as.matrix(spatial_coords)

    if(ask) {
        old_ask <- grDevices::devAskNewPage(TRUE)
        on.exit(grDevices::devAskNewPage(old_ask), add = TRUE)
        }

    .sagee_add_smoother <- function(xx, yy, col = "grey3") {
        smooth_data <- data.frame(xx = xx, yy = yy)
        fit <- try(mgcv::gam(yy ~ s(xx), method = "REML", data = smooth_data), silent = TRUE)
        if(inherits(fit, "try-error"))
            return(invisible(NULL))

        grid <- data.frame(xx = seq(min(smooth_data$xx), max(smooth_data$xx), length.out = 200))
        pred <- stats::predict(fit, newdata = grid, se.fit = TRUE)
        graphics::lines(grid$xx, pred$fit, col = col, lwd = 2)
        graphics::lines(grid$xx, pred$fit + 2 * pred$se.fit, col = col, lty = 2)
        graphics::lines(grid$xx, pred$fit - 2 * pred$se.fit, col = col, lty = 2)
        return(invisible(NULL))
        }


    ##---------------------------------
    #' # Plot 1: Residuals vs. linear predictors
    ##---------------------------------
    if(1 %in% which_plot) {
        graphics::plot(x = linpred[ok],
                       y = getres[ok],
                       col = point_cols[ok],
                       pch = 16,
                       cex = 0.6,
                       xlab = "Linear predictors",
                       ylab = residual_label,
                       main = titles[1],
                       ...)
        graphics::abline(h = 0, lty = 2, col = "black")
        if(smooth)
            .sagee_add_smoother(xx = linpred[ok], yy = getres[ok])
        }


    ##---------------------------------
    #' # Plot 2: Quantile-quantile plot
    ##---------------------------------
    if(2 %in% which_plot) {
        getres_vec <- getres[ok]
        n <- length(getres_vec)
        ord <- order(getres_vec)
        theoretical <- stats::qqnorm(getres_vec, plot.it = FALSE)$x[ord]
        sorted_getres <- getres_vec[ord]

        ylim <- range(sorted_getres)
        if(envelope) {
            #' For large numbers of residuals, compute the envelope at a thinned set of order statistics, keeping both tails in full
            idx <- if(n > 2500) unique(round(c(1:200, seq(1, n, length = 2100), (n - 199):n))) else 1:n
            sims <- sapply(1:100, function(r) sort(stats::rnorm(n))[idx])
            env <- apply(sims, 1, stats::quantile, probs = c(0.025, 0.975))
            ylim <- range(ylim, env)
            }

        graphics::plot(x = theoretical,
                       y = sorted_getres,
                       type = "n",
                       ylim = ylim,
                       xlab = "Theoretical quantiles",
                       ylab = residual_label,
                       main = titles[2],
                       ...)
        if(envelope) {
            graphics::polygon(c(theoretical[idx], rev(theoretical[idx])),
                              c(env[1, ], rev(env[2, ])),
                              col = envelope_color[2], border = NA)
            graphics::lines(theoretical[idx], env[1, ], col = envelope_color[1], lty = 3)
            graphics::lines(theoretical[idx], env[2, ], col = envelope_color[1], lty = 3)
            }
        graphics::abline(0, 1, col = "black")
        graphics::points(x = theoretical,
                         y = sorted_getres,
                         col = point_cols[ok][ord],
                         pch = 16,
                         cex = 0.6)
        }


    ##---------------------------------
    #' # Plot 3: Scale-location plot
    ##---------------------------------
    if(3 %in% which_plot) {
        sqrt_abs_getres <- sqrt(abs(getres[ok]))
        graphics::plot(x = linpred[ok],
                       y = sqrt_abs_getres,
                       col = point_cols[ok],
                       pch = 16,
                       cex = 0.6,
                       xlab = "Linear predictors",
                       ylab = bquote(sqrt("|" * .(residual_label) * "|")),
                       main = titles[3], ...)
        if(smooth)
            .sagee_add_smoother(xx = linpred[ok], yy = sqrt_abs_getres)
        }


    ##---------------------------------
    #' # Plot 4: Residuals boxplots vs. species, grouped by archetype
    ##---------------------------------
    if(4 %in% which_plot) {
        species_order <- order(species_groupings, 1:num_species)
        getres_list <- lapply(species_order, function(j) getres[ok[, j], j])
        show_names <- num_species <= 50

        graphics::boxplot(getres_list,
                          col = species_cols[species_order],
                          names = if(show_names) colnames(getres)[species_order] else rep("", num_species),
                          las = 2,
                          cex.axis = if(num_species > 40) 0.5 else 0.7,
                          pch = 16,
                          cex = 0.4,
                          outcol = "lightgrey",
                          xlab = if(show_names) "" else "Species (ordered by archetype)",
                          ylab = residual_label,
                          main = titles[4])
        graphics::abline(h = 0, lty = 2, col = "black")
        }


    #' ##---------------------------------
    #' #' # Plot 5: Map of site-averaged residuals
    #' ##---------------------------------
    #' if(5 %in% which_plot) {
    #'     site_z <- rowMeans(getres, na.rm = TRUE) * sqrt(rowSums(ok))
    #'     site_cols <- grDevices::adjustcolor(ifelse(site_z > 0, "firebrick", "steelblue"), alpha.f = 0.6)
    #'     site_cex <- pmin(0.4 + abs(site_z) / 1.5, 3)
    #'
    #'     graphics::plot(spatial_coords,
    #'                    col = site_cols,
    #'                    pch = 16,
    #'                    cex = site_cex,
    #'                    asp = 1,
    #'                    xlab = if(is.null(colnames(spatial_coords))) "Coordinate 1" else colnames(spatial_coords)[1],
    #'                    ylab = if(is.null(colnames(spatial_coords))) "Coordinate 2" else colnames(spatial_coords)[2],
    #'                    main = titles[5],
    #'                    ...)
    #'     graphics::legend("topright", legend = c("z = +2", "z = -2"), pch = 16, pt.cex = 0.4 + 2 / 1.5,
    #'                      col = grDevices::adjustcolor(c("firebrick", "steelblue"), alpha.f = 0.6),
    #'                      bg = grDevices::adjustcolor("white", alpha.f = 0.7), cex = 0.7)
    #'     }


    return(invisible(x))
    }


