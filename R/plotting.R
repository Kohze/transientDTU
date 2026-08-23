#' Plot stage-level evidence for a feature and focal group
#'
#' @param x A [TransientDTUResult-class] or stage table from
#'   [makeStageDTU()].
#' @param feature_id Feature identifier.
#' @param focal_group Focal biological group.
#' @param y Stage summary to plot: signed mean effect, minimum absolute effect,
#'   or maximum absolute effect.
#' @param show_threshold Add the signed effect threshold when available.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return The plotted data, invisibly.
#'
#' @examples
#' simulated <- simulateTransientDTU(n_genes = 12, seed = 13)
#' fit <- runTransientDTU(simulated$pairwise, simulated$stage_order)
#' plotTransientTrajectory(
#'     fit, simulated$truth$feature_id[[1]], simulated$groups[[1]]
#' )
#'
#' @export
plotTransientTrajectory <- function(
    x,
    feature_id,
    focal_group,
    y = c("mean_effect", "min_abs_effect", "max_abs_effect"),
    show_threshold = TRUE,
    ...
) {
    y <- match.arg(y)
    stage_data <- if (methods::is(x, "TransientDTUResult")) {
        as.data.frame(stageTable(x))
    } else {
        .assert_data_frame(x, "x")
    }
    .assert_columns(
        stage_data,
        c("feature_id", "focal_group", "stage", "stage_index", y),
        "stage data"
    )
    selected <- stage_data[
        as.character(stage_data$feature_id) == feature_id &
            as.character(stage_data$focal_group) == focal_group,
        , drop = FALSE
    ]
    if (!nrow(selected)) {
        stop("No stage data match 'feature_id' and 'focal_group'.",
            call. = FALSE)
    }
    selected <- selected[order(selected$stage_index), , drop = FALSE]
    graphics::plot(
        selected$stage_index, selected[[y]], type = "b", xaxt = "n",
        xlab = "Ordered stage", ylab = gsub("_", " ", y), ...
    )
    graphics::axis(
        1, at = selected$stage_index, labels = as.character(selected$stage)
    )
    graphics::abline(h = 0, col = "grey70", lty = 3)
    if (isTRUE(show_threshold) && methods::is(x, "TransientDTUResult")) {
        threshold <- decisionParameters(x)$effect_threshold
        if (!is.null(threshold) && is.finite(threshold)) {
            if (y == "mean_effect") {
                graphics::abline(h = c(-threshold, threshold),
                    col = "grey55", lty = 2)
            } else {
                graphics::abline(h = threshold, col = "grey55", lty = 2)
            }
        }
    }
    invisible(S4Vectors::DataFrame(selected))
}

#' Plot candidate selection frequency across a threshold grid
#'
#' @param stability Result from [thresholdStability()].
#' @param top Number of genes shown.
#' @param ... Additional arguments passed to [graphics::barplot()].
#'
#' @return The plotted selection-frequency table, invisibly.
#'
#' @examples
#' simulated <- simulateTransientDTU(n_genes = 12, seed = 14)
#' stage <- makeStageDTU(simulated$pairwise, simulated$stage_order)
#' stability <- thresholdStability(
#'     stage, q_threshold = c(0.01, 0.05),
#'     effect_threshold = c(0.1, 0.15)
#' )
#' plotThresholdStability(stability)
#'
#' @export
plotThresholdStability <- function(stability, top = 20L, ...) {
    .assert_count(top, "top")
    if (!is.list(stability) || is.null(stability$selection_frequency)) {
        stop("'stability' must be returned by thresholdStability().",
            call. = FALSE)
    }
    frequency <- as.data.frame(stability$selection_frequency)
    if (!nrow(frequency)) {
        stop("No genes were selected in the threshold grid.", call. = FALSE)
    }
    frequency <- utils::head(frequency, top)
    graphics::barplot(
        height = frequency$selection_frequency,
        names.arg = frequency$gene_id,
        las = 2,
        ylab = "Selection frequency",
        ...
    )
    invisible(S4Vectors::DataFrame(frequency))
}
