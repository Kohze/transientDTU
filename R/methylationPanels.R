#' Simulate methylation and transcript-fraction measurements
#'
#' Generates a compact, reproducible example with developmental stages,
#' biological replicates, two complementary transcript fractions, and
#' CpG-level methylation measurements. The generator is intended for plotting
#' examples and interface testing; it is not a sequencing or causal model.
#'
#' @param n_stages Number of ordered developmental stages; at least three.
#' @param n_replicates Replicates per stage; at least two.
#' @param n_cpg Number of CpG positions; at least two.
#' @param stage_labels Optional unique labels of length `n_stages`.
#' @param transcripts Two unique transcript identifiers.
#' @param noise_sd Gaussian observation noise on the proportion scale.
#' @param seed Integer random seed.
#'
#' @return A list containing `observations`, `cpg`, `stage_order`, and
#'   `simulation_parameters`.
#'
#' @details
#' Methylation and transcript fractions are generated on the unit interval.
#' The first transcript is positively associated with methylation and the
#' second is its compositional complement. This known relationship makes the
#' output useful for checking descriptive plots, not for evaluating biological
#' or causal inference.
#'
#' @examples
#' methylation_example <- simulateMethylationExpression(seed = 22)
#' head(methylation_example$observations)
#' head(methylation_example$cpg)
#'
#' @export
simulateMethylationExpression <- function(
    n_stages = 8L,
    n_replicates = 2L,
    n_cpg = 24L,
    stage_labels = NULL,
    transcripts = c("transcript_1", "transcript_2"),
    noise_sd = 0.025,
    seed = 1L
) {
    .assert_count(n_stages, "n_stages", 3L)
    .assert_count(n_replicates, "n_replicates", 2L)
    .assert_count(n_cpg, "n_cpg", 2L)
    .assert_number(noise_sd, "noise_sd", .Machine$double.eps, 0.25)
    if (length(seed) != 1L || !is.numeric(seed) || is.na(seed) ||
            seed != as.integer(seed)) {
        stop("'seed' must be one integer.", call. = FALSE)
    }
    if (is.null(stage_labels)) {
        stage_labels <- paste0("stage", seq_len(n_stages))
    }
    stage_labels <- .as_id(stage_labels, "stage_labels")
    if (length(stage_labels) != n_stages || anyDuplicated(stage_labels)) {
        stop(
            "'stage_labels' must contain one unique label per stage.",
            call. = FALSE
        )
    }
    transcripts <- .as_id(transcripts, "transcripts")
    if (length(transcripts) != 2L || anyDuplicated(transcripts)) {
        stop("'transcripts' must contain two unique identifiers.",
            call. = FALSE)
    }

    withr::with_seed(seed, {
        samples <- expand.grid(
            replicate = seq_len(n_replicates),
            stage_index = seq_len(n_stages),
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
        samples <- samples[order(
            samples$stage_index, samples$replicate
        ), , drop = FALSE]
        samples$stage <- stage_labels[samples$stage_index]
        samples$sample_id <- paste0(
            samples$stage, "_rep", samples$replicate
        )
        developmental_time <- (samples$stage_index - 1) / (n_stages - 1)
        methylation_mean <- 0.74 - 0.40 * developmental_time +
            0.05 * sin(2 * pi * developmental_time)
        samples$methylation <- .clamp_unit(
            methylation_mean + stats::rnorm(nrow(samples), 0, noise_sd)
        )
        associated_fraction <- 0.18 +
            0.64 * .rescale_unit(samples$methylation) +
            stats::rnorm(nrow(samples), 0, noise_sd)
        associated_fraction <- .clamp_unit(associated_fraction)

        observations <- rbind(
            data.frame(
                samples,
                transcript_id = transcripts[[1L]],
                transcript_fraction = associated_fraction,
                stringsAsFactors = FALSE
            ),
            data.frame(
                samples,
                transcript_id = transcripts[[2L]],
                transcript_fraction = 1 - associated_fraction,
                stringsAsFactors = FALSE
            )
        )
        observations <- observations[order(
            observations$stage_index,
            observations$replicate,
            observations$transcript_id
        ), , drop = FALSE]
        rownames(observations) <- NULL

        positions <- 1000L + 25L * (seq_len(n_cpg) - 1L)
        cpg <- merge(
            samples,
            data.frame(
                cpg_index = seq_len(n_cpg),
                position = positions
            ),
            by = NULL
        )
        position_effect <- stats::rnorm(n_cpg, 0, 0.08)
        cpg$methylation <- .clamp_unit(
            cpg$methylation + position_effect[cpg$cpg_index] +
                stats::rnorm(nrow(cpg), 0, noise_sd)
        )
        cpg <- cpg[order(
            cpg$stage_index, cpg$replicate, cpg$cpg_index
        ), , drop = FALSE]
        rownames(cpg) <- NULL

        list(
            observations = observations,
            cpg = cpg,
            stage_order = stage_labels,
            simulation_parameters = list(
                n_stages = n_stages,
                n_replicates = n_replicates,
                n_cpg = n_cpg,
                noise_sd = noise_sd,
                seed = as.integer(seed)
            )
        )
    })
}

#' Plot methylation and transcript measurements across development
#'
#' Produces a descriptive multi-panel view modeled on the companion thesis
#' analysis: stage-mean methylation versus transcript fraction, scaled
#' developmental trajectories, and an optional CpG-by-stage methylation map.
#'
#' @param observations A data frame containing methylation and transcript
#'   measurements in long format.
#' @param cpg Optional data frame containing CpG-level methylation measurements.
#' @param stage_order Optional complete stage order. Factor levels or first
#'   appearance are used when omitted.
#' @param stage_col,transcript_col,methylation_col,value_col Column names in
#'   `observations`.
#' @param cpg_stage_col,cpg_position_col,cpg_methylation_col Column names in
#'   `cpg`.
#' @param main Optional overall figure title.
#' @param colors Optional transcript colors.
#' @param show_fits Add transcript-specific least-squares trend lines.
#'
#' @return Invisibly, a list of plotting summaries: `stage_summary`,
#'   `correlations`, `trajectories`, and `cpg_summary`.
#'
#' @details
#' The function is descriptive. Correlations and fitted lines do not establish
#' methylation causality, and input measurements should already be normalized
#' and matched to the intended biological units. Replicates are averaged within
#' stage for display; upstream uncertainty and multiplicity control remain the
#' responsibility of the analysis supplying the data.
#'
#' @examples
#' methylation_example <- simulateMethylationExpression(seed = 22)
#' panel <- plotMethylationExpression(
#'     methylation_example$observations,
#'     cpg = methylation_example$cpg,
#'     stage_order = methylation_example$stage_order,
#'     main = "Example methylation and transcript-fraction panel"
#' )
#' panel$correlations
#'
#' @export
plotMethylationExpression <- function(
    observations,
    cpg = NULL,
    stage_order = NULL,
    stage_col = "stage",
    transcript_col = "transcript_id",
    methylation_col = "methylation",
    value_col = "transcript_fraction",
    cpg_stage_col = "stage",
    cpg_position_col = "position",
    cpg_methylation_col = "methylation",
    main = NULL,
    colors = NULL,
    show_fits = TRUE
) {
    observations <- .assert_data_frame(observations, "observations")
    .assert_columns(
        observations,
        c(stage_col, transcript_col, methylation_col, value_col),
        "observations"
    )
    display_data <- data.frame(
        stage = .as_id(observations[[stage_col]], stage_col),
        transcript_id = .as_id(
            observations[[transcript_col]], transcript_col
        ),
        methylation = observations[[methylation_col]],
        transcript_value = observations[[value_col]],
        stringsAsFactors = FALSE
    )
    .assert_finite_numeric(
        display_data$methylation, methylation_col, "observations"
    )
    .assert_finite_numeric(
        display_data$transcript_value, value_col, "observations"
    )
    if (!nrow(display_data)) {
        stop("'observations' must contain at least one row.", call. = FALSE)
    }

    stage_order <- .resolve_panel_stage_order(
        observations[[stage_col]], stage_order
    )
    if (any(!display_data$stage %in% stage_order)) {
        stop("'stage_order' omits observed stages.", call. = FALSE)
    }
    transcripts <- unique(display_data$transcript_id)
    if (length(transcripts) < 1L) {
        stop("At least one transcript is required.", call. = FALSE)
    }

    stage_summary <- stats::aggregate(
        display_data[c("methylation", "transcript_value")],
        by = display_data[c("stage", "transcript_id")],
        FUN = mean
    )
    stage_summary$stage_index <- match(stage_summary$stage, stage_order)
    stage_summary <- stage_summary[order(
        stage_summary$stage_index, stage_summary$transcript_id
    ), , drop = FALSE]

    correlations <- do.call(rbind, lapply(transcripts, function(transcript) {
        selected <- stage_summary[
            stage_summary$transcript_id == transcript, , drop = FALSE
        ]
        informative <- nrow(selected) >= 3L &&
            length(unique(selected$methylation)) >= 2L &&
            length(unique(selected$transcript_value)) >= 2L
        correlation <- if (informative) {
            stats::cor(selected$methylation, selected$transcript_value)
        } else {
            NA_real_
        }
        data.frame(
            transcript_id = transcript,
            n_stages = nrow(selected),
            correlation = correlation,
            r_squared = correlation^2,
            stringsAsFactors = FALSE
        )
    }))
    rownames(correlations) <- NULL

    methylation_trajectory <- stats::aggregate(
        display_data["methylation"],
        by = display_data["stage"],
        FUN = mean
    )
    names(methylation_trajectory)[[2L]] <- "value"
    methylation_trajectory$series <- "methylation"
    transcript_trajectory <- stage_summary[
        c("stage", "transcript_id", "transcript_value")
    ]
    names(transcript_trajectory) <- c("stage", "series", "value")
    trajectories <- rbind(
        methylation_trajectory[c("stage", "series", "value")],
        transcript_trajectory
    )
    trajectories$stage_index <- match(trajectories$stage, stage_order)
    trajectories$scaled_value <- stats::ave(
        trajectories$value,
        trajectories$series,
        FUN = .rescale_unit
    )
    trajectories <- trajectories[order(
        match(trajectories$series, c("methylation", transcripts)),
        trajectories$stage_index
    ), , drop = FALSE]

    cpg_summary <- NULL
    if (!is.null(cpg)) {
        cpg <- .assert_data_frame(cpg, "cpg")
        .assert_columns(
            cpg,
            c(cpg_stage_col, cpg_position_col, cpg_methylation_col),
            "cpg"
        )
        cpg_values <- data.frame(
            stage = .as_id(cpg[[cpg_stage_col]], cpg_stage_col),
            position = cpg[[cpg_position_col]],
            methylation = cpg[[cpg_methylation_col]],
            stringsAsFactors = FALSE
        )
        .assert_finite_numeric(
            cpg_values$methylation, cpg_methylation_col, "cpg"
        )
        if (!nrow(cpg_values) || any(!cpg_values$stage %in% stage_order)) {
            stop("'cpg' must contain measurements at declared stages.",
                call. = FALSE)
        }
        if (anyNA(cpg_values$position)) {
            stop("CpG positions cannot be missing.", call. = FALSE)
        }
        cpg_summary <- stats::aggregate(
            cpg_values["methylation"],
            by = cpg_values[c("stage", "position")],
            FUN = mean
        )
        cpg_summary$stage_index <- match(cpg_summary$stage, stage_order)
    }

    if (is.null(colors)) {
        colors <- grDevices::hcl.colors(length(transcripts), "Dark 3")
    }
    if (length(colors) < length(transcripts)) {
        stop("'colors' must provide one color per transcript.",
            call. = FALSE)
    }
    colors <- colors[seq_along(transcripts)]
    names(colors) <- transcripts

    old_parameters <- graphics::par(no.readonly = TRUE)
    on.exit({
        graphics::layout(1L)
        graphics::par(old_parameters)
    }, add = TRUE)
    has_cpg <- !is.null(cpg_summary)
    panel_count <- if (has_cpg) 3L else 2L
    graphics::layout(
        matrix(seq_len(panel_count), ncol = 1L),
        heights = if (has_cpg) c(1, 1, 1.2) else c(1, 1)
    )
    graphics::par(mar = c(4.1, 4.3, 2.5, 1.0), oma = c(0, 0, 2.2, 0))

    graphics::plot(
        NA_real_, NA_real_,
        xlim = .expand_panel_range(stage_summary$methylation),
        ylim = .expand_panel_range(stage_summary$transcript_value),
        xlab = "Mean methylation",
        ylab = "Mean transcript fraction",
        main = "Methylation-transcript association"
    )
    legend_labels <- character(length(transcripts))
    for (index in seq_along(transcripts)) {
        transcript <- transcripts[[index]]
        selected <- stage_summary[
            stage_summary$transcript_id == transcript, , drop = FALSE
        ]
        graphics::points(
            selected$methylation, selected$transcript_value,
            pch = 19, col = colors[[transcript]]
        )
        informative <- nrow(selected) >= 2L &&
            length(unique(selected$methylation)) >= 2L
        if (isTRUE(show_fits) && informative) {
            graphics::abline(
                stats::lm(transcript_value ~ methylation, data = selected),
                col = colors[[transcript]], lwd = 2
            )
        }
        r_squared <- correlations$r_squared[
            correlations$transcript_id == transcript
        ]
        legend_labels[[index]] <- if (is.finite(r_squared)) {
            sprintf("%s (R2 = %.2f)", transcript, r_squared)
        } else {
            transcript
        }
    }
    graphics::legend(
        "topright", legend = legend_labels, col = colors,
        pch = 19, lty = if (isTRUE(show_fits)) 1 else NA,
        bty = "n", cex = 0.8
    )

    series_order <- c("methylation", transcripts)
    trajectory_colors <- c(methylation = "grey25", colors)
    trajectory_types <- c(2L, rep(1L, length(transcripts)))
    graphics::plot(
        NA_real_, NA_real_,
        xlim = c(1, length(stage_order)), ylim = c(0, 1), xaxt = "n",
        xlab = "Developmental stage", ylab = "Scaled mean value",
        main = "Developmental trajectories"
    )
    graphics::axis(1, at = seq_along(stage_order), labels = stage_order)
    for (index in seq_along(series_order)) {
        series <- series_order[[index]]
        selected <- trajectories[
            trajectories$series == series, , drop = FALSE
        ]
        selected <- selected[order(selected$stage_index), , drop = FALSE]
        graphics::lines(
            selected$stage_index, selected$scaled_value,
            type = "b", pch = 19, col = trajectory_colors[[series]],
            lty = trajectory_types[[index]], lwd = 2
        )
    }
    graphics::legend(
        "topright", legend = series_order,
        col = trajectory_colors[series_order],
        lty = trajectory_types, pch = 19, bty = "n", cex = 0.8
    )

    if (has_cpg) {
        positions <- sort(unique(cpg_summary$position))
        methylation_matrix <- matrix(
            NA_real_, nrow = length(positions), ncol = length(stage_order),
            dimnames = list(as.character(positions), stage_order)
        )
        row_index <- match(cpg_summary$position, positions)
        column_index <- match(cpg_summary$stage, stage_order)
        methylation_matrix[cbind(row_index, column_index)] <-
            cpg_summary$methylation
        palette <- grDevices::colorRampPalette(
            c("#2166AC", "#F7F7F7", "#B2182B")
        )(100L)
        finite_values <- methylation_matrix[is.finite(methylation_matrix)]
        color_range <- range(finite_values)
        if (diff(color_range) <= sqrt(.Machine$double.eps)) {
            color_range <- color_range + c(-0.5, 0.5)
        }
        graphics::image(
            x = seq_along(stage_order),
            y = seq_along(positions),
            z = t(methylation_matrix),
            col = palette, zlim = color_range,
            xaxt = "n", yaxt = "n",
            xlab = "Developmental stage", ylab = "CpG position",
            main = "CpG methylation states"
        )
        graphics::axis(1, at = seq_along(stage_order), labels = stage_order)
        position_ticks <- unique(round(seq(
            1, length(positions), length.out = min(6L, length(positions))
        )))
        graphics::axis(
            2, at = position_ticks, labels = positions[position_ticks],
            las = 1, cex.axis = 0.75
        )
    }
    if (!is.null(main) && length(main) == 1L && !is.na(main) && nzchar(main)) {
        graphics::mtext(main, side = 3, outer = TRUE, line = 0.5,
            font = 2)
    }

    invisible(list(
        stage_summary = S4Vectors::DataFrame(stage_summary),
        correlations = S4Vectors::DataFrame(correlations),
        trajectories = S4Vectors::DataFrame(trajectories),
        cpg_summary = if (is.null(cpg_summary)) NULL else
            S4Vectors::DataFrame(cpg_summary)
    ))
}

.clamp_unit <- function(x) {
    pmin(0.995, pmax(0.005, x))
}

.rescale_unit <- function(x) {
    value_range <- range(x, na.rm = TRUE)
    if (!all(is.finite(value_range)) ||
            diff(value_range) <= sqrt(.Machine$double.eps)) {
        return(rep(0.5, length(x)))
    }
    (x - value_range[[1L]]) / diff(value_range)
}

.assert_finite_numeric <- function(x, column, data_name) {
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
        stop(
            "'", column, "' in '", data_name,
            "' must contain only finite numeric values.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.resolve_panel_stage_order <- function(stage, stage_order) {
    observed <- as.character(stage)
    if (is.null(stage_order)) {
        stage_order <- if (is.factor(stage)) levels(stage) else unique(observed)
    }
    stage_order <- .as_id(stage_order, "stage_order")
    if (anyDuplicated(stage_order)) {
        stop("'stage_order' must contain unique labels.", call. = FALSE)
    }
    stage_order
}

.expand_panel_range <- function(x) {
    value_range <- range(x, na.rm = TRUE)
    if (diff(value_range) <= sqrt(.Machine$double.eps)) {
        padding <- max(0.04, abs(value_range[[1L]]) * 0.04)
    } else {
        padding <- diff(value_range) * 0.04
    }
    value_range + c(-padding, padding)
}
