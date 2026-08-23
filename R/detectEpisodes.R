#' Detect bounded diverge-reconverge episodes
#'
#' Applies a frozen post-inference decision rule to standardized stage-level
#' evidence from [makeStageDTU()]. Consecutive eligible stages are collapsed
#' into a single episode and retained only when the requested flanks satisfy
#' the reconvergence tolerance.
#'
#' @param stage_data Output from [makeStageDTU()] or an equivalent table with
#'   the documented columns.
#' @param stage_order Optional complete stage order. If omitted, the
#'   `stage_order` attribute created by [makeStageDTU()] is used.
#' @param q_threshold Strict upper bound for every required component adjusted
#'   p-value.
#' @param effect_threshold Inclusive lower bound for every absolute focal versus
#'   comparator effect.
#' @param comparator_tolerance Strict upper bound for the range of comparator
#'   means at an episode stage.
#' @param flank_tolerance Strict upper bound for every absolute focal versus
#'   comparator effect at both immediate flanks.
#' @param gene_q_threshold Optional strict upper bound for gene-level adjusted
#'   p-values. Use `NULL` to disable gene-level screening.
#' @param min_comparators Minimum number of comparator groups required.
#' @param min_episode_stages,max_episode_stages Inclusive bounds on the number
#'   of consecutive eligible stages in an episode. `Inf` disables the upper
#'   bound.
#' @param flank_width Number of observed stages required before and after an
#'   episode. The paper-compatible default is one on each side.
#' @param stage_coordinates Optional strictly increasing numeric coordinates in
#'   `stage_order`, or a named numeric vector keyed by stage. These annotate
#'   irregularly spaced designs without changing adjacency or threshold rules.
#' @param flank_missing `"available"` reproduces the original rule by using
#'   all finite requested-flank effects and requiring at least one finite
#'   flank. `"complete"` requires finite summaries at every requested flank.
#'   Entirely absent stage rows always break an episode under either setting.
#'
#' @return A [S4Vectors::DataFrame] with one row per detected feature-level
#'   episode. Empty analyses return the same typed columns with zero rows.
#'
#' @details
#' Missing stage rows are inserted as ineligible states and therefore break a
#' run. Runs touching the first or last stage are excluded because two observed
#' reconverged flanks cannot be demonstrated. Threshold strictness deliberately
#' matches the published decision rule: q-values, comparator spread, and flank
#' differences use `<`; episode effects use `>=`.
#'
#' @examples
#' simulated <- simulateTransientDTU(n_genes = 12, seed = 4)
#' stage <- makeStageDTU(simulated$pairwise, simulated$stage_order)
#' episodes <- detectEpisodes(stage)
#' episodes
#'
#' @export
detectEpisodes <- function(
    stage_data,
    stage_order = NULL,
    q_threshold = 0.05,
    effect_threshold = 0.10,
    comparator_tolerance = effect_threshold,
    flank_tolerance = effect_threshold,
    gene_q_threshold = NULL,
    min_comparators = 2L,
    min_episode_stages = 1L,
    max_episode_stages = Inf,
    flank_width = 1L,
    stage_coordinates = NULL,
    flank_missing = c("available", "complete")
) {
    stage_attribute <- attr(stage_data, "stage_order", exact = TRUE)
    expected_groups <- attr(stage_data, "expected_groups", exact = TRUE)
    stage_data <- .assert_data_frame(stage_data, "stage_data")
    .assert_number(q_threshold, "q_threshold", 0, 1)
    .assert_number(effect_threshold, "effect_threshold", 0, Inf)
    .assert_number(
        comparator_tolerance, "comparator_tolerance", 0, Inf
    )
    .assert_number(flank_tolerance, "flank_tolerance", 0, Inf)
    .assert_count(min_comparators, "min_comparators")
    .assert_count(min_episode_stages, "min_episode_stages")
    .assert_count_or_inf(max_episode_stages, "max_episode_stages")
    .assert_count(flank_width, "flank_width")
    if (max_episode_stages < min_episode_stages) {
        stop(
            "'max_episode_stages' cannot be smaller than ",
            "'min_episode_stages'.",
            call. = FALSE
        )
    }
    flank_missing <- match.arg(flank_missing)
    if (!is.null(gene_q_threshold)) {
        .assert_number(gene_q_threshold, "gene_q_threshold", 0, 1)
    }

    required <- c(
        "feature_id", "gene_id", "gene_name", "focal_group", "stage",
        "stage_index", "n_comparators", "complete_comparisons",
        "direction_coherent", "direction_sign", "mean_effect",
        "min_abs_effect", "max_abs_effect", "comparator_spread", "worst_q",
        "gene_q"
    )
    .assert_columns(stage_data, required, "stage_data")

    if (is.null(stage_order)) stage_order <- stage_attribute
    if (is.null(stage_order)) {
        indexed <- unique(stage_data[c("stage", "stage_index")])
        indexed <- indexed[order(indexed$stage_index), , drop = FALSE]
        stage_order <- indexed$stage
    }
    stage_order <- as.character(stage_order)
    if (length(stage_order) < 3L || anyNA(stage_order) ||
            anyDuplicated(stage_order)) {
        stop("A unique stage order with at least three stages is required.",
            call. = FALSE)
    }
    if (any(!as.character(stage_data$stage) %in% stage_order)) {
        stop("'stage_order' omits stages in 'stage_data'.", call. = FALSE)
    }
    if (is.null(stage_coordinates)) {
        coordinates <- rep(NA_real_, length(stage_order))
        names(coordinates) <- stage_order
    } else {
        if (!is.numeric(stage_coordinates) || anyNA(stage_coordinates) ||
                any(!is.finite(stage_coordinates))) {
            stop("'stage_coordinates' must contain finite numbers.",
                call. = FALSE)
        }
        if (!is.null(names(stage_coordinates)) &&
                all(nzchar(names(stage_coordinates)))) {
            if (!setequal(names(stage_coordinates), stage_order)) {
                stop(
                    "Named 'stage_coordinates' must match 'stage_order'.",
                    call. = FALSE
                )
            }
            coordinates <- as.numeric(stage_coordinates[stage_order])
        } else {
            if (length(stage_coordinates) != length(stage_order)) {
                stop(
                    "Unnamed 'stage_coordinates' must match 'stage_order' ",
                    "length.",
                    call. = FALSE
                )
            }
            coordinates <- as.numeric(stage_coordinates)
        }
        if (any(diff(coordinates) <= 0)) {
            stop("'stage_coordinates' must be strictly increasing.",
                call. = FALSE)
        }
        names(coordinates) <- stage_order
    }

    unique_key <- .stable_key(
        stage_data$feature_id, stage_data$focal_group, stage_data$stage
    )
    if (anyDuplicated(unique_key)) {
        stop(
            "'stage_data' must have one row per feature/focal group/stage.",
            call. = FALSE
        )
    }

    gene_ok <- if (is.null(gene_q_threshold)) {
        rep(TRUE, nrow(stage_data))
    } else {
        is.finite(stage_data$gene_q) &
            stage_data$gene_q < gene_q_threshold
    }
    eligible <- stage_data$complete_comparisons &
        stage_data$n_comparators >= min_comparators &
        stage_data$direction_coherent &
        stage_data$direction_sign != 0 &
        is.finite(stage_data$worst_q) &
        stage_data$worst_q < q_threshold &
        is.finite(stage_data$min_abs_effect) &
        .greater_equal_numeric(
            stage_data$min_abs_effect, effect_threshold
        ) &
        is.finite(stage_data$comparator_spread) &
        stage_data$comparator_spread < comparator_tolerance &
        gene_ok
    stage_data$state <- ifelse(
        eligible, as.integer(stage_data$direction_sign), 0L
    )

    split_index <- split(
        seq_len(nrow(stage_data)),
        .stable_key(stage_data$feature_id, stage_data$focal_group)
    )
    output <- vector("list", 0L)
    output_index <- 0L

    for (index in split_index) {
        observed <- stage_data[index, , drop = FALSE]
        observed <- observed[order(observed$stage_index), , drop = FALSE]
        template <- data.frame(
            stage = stage_order,
            stage_index = seq_along(stage_order),
            stringsAsFactors = FALSE
        )
        matched <- match(template$stage, as.character(observed$stage))
        template$observed <- !is.na(matched)
        template$state <- 0L
        template$max_abs_effect <- NA_real_
        template$mean_effect <- NA_real_
        template$worst_q <- NA_real_
        template$gene_q <- NA_real_
        present <- which(template$observed)
        template$state[present] <- observed$state[matched[present]]
        template$max_abs_effect[present] <-
            observed$max_abs_effect[matched[present]]
        template$mean_effect[present] <- observed$mean_effect[matched[present]]
        template$worst_q[present] <- observed$worst_q[matched[present]]
        template$gene_q[present] <- observed$gene_q[matched[present]]

        runs <- rle(template$state)
        starts <- cumsum(c(1L, utils::head(runs$lengths, -1L)))
        ends <- cumsum(runs$lengths)
        nonzero <- which(runs$values != 0L)
        if (!length(nonzero)) next

        for (run_index in nonzero) {
            start <- starts[[run_index]]
            end <- ends[[run_index]]
            episode_length <- end - start + 1L
            if (episode_length < min_episode_stages ||
                    episode_length > max_episode_stages) {
                next
            }
            if (start <= flank_width ||
                    end + flank_width > length(stage_order)) {
                next
            }
            flanks <- c(
                seq.int(start - flank_width, start - 1L),
                seq.int(end + 1L, end + flank_width)
            )
            if (!all(template$observed[flanks])) next
            flank_values <- template$max_abs_effect[flanks]
            finite_flanks <- is.finite(flank_values)
            if (flank_missing == "complete" && !all(finite_flanks)) next
            if (!any(finite_flanks)) next
            flanking_similarity <- max(flank_values[finite_flanks])
            if (flanking_similarity >= flank_tolerance) next

            episode_rows <- seq.int(start, end)
            output_index <- output_index + 1L
            gene_q <- template$gene_q[episode_rows]
            output[[output_index]] <- data.frame(
                feature_id = as.character(observed$feature_id[[1L]]),
                gene_id = as.character(observed$gene_id[[1L]]),
                gene_name = as.character(observed$gene_name[[1L]]),
                focal_group = as.character(observed$focal_group[[1L]]),
                start_stage = stage_order[[start]],
                end_stage = stage_order[[end]],
                start_index = start,
                end_index = end,
                n_stages = episode_length,
                start_coordinate = coordinates[[start]],
                end_coordinate = coordinates[[end]],
                coordinate_span = coordinates[[end]] - coordinates[[start]],
                direction = if (runs$values[[run_index]] > 0L) {
                    "higher"
                } else {
                    "lower"
                },
                max_abs_usage_difference = max(
                    abs(template$mean_effect[episode_rows])
                ),
                mean_abs_usage_difference = mean(
                    abs(template$mean_effect[episode_rows])
                ),
                weakest_component_q = max(
                    template$worst_q[episode_rows]
                ),
                flanking_max_abs_difference = flanking_similarity,
                n_required_flanks = length(flanks),
                n_finite_flanks = sum(finite_flanks),
                flanking_complete = all(finite_flanks),
                gene_q = if (all(is.na(gene_q))) {
                    NA_real_
                } else {
                    max(gene_q, na.rm = TRUE)
                },
                replicate_separation = NA_real_,
                replicate_consistent = NA,
                reciprocal_exchange = FALSE,
                stringsAsFactors = FALSE
            )
        }
    }

    if (!length(output)) {
        answer <- .empty_episode_table()
    } else {
        answer <- do.call(rbind, output)
        rownames(answer) <- NULL
        answer <- answer[order(
            answer$weakest_component_q,
            -answer$max_abs_usage_difference,
            -answer$n_stages,
            answer$gene_id,
            answer$feature_id,
            answer$focal_group,
            answer$start_index,
            answer$end_index,
            answer$direction
        ), , drop = FALSE]
    }
    result <- S4Vectors::DataFrame(answer)
    attr(result, "stage_order") <- stage_order
    attr(result, "expected_groups") <- expected_groups
    attr(result, "decision_parameters") <- list(
        q_threshold = q_threshold,
        effect_threshold = effect_threshold,
        comparator_tolerance = comparator_tolerance,
        flank_tolerance = flank_tolerance,
        gene_q_threshold = gene_q_threshold,
        min_comparators = as.integer(min_comparators),
        min_episode_stages = as.integer(min_episode_stages),
        max_episode_stages = max_episode_stages,
        flank_width = as.integer(flank_width),
        stage_coordinates = coordinates,
        flank_missing = flank_missing
    )
    result
}
