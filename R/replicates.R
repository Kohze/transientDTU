.extract_usage <- function(
    usage,
    col_data = NULL,
    assay_name = "usage",
    sample_col = NULL,
    group_col = "group",
    stage_col = "stage"
) {
    if (methods::is(usage, "SummarizedExperiment")) {
        if (!assay_name %in% SummarizedExperiment::assayNames(usage)) {
            stop(
                "Assay '", assay_name, "' was not found in 'usage'.",
                call. = FALSE
            )
        }
        matrix_value <- as.matrix(
            SummarizedExperiment::assay(usage, assay_name)
        )
        metadata <- as.data.frame(SummarizedExperiment::colData(usage))
        sample_ids <- colnames(usage)
    } else {
        matrix_value <- as.matrix(usage)
        if (is.null(col_data)) {
            stop("'col_data' is required when 'usage' is a matrix.",
                call. = FALSE)
        }
        metadata <- .assert_data_frame(col_data, "col_data")
        sample_ids <- if (is.null(sample_col)) {
            rownames(metadata)
        } else {
            .assert_columns(metadata, sample_col, "col_data")
            as.character(metadata[[sample_col]])
        }
    }
    .assert_columns(metadata, c(group_col, stage_col), "col_data")
    if (is.null(rownames(matrix_value)) ||
            anyDuplicated(rownames(matrix_value))) {
        stop("'usage' needs unique feature identifiers as row names.",
            call. = FALSE)
    }
    if (is.null(colnames(matrix_value)) ||
            anyDuplicated(colnames(matrix_value))) {
        stop("'usage' needs unique sample identifiers as column names.",
            call. = FALSE)
    }
    if (is.null(sample_ids) || anyNA(sample_ids) || anyDuplicated(sample_ids)) {
        stop("Sample identifiers in 'col_data' must be unique and non-missing.",
            call. = FALSE)
    }
    matched <- match(colnames(matrix_value), sample_ids)
    if (anyNA(matched)) {
        stop(
            "Every usage column must match one 'col_data' sample.",
            call. = FALSE
        )
    }
    metadata <- metadata[matched, , drop = FALSE]
    rownames(metadata) <- colnames(matrix_value)
    storage.mode(matrix_value) <- "numeric"
    invalid <- !is.na(matrix_value) &
        (!is.finite(matrix_value) | matrix_value < 0 | matrix_value > 1)
    if (any(invalid)) {
        stop("Usage values must be finite proportions in [0, 1] or NA.",
            call. = FALSE)
    }
    list(
        usage = matrix_value,
        col_data = metadata,
        group_col = group_col,
        stage_col = stage_col
    )
}

#' Check complete replicate separation within detected episodes
#'
#' @param episodes Episode table returned by [detectEpisodes()].
#' @param usage A numeric feature-by-sample matrix or a
#'   [SummarizedExperiment::SummarizedExperiment] containing usage proportions.
#' @param col_data Sample metadata for matrix input. Row names, or `sample_col`,
#'   must match matrix columns.
#' @param assay_name Assay to use for `SummarizedExperiment` input.
#' @param sample_col Optional sample-identifier column in `col_data`.
#' @param group_col,stage_col Metadata columns defining groups and stages.
#' @param min_separation Strict lower bound for consistency. Zero requires
#'   non-overlap in the episode direction.
#' @param incomplete One of `"keep"`, `"drop"`, or `"error"`, controlling
#'   episodes lacking complete finite replicate values.
#' @param keep One of `"all"` or `"consistent"`.
#'
#' @return The episode [S4Vectors::DataFrame] with `replicate_separation` and
#'   `replicate_consistent` updated.
#'
#' @examples
#' simulated <- simulateTransientDTU(n_genes = 12, seed = 5)
#' stage <- makeStageDTU(simulated$pairwise, simulated$stage_order)
#' episodes <- detectEpisodes(stage)
#' checked <- checkReplicateSeparation(
#'     episodes, simulated$usage, simulated$col_data
#' )
#' checked
#'
#' @export
checkReplicateSeparation <- function(
    episodes,
    usage,
    col_data = NULL,
    assay_name = "usage",
    sample_col = NULL,
    group_col = "group",
    stage_col = "stage",
    min_separation = 0,
    incomplete = c("keep", "drop", "error"),
    keep = c("all", "consistent")
) {
    stage_order <- attr(episodes, "stage_order", exact = TRUE)
    expected_groups <- attr(episodes, "expected_groups", exact = TRUE)
    parameters <- attr(episodes, "decision_parameters", exact = TRUE)
    episodes <- .assert_data_frame(episodes, "episodes")
    .assert_columns(
        episodes,
        c(
            "feature_id", "focal_group", "start_index", "end_index",
            "direction", "replicate_separation", "replicate_consistent"
        ),
        "episodes"
    )
    .assert_number(min_separation, "min_separation", -Inf, Inf)
    incomplete <- match.arg(incomplete)
    keep <- match.arg(keep)
    extracted <- .extract_usage(
        usage, col_data, assay_name, sample_col, group_col, stage_col
    )
    if (is.null(stage_order)) {
        stop("Episode stage order is missing.", call. = FALSE)
    }
    metadata <- extracted$col_data
    matrix_value <- extracted$usage
    groups <- as.character(metadata[[group_col]])
    stages <- as.character(metadata[[stage_col]])
    if (is.null(expected_groups)) expected_groups <- unique(groups)
    if (any(!expected_groups %in% groups)) {
        stop("Replicate metadata omits an expected biological group.",
            call. = FALSE)
    }

    separation <- vapply(seq_len(nrow(episodes)), function(index) {
        feature <- as.character(episodes$feature_id[[index]])
        feature_index <- match(feature, rownames(matrix_value))
        if (is.na(feature_index)) return(NA_real_)
        focal <- as.character(episodes$focal_group[[index]])
        episode_stages <- stage_order[
            episodes$start_index[[index]]:episodes$end_index[[index]]
        ]
        stage_separation <- vapply(episode_stages, function(stage) {
            focal_columns <- groups == focal & stages == stage
            other_columns <- groups %in% setdiff(expected_groups, focal) &
                stages == stage
            focal_values <- matrix_value[feature_index, focal_columns]
            other_values <- matrix_value[feature_index, other_columns]
            if (!length(focal_values) || !length(other_values) ||
                    any(!is.finite(focal_values)) ||
                    any(!is.finite(other_values))) {
                return(NA_real_)
            }
            if (episodes$direction[[index]] == "higher") {
                min(focal_values) - max(other_values)
            } else {
                min(other_values) - max(focal_values)
            }
        }, numeric(1))
        if (anyNA(stage_separation)) NA_real_ else min(stage_separation)
    }, numeric(1))

    episodes$replicate_separation <- separation
    episodes$replicate_consistent <- is.finite(separation) &
        separation > min_separation
    incomplete_rows <- is.na(separation)
    if (any(incomplete_rows) && incomplete == "error") {
        stop(
            sum(incomplete_rows),
            " episodes lack complete replicate information.",
            call. = FALSE
        )
    }
    if (incomplete == "drop") {
        episodes <- episodes[!incomplete_rows, , drop = FALSE]
    }
    if (keep == "consistent") {
        episodes <- episodes[
            episodes$replicate_consistent %in% TRUE, , drop = FALSE
        ]
    }
    result <- S4Vectors::DataFrame(episodes)
    attr(result, "stage_order") <- stage_order
    attr(result, "expected_groups") <- expected_groups
    attr(result, "decision_parameters") <- parameters
    result
}
