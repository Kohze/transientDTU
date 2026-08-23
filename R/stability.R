#' Evaluate a joint decision-threshold grid
#'
#' Re-runs the complete frozen detector across all joint combinations of the
#' principal thresholds and records episode counts and panel membership.
#'
#' @param stage_data Standardized stage evidence from [makeStageDTU()].
#' @param q_threshold,effect_threshold,comparator_tolerance,flank_tolerance
#'   Numeric vectors defining a Cartesian threshold grid.
#' @param gene_q_threshold List or numeric vector of gene-level thresholds.
#'   Include `NULL` in a list to evaluate no gene-level screen.
#' @param usage,col_data Optional replicate usage and metadata.
#' @param panel_size Candidate panel size.
#' @param reciprocal_only Restrict rankings to reciprocal groups.
#' @param min_comparators Minimum number of comparator groups.
#' @param replicate_action `"none"`, `"annotate"`, or `"filter"`.
#' @param ... Additional arguments passed to [checkReplicateSeparation()].
#'
#' @return A list containing `summary`, `selections`, and `selection_frequency`
#'   [S4Vectors::DataFrame] objects.
#'
#' @examples
#' simulated <- simulateTransientDTU(n_genes = 16, seed = 9)
#' stage <- makeStageDTU(simulated$pairwise, simulated$stage_order)
#' stability <- thresholdStability(
#'     stage,
#'     q_threshold = c(0.01, 0.05),
#'     effect_threshold = c(0.10, 0.15),
#'     usage = simulated$usage,
#'     col_data = simulated$col_data,
#'     replicate_action = "filter"
#' )
#' stability$selection_frequency
#'
#' @export
thresholdStability <- function(
    stage_data,
    q_threshold = c(0.01, 0.05),
    effect_threshold = c(0.10, 0.15, 0.20),
    comparator_tolerance = effect_threshold,
    flank_tolerance = effect_threshold,
    gene_q_threshold = list(NULL),
    usage = NULL,
    col_data = NULL,
    panel_size = 6L,
    reciprocal_only = TRUE,
    min_comparators = 2L,
    replicate_action = c("none", "annotate", "filter"),
    ...
) {
    stage_order <- attr(stage_data, "stage_order", exact = TRUE)
    replicate_action <- match.arg(replicate_action)
    if (replicate_action != "none" && is.null(usage)) {
        stop("'usage' is required for replicate checking.", call. = FALSE)
    }
    for (name in c(
        "q_threshold", "effect_threshold", "comparator_tolerance",
        "flank_tolerance"
    )) {
        value <- get(name)
        if (!is.numeric(value) || !length(value) || anyNA(value) ||
                any(!is.finite(value)) || any(value < 0)) {
            stop("'", name, "' must contain finite non-negative numbers.",
                call. = FALSE)
        }
    }
    if (!is.list(gene_q_threshold)) {
        gene_q_threshold <- as.list(gene_q_threshold)
    }
    gene_labels <- vapply(gene_q_threshold, function(value) {
        if (is.null(value)) "none" else format(value, scientific = FALSE)
    }, character(1))
    grid <- expand.grid(
        q_threshold = q_threshold,
        effect_threshold = effect_threshold,
        comparator_tolerance = comparator_tolerance,
        flank_tolerance = flank_tolerance,
        gene_q_index = seq_along(gene_q_threshold),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    grid$setting_id <- sprintf("setting_%04d", seq_len(nrow(grid)))
    grid$gene_q_threshold <- gene_labels[grid$gene_q_index]
    summaries <- vector("list", nrow(grid))
    selections <- vector("list", nrow(grid))

    for (index in seq_len(nrow(grid))) {
        setting <- grid[index, , drop = FALSE]
        episodes <- detectEpisodes(
            stage_data,
            stage_order = stage_order,
            q_threshold = setting$q_threshold,
            effect_threshold = setting$effect_threshold,
            comparator_tolerance = setting$comparator_tolerance,
            flank_tolerance = setting$flank_tolerance,
            gene_q_threshold = gene_q_threshold[[setting$gene_q_index]],
            min_comparators = min_comparators
        )
        if (replicate_action != "none") {
            episodes <- checkReplicateSeparation(
                episodes, usage, col_data,
                keep = if (replicate_action == "filter") {
                    "consistent"
                } else {
                    "all"
                },
                ...
            )
        }
        episodes <- annotateReciprocal(episodes)
        candidates <- rankCandidates(
            episodes,
            n = panel_size,
            reciprocal_only = reciprocal_only,
            unique_genes = TRUE
        )
        summaries[[index]] <- data.frame(
            setting_id = setting$setting_id,
            q_threshold = setting$q_threshold,
            effect_threshold = setting$effect_threshold,
            comparator_tolerance = setting$comparator_tolerance,
            flank_tolerance = setting$flank_tolerance,
            gene_q_threshold = setting$gene_q_threshold,
            episodes = nrow(episodes),
            genes = length(unique(as.character(episodes$gene_id))),
            reciprocal_episodes = sum(
                episodes$reciprocal_exchange %in% TRUE
            ),
            panel_genes = nrow(candidates),
            stringsAsFactors = FALSE
        )
        candidate_frame <- as.data.frame(candidates)
        if (nrow(candidate_frame)) {
            candidate_frame$setting_id <- setting$setting_id
            selections[[index]] <- candidate_frame
        }
    }
    summary <- do.call(rbind, summaries)
    nonempty <- selections[vapply(selections, Negate(is.null), logical(1))]
    selection_table <- if (length(nonempty)) {
        do.call(rbind, nonempty)
    } else {
        empty <- .empty_candidate_table()
        empty$setting_id <- character()
        empty
    }
    if (nrow(selection_table)) {
        frequency <- as.data.frame(table(selection_table$gene_id),
            stringsAsFactors = FALSE)
        colnames(frequency) <- c("gene_id", "selected_settings")
        frequency$selection_frequency <-
            frequency$selected_settings / nrow(grid)
        frequency <- frequency[order(
            -frequency$selection_frequency, frequency$gene_id
        ), , drop = FALSE]
    } else {
        frequency <- data.frame(
            gene_id = character(), selected_settings = integer(),
            selection_frequency = numeric(), stringsAsFactors = FALSE
        )
    }
    list(
        summary = S4Vectors::DataFrame(summary),
        selections = S4Vectors::DataFrame(selection_table),
        selection_frequency = S4Vectors::DataFrame(frequency)
    )
}
