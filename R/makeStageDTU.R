#' Standardize pairwise DTU evidence by focal group and stage
#'
#' `makeStageDTU()` converts a generic long pairwise contrast table into the
#' stage-level evidence used by [detectEpisodes()]. Upstream inference may come
#' from any method, provided each row supplies a signed focal-minus-comparator
#' usage effect and an adjusted p-value.
#'
#' @param pairwise A `data.frame` or [S4Vectors::DataFrame].
#' @param stage_order Character vector giving the complete biological order of
#'   stages. Values must be unique and must include every observed stage.
#' @param feature_col,gene_col,focal_col,comparator_col Column names defining
#'   feature, gene, focal-group, and comparator-group identifiers.
#' @param stage_col,effect_col,q_col Column names defining stage, signed effect,
#'   and adjusted component p-value values.
#' @param gene_name_col Optional display-name column.
#' @param gene_q_col Optional gene-level adjusted p-value column. Repeated
#'   values within a feature, focal group, and stage are conservatively reduced
#'   to their maximum.
#' @param expected_groups Optional complete set of biological groups. By
#'   default it is inferred from the union of focal and comparator values.
#' @param min_comparators Minimum number of distinct comparison groups required
#'   for a stage to be eligible. The default of two implements a focal group
#'   separating from at least two alternatives.
#' @param require_all_comparators Whether every non-focal `expected_groups`
#'   member must be observed. Recommended and `TRUE` by default.
#'
#' @return A [S4Vectors::DataFrame] with one row per observed feature, focal
#'   group, and stage. It contains conservative component evidence, effect
#'   summaries, comparator coherence, completeness, and the ordered stage
#'   index.
#'
#' @details
#' The comparator spread is `max(effect) - min(effect)`. Because all effects
#' use the same focal-minus-comparator orientation, this is exactly the range
#' of comparator-group means. Duplicate atomic contrasts are rejected rather
#' than silently reduced.
#'
#' @examples
#' simulated <- simulateTransientDTU(n_genes = 8, seed = 2)
#' stage <- makeStageDTU(simulated$pairwise, simulated$stage_order)
#' stage
#'
#' @export
makeStageDTU <- function(
    pairwise,
    stage_order,
    feature_col = "feature_id",
    gene_col = "gene_id",
    focal_col = "focal_group",
    comparator_col = "comparator_group",
    stage_col = "stage",
    effect_col = "effect",
    q_col = "q_value",
    gene_name_col = NULL,
    gene_q_col = NULL,
    expected_groups = NULL,
    min_comparators = 2L,
    require_all_comparators = TRUE
) {
    pairwise <- .assert_data_frame(pairwise, "pairwise")
    .assert_count(min_comparators, "min_comparators")
    if (length(stage_order) < 3L || anyNA(stage_order) ||
            anyDuplicated(stage_order)) {
        stop(
            "'stage_order' must contain at least three unique, ",
            "non-missing stages.",
            call. = FALSE
        )
    }
    stage_order <- as.character(stage_order)

    required <- c(
        feature_col, gene_col, focal_col, comparator_col, stage_col,
        effect_col, q_col
    )
    optional <- c(gene_name_col, gene_q_col)
    all_columns <- c(required, optional[!is.null(optional)])
    .assert_columns(pairwise, all_columns, "pairwise")

    canonical <- data.frame(
        feature_id = .as_id(pairwise[[feature_col]], feature_col),
        gene_id = .as_id(pairwise[[gene_col]], gene_col),
        focal_group = .as_id(pairwise[[focal_col]], focal_col),
        comparator_group = .as_id(
            pairwise[[comparator_col]], comparator_col
        ),
        stage = .as_id(pairwise[[stage_col]], stage_col),
        effect = as.numeric(pairwise[[effect_col]]),
        q_value = as.numeric(pairwise[[q_col]]),
        stringsAsFactors = FALSE
    )
    canonical$gene_name <- if (is.null(gene_name_col)) {
        canonical$gene_id
    } else {
        .as_id(pairwise[[gene_name_col]], gene_name_col, allow_na = TRUE)
    }
    canonical$gene_q <- if (is.null(gene_q_col)) {
        NA_real_
    } else {
        as.numeric(pairwise[[gene_q_col]])
    }

    if (any(canonical$focal_group == canonical$comparator_group)) {
        stop("A focal group cannot be compared with itself.", call. = FALSE)
    }
    if (any(!canonical$stage %in% stage_order)) {
        bad <- unique(canonical$stage[!canonical$stage %in% stage_order])
        stop(
            "Observed stages are absent from 'stage_order': ",
            paste(bad, collapse = ", "), ".",
            call. = FALSE
        )
    }
    if (any(!is.na(canonical$q_value) &
            (canonical$q_value < 0 | canonical$q_value > 1))) {
        stop("Adjusted p-values must lie in [0, 1] or be NA.", call. = FALSE)
    }
    if (any(!is.na(canonical$gene_q) &
            (canonical$gene_q < 0 | canonical$gene_q > 1))) {
        stop("Gene-level adjusted p-values must lie in [0, 1] or be NA.",
            call. = FALSE)
    }
    genes_per_feature <- tapply(
        canonical$gene_id, canonical$feature_id,
        function(value) length(unique(value))
    )
    if (any(genes_per_feature != 1L)) {
        stop("Each feature must map to exactly one gene identifier.",
            call. = FALSE)
    }

    atomic_key <- .stable_key(
        canonical$feature_id, canonical$focal_group,
        canonical$comparator_group, canonical$stage
    )
    if (anyDuplicated(atomic_key)) {
        stop(
            "Duplicate feature/focal/comparator/stage contrasts ",
            "are not allowed.",
            call. = FALSE
        )
    }

    groups <- if (is.null(expected_groups)) {
        sort(unique(c(canonical$focal_group, canonical$comparator_group)))
    } else {
        unique(.as_id(expected_groups, "expected_groups"))
    }
    if (length(groups) < min_comparators + 1L) {
        stop(
            "The design needs at least min_comparators + 1 groups.",
            call. = FALSE
        )
    }
    observed_groups <- c(
        canonical$focal_group, canonical$comparator_group
    )
    if (any(!observed_groups %in% groups)) {
        stop("'expected_groups' omits an observed group.", call. = FALSE)
    }

    split_index <- split(
        seq_len(nrow(canonical)),
        .stable_key(
            canonical$feature_id, canonical$focal_group, canonical$stage
        )
    )
    first_rows <- vapply(
        split_index, function(index) index[[1L]], integer(1)
    )
    summaries <- vapply(split_index, function(index) {
        value <- canonical[index, , drop = FALSE]
        effects <- value$effect
        q_values <- value$q_value
        finite_effects <- all(is.finite(effects))
        c(
            mean_effect = if (finite_effects) mean(effects) else NA_real_,
            min_abs_effect = if (finite_effects) {
                min(abs(effects))
            } else {
                NA_real_
            },
            max_abs_effect = if (finite_effects) {
                max(abs(effects))
            } else {
                NA_real_
            },
            min_effect = if (finite_effects) min(effects) else NA_real_,
            max_effect = if (finite_effects) max(effects) else NA_real_,
            comparator_spread = if (finite_effects) {
                diff(range(effects))
            } else {
                NA_real_
            },
            worst_q = if (all(is.finite(q_values))) {
                max(q_values)
            } else {
                NA_real_
            },
            gene_q = if (all(is.na(value$gene_q))) {
                NA_real_
            } else {
                max(value$gene_q, na.rm = TRUE)
            }
        )
    }, numeric(8))
    coherent <- is.finite(summaries["min_effect", ]) &
        (summaries["min_effect", ] > 0 |
            summaries["max_effect", ] < 0)
    comparator_count <- lengths(split_index)
    complete <- if (isTRUE(require_all_comparators)) {
        comparator_count == length(groups) - 1L
    } else {
        comparator_count >= min_comparators
    }
    answer <- data.frame(
        feature_id = canonical$feature_id[first_rows],
        gene_id = canonical$gene_id[first_rows],
        gene_name = canonical$gene_name[first_rows],
        focal_group = canonical$focal_group[first_rows],
        stage = canonical$stage[first_rows],
        stage_index = match(canonical$stage[first_rows], stage_order),
        n_comparators = as.integer(comparator_count),
        complete_comparisons = complete,
        direction_coherent = coherent,
        direction_sign = ifelse(
            coherent, sign(summaries["mean_effect", ]), 0
        ),
        mean_effect = summaries["mean_effect", ],
        min_abs_effect = summaries["min_abs_effect", ],
        max_abs_effect = summaries["max_abs_effect", ],
        comparator_spread = summaries["comparator_spread", ],
        worst_q = summaries["worst_q", ],
        gene_q = summaries["gene_q", ],
        stringsAsFactors = FALSE
    )
    rownames(answer) <- NULL
    answer <- answer[order(
        answer$feature_id, answer$focal_group, answer$stage_index
    ), , drop = FALSE]
    result <- S4Vectors::DataFrame(answer)
    attr(result, "stage_order") <- stage_order
    attr(result, "expected_groups") <- groups
    attr(result, "min_comparators") <- as.integer(min_comparators)
    attr(result, "require_all_comparators") <- require_all_comparators
    result
}
