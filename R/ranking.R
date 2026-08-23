#' Annotate reciprocal transcript-usage exchanges
#'
#' A reciprocal exchange contains at least one higher and one lower feature
#' episode for the same gene, focal group, and stage interval.
#'
#' @param episodes Episode table from [detectEpisodes()] or
#'   [checkReplicateSeparation()].
#' @param min_directions Number of distinct directions required. The only
#'   currently meaningful value is two.
#'
#' @return The episode [S4Vectors::DataFrame] with `reciprocal_exchange`
#'   replaced by the deterministic group annotation.
#'
#' @examples
#' simulated <- simulateTransientDTU(n_genes = 12, seed = 8)
#' episodes <- detectEpisodes(
#'     makeStageDTU(simulated$pairwise, simulated$stage_order)
#' )
#' annotateReciprocal(episodes)
#'
#' @export
annotateReciprocal <- function(episodes, min_directions = 2L) {
    stage_order <- attr(episodes, "stage_order", exact = TRUE)
    expected_groups <- attr(episodes, "expected_groups", exact = TRUE)
    parameters <- attr(episodes, "decision_parameters", exact = TRUE)
    episodes <- .assert_data_frame(episodes, "episodes")
    .assert_count(min_directions, "min_directions", 2L)
    if (min_directions != 2L) {
        stop("'min_directions' must be two for higher/lower episodes.",
            call. = FALSE)
    }
    .assert_columns(
        episodes,
        c(
            "gene_id", "focal_group", "start_stage", "end_stage",
            "direction", "reciprocal_exchange"
        ),
        "episodes"
    )
    if (!nrow(episodes)) {
        result <- S4Vectors::DataFrame(episodes)
        attr(result, "stage_order") <- stage_order
        attr(result, "expected_groups") <- expected_groups
        attr(result, "decision_parameters") <- parameters
        return(result)
    }
    key <- .stable_key(
        episodes$gene_id, episodes$focal_group,
        episodes$start_stage, episodes$end_stage
    )
    direction_count <- tapply(
        as.character(episodes$direction), key,
        function(value) length(unique(value))
    )
    episodes$reciprocal_exchange <- unname(direction_count[key] >= 2L)
    result <- S4Vectors::DataFrame(episodes)
    attr(result, "stage_order") <- stage_order
    attr(result, "expected_groups") <- expected_groups
    attr(result, "decision_parameters") <- parameters
    result
}

#' Rank transient transcript-usage candidate groups
#'
#' Reduces feature-level episodes to gene, focal-group, and interval groups,
#' then applies the fully specified stable ordering in `inst/ALGORITHM.md`.
#'
#' @param episodes Annotated episode table.
#' @param n Maximum number of candidates. Use `Inf` for all candidates.
#' @param reciprocal_only Restrict to reciprocal exchanges.
#' @param replicate_consistent_only Restrict to episodes passing replicate
#'   separation. An informative error is raised if no replicate check has been
#'   performed.
#' @param unique_genes Keep only the best-ranked interval for each gene.
#'
#' @return A deterministically ordered [S4Vectors::DataFrame] with a `rank`
#'   column.
#'
#' @examples
#' simulated <- simulateTransientDTU(n_genes = 20, seed = 3)
#' episodes <- annotateReciprocal(detectEpisodes(
#'     makeStageDTU(simulated$pairwise, simulated$stage_order)
#' ))
#' rankCandidates(episodes, n = 6)
#'
#' @export
rankCandidates <- function(
    episodes,
    n = Inf,
    reciprocal_only = TRUE,
    replicate_consistent_only = FALSE,
    unique_genes = TRUE
) {
    episodes <- .assert_data_frame(episodes, "episodes")
    required <- c(
        "feature_id", "gene_id", "gene_name", "focal_group",
        "start_stage", "end_stage", "start_index", "end_index", "n_stages",
        "start_coordinate", "end_coordinate", "coordinate_span",
        "direction", "weakest_component_q", "max_abs_usage_difference",
        "replicate_separation", "replicate_consistent", "reciprocal_exchange"
    )
    .assert_columns(episodes, required, "episodes")
    if (!(is.numeric(n) && length(n) == 1L && !is.na(n) &&
            (is.infinite(n) || (n >= 1 && n == as.integer(n))))) {
        stop("'n' must be a positive integer or Inf.", call. = FALSE)
    }
    if (isTRUE(reciprocal_only)) {
        episodes <- episodes[
            episodes$reciprocal_exchange %in% TRUE, , drop = FALSE
        ]
    }
    if (isTRUE(replicate_consistent_only)) {
        if (nrow(episodes) && all(is.na(episodes$replicate_consistent))) {
            stop(
                "Replicate consistency was requested but has not been checked.",
                call. = FALSE
            )
        }
        episodes <- episodes[
            episodes$replicate_consistent %in% TRUE, , drop = FALSE
        ]
    }
    if (!nrow(episodes)) return(S4Vectors::DataFrame(.empty_candidate_table()))

    key <- .stable_key(
        episodes$gene_id, episodes$focal_group,
        episodes$start_stage, episodes$end_stage
    )
    split_index <- split(seq_len(nrow(episodes)), key)
    rows <- lapply(split_index, function(index) {
        value <- episodes[index, , drop = FALSE]
        feature_ids <- sort(unique(as.character(value$feature_id)))
        directions <- sort(unique(as.character(value$direction)))
        separation <- value$replicate_separation
        data.frame(
            rank = NA_integer_,
            gene_id = as.character(value$gene_id[[1L]]),
            gene_name = as.character(value$gene_name[[1L]]),
            focal_group = as.character(value$focal_group[[1L]]),
            start_stage = as.character(value$start_stage[[1L]]),
            end_stage = as.character(value$end_stage[[1L]]),
            start_index = as.integer(value$start_index[[1L]]),
            end_index = as.integer(value$end_index[[1L]]),
            n_stages = max(value$n_stages),
            start_coordinate = as.numeric(value$start_coordinate[[1L]]),
            end_coordinate = as.numeric(value$end_coordinate[[1L]]),
            coordinate_span = as.numeric(value$coordinate_span[[1L]]),
            n_features = length(feature_ids),
            feature_ids = paste(feature_ids, collapse = ";"),
            directions = paste(directions, collapse = ";"),
            reciprocal_exchange = any(value$reciprocal_exchange %in% TRUE),
            max_weakest_component_q = max(value$weakest_component_q),
            max_abs_usage_difference = max(
                value$max_abs_usage_difference
            ),
            min_replicate_separation = if (all(is.na(separation))) {
                NA_real_
            } else {
                min(separation, na.rm = TRUE)
            },
            stringsAsFactors = FALSE
        )
    })
    answer <- do.call(rbind, rows)
    rownames(answer) <- NULL
    separation_order <- .order_na_last_desc(answer$min_replicate_separation)
    answer <- answer[order(
        answer$max_weakest_component_q,
        -answer$max_abs_usage_difference,
        -separation_order,
        -answer$n_stages,
        answer$gene_id,
        answer$focal_group,
        answer$start_index,
        answer$end_index,
        answer$feature_ids
    ), , drop = FALSE]
    if (isTRUE(unique_genes)) {
        answer <- answer[!duplicated(answer$gene_id), , drop = FALSE]
    }
    if (is.finite(n)) {
        answer <- utils::head(answer, as.integer(n))
    }
    answer$rank <- seq_len(nrow(answer))
    rownames(answer) <- NULL
    S4Vectors::DataFrame(answer)
}
