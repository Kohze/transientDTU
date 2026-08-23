#' Run the complete transient-DTU decision workflow
#'
#' `runTransientDTU()` standardizes generic pairwise evidence, detects bounded
#' episodes, optionally enforces complete replicate separation, annotates
#' reciprocal exchanges, and creates a deterministic candidate ranking.
#'
#' @param pairwise Generic pairwise DTU result table; see [makeStageDTU()].
#' @param stage_order Complete biological stage order.
#' @param usage Optional feature-by-sample usage matrix or
#'   [SummarizedExperiment::SummarizedExperiment].
#' @param col_data Sample metadata required for matrix `usage`.
#' @param replicate_action `"none"`, `"annotate"`, or `"filter"`. If `NULL`,
#'   it defaults to `"filter"` when `usage` is supplied and `"none"`
#'   otherwise.
#' @param incomplete_replicates Action for incomplete replicate information;
#'   see [checkReplicateSeparation()].
#' @param panel_size Number of ranked candidates, or `Inf`.
#' @param reciprocal_only,unique_genes Candidate-ranking controls.
#' @param ... Arguments passed to [makeStageDTU()] and [detectEpisodes()].
#'   Names shared by both functions are routed to both when appropriate.
#'
#' @return A validated [TransientDTUResult-class] object.
#'
#' @examples
#' simulated <- simulateTransientDTU(n_genes = 20, seed = 6)
#' result <- runTransientDTU(
#'     simulated$pairwise,
#'     stage_order = simulated$stage_order,
#'     usage = simulated$usage,
#'     col_data = simulated$col_data,
#'     panel_size = 6
#' )
#' result
#' candidateTable(result)
#'
#' @export
runTransientDTU <- function(
    pairwise,
    stage_order,
    usage = NULL,
    col_data = NULL,
    replicate_action = NULL,
    incomplete_replicates = c("drop", "keep", "error"),
    panel_size = 6L,
    reciprocal_only = TRUE,
    unique_genes = TRUE,
    ...
) {
    call <- match.call()
    dots <- list(...)
    if (length(dots) &&
            (is.null(names(dots)) || any(!nzchar(names(dots))))) {
        stop("All arguments in '...' must be named.", call. = FALSE)
    }
    incomplete_replicates <- match.arg(incomplete_replicates)
    if (is.null(replicate_action)) {
        replicate_action <- if (is.null(usage)) "none" else "filter"
    }
    replicate_action <- match.arg(
        replicate_action, c("none", "annotate", "filter")
    )
    if (replicate_action != "none" && is.null(usage)) {
        stop("'usage' is required for replicate checking.", call. = FALSE)
    }

    make_names <- names(formals(makeStageDTU))
    detect_names <- names(formals(detectEpisodes))
    replicate_names <- names(formals(checkReplicateSeparation))
    reserved <- c("pairwise", "stage_order", "stage_data", "episodes", "usage")
    unknown <- setdiff(names(dots), union(
        union(make_names, detect_names), replicate_names
    ))
    if (length(unknown)) {
        stop(
            "Unknown workflow arguments: ", paste(unknown, collapse = ", "),
            ".", call. = FALSE
        )
    }
    make_args <- dots[intersect(names(dots), setdiff(make_names, reserved))]
    detect_args <- dots[intersect(names(dots), setdiff(detect_names, reserved))]
    replicate_args <- dots[
        intersect(names(dots), setdiff(replicate_names, reserved))
    ]

    stage_data <- do.call(
        makeStageDTU,
        c(list(pairwise = pairwise, stage_order = stage_order), make_args)
    )
    episodes_before_replicates <- do.call(
        detectEpisodes,
        c(list(stage_data = stage_data, stage_order = stage_order), detect_args)
    )
    episodes <- episodes_before_replicates
    incomplete_count <- 0L
    if (replicate_action != "none") {
        annotated <- do.call(
            checkReplicateSeparation,
            c(
                list(
                    episodes = episodes,
                    usage = usage,
                    col_data = col_data,
                    incomplete = incomplete_replicates,
                    keep = "all"
                ),
                replicate_args
            )
        )
        incomplete_count <- sum(is.na(annotated$replicate_separation))
        episodes <- if (replicate_action == "filter") {
            annotated[annotated$replicate_consistent %in% TRUE, , drop = FALSE]
        } else {
            annotated
        }
        attr(episodes, "stage_order") <- stage_order
        attr(episodes, "expected_groups") <-
            attr(annotated, "expected_groups", exact = TRUE)
        attr(episodes, "decision_parameters") <-
            attr(annotated, "decision_parameters", exact = TRUE)
    }
    episodes <- annotateReciprocal(episodes)
    candidates <- rankCandidates(
        episodes,
        n = panel_size,
        reciprocal_only = reciprocal_only,
        replicate_consistent_only = FALSE,
        unique_genes = unique_genes
    )

    parameters <- c(
        attr(episodes_before_replicates, "decision_parameters", exact = TRUE),
        list(
            stage_order = as.character(stage_order),
            replicate_action = replicate_action,
            incomplete_replicates = incomplete_replicates,
            panel_size = panel_size,
            reciprocal_only = reciprocal_only,
            unique_genes = unique_genes,
            make_stage_arguments = make_args,
            replicate_arguments = replicate_args
        )
    )
    diagnostics <- data.frame(
        metric = c(
            "pairwise_rows", "stage_rows", "episodes_before_replicates",
            "episodes_incomplete_replicates", "episodes_retained",
            "genes_retained", "reciprocal_episodes", "candidate_groups"
        ),
        value = c(
            nrow(pairwise), nrow(stage_data), nrow(episodes_before_replicates),
            incomplete_count, nrow(episodes),
            length(unique(as.character(episodes$gene_id))),
            sum(episodes$reciprocal_exchange %in% TRUE), nrow(candidates)
        ),
        stringsAsFactors = FALSE
    )
    .new_transient_result(
        stage_data, episodes, candidates, diagnostics, parameters, call
    )
}
