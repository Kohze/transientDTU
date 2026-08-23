#' Results from a transient transcript-usage analysis
#'
#' `TransientDTUResult` stores the standardized stage-level evidence, detected
#' episodes, ranked candidates, diagnostics, parameters, and the matched call.
#' Tables are stored as [S4Vectors::DataFrame] objects so they interoperate with
#' the Bioconductor ecosystem without coercing identifiers to factors.
#'
#' @slot stageData Standardized feature-by-focal-group-by-stage evidence.
#' @slot episodes Detected and annotated episode table.
#' @slot candidates Deterministically ranked candidate groups.
#' @slot diagnostics Named analysis counts and messages.
#' @slot parameters Complete decision-rule parameters.
#' @slot call The matched analysis call.
#'
#' @return A `TransientDTUResult` class definition.
#'
#' @name TransientDTUResult-class
#' @rdname TransientDTUResult-class
#' @importClassesFrom S4Vectors DFrame
#' @exportClass TransientDTUResult
methods::setClass(
    "TransientDTUResult",
    slots = c(
        stageData = "DFrame",
        episodes = "DFrame",
        candidates = "DFrame",
        diagnostics = "DFrame",
        parameters = "list",
        call = "language"
    )
)

methods::setValidity("TransientDTUResult", function(object) {
    stage_required <- c(
        "feature_id", "gene_id", "focal_group", "stage", "stage_index"
    )
    episode_required <- c(
        "feature_id", "gene_id", "focal_group", "start_stage",
        "end_stage", "direction"
    )
    candidate_required <- c("gene_id", "rank")

    missing_stage <- setdiff(stage_required, colnames(object@stageData))
    missing_episode <- setdiff(episode_required, colnames(object@episodes))
    missing_candidate <- setdiff(
        candidate_required, colnames(object@candidates)
    )
    messages <- character()
    if (length(missing_stage)) {
        messages <- c(
            messages,
            paste("stageData lacks:", paste(missing_stage, collapse = ", "))
        )
    }
    if (length(missing_episode)) {
        messages <- c(
            messages,
            paste("episodes lacks:", paste(missing_episode, collapse = ", "))
        )
    }
    if (length(missing_candidate)) {
        messages <- c(
            messages,
            paste(
                "candidates lacks:", paste(missing_candidate, collapse = ", ")
            )
        )
    }
    if (length(messages)) messages else TRUE
})

.new_transient_result <- function(
    stage_data,
    episodes,
    candidates,
    diagnostics,
    parameters,
    call
) {
    methods::new(
        "TransientDTUResult",
        stageData = S4Vectors::DataFrame(stage_data),
        episodes = S4Vectors::DataFrame(episodes),
        candidates = S4Vectors::DataFrame(candidates),
        diagnostics = S4Vectors::DataFrame(diagnostics),
        parameters = parameters,
        call = call
    )
}

#' Extract result tables and parameters
#'
#' @param x A [TransientDTUResult-class] object.
#'
#' @return `stageTable()`, `episodeTable()`, `candidateTable()`, and
#'   `diagnosticTable()` return [S4Vectors::DataFrame] objects.
#'   `decisionParameters()` returns a named list.
#'
#' @examples
#' example <- simulateTransientDTU(n_genes = 8, seed = 10)
#' fit <- runTransientDTU(example$pairwise, stage_order = example$stage_order)
#' head(episodeTable(fit))
#' decisionParameters(fit)$q_threshold
#'
#' @name result-accessors
NULL

#' @rdname result-accessors
#' @export
stageTable <- function(x) {
    if (!methods::is(x, "TransientDTUResult")) {
        stop("'x' must be a TransientDTUResult.", call. = FALSE)
    }
    x@stageData
}

#' @rdname result-accessors
#' @export
episodeTable <- function(x) {
    if (!methods::is(x, "TransientDTUResult")) {
        stop("'x' must be a TransientDTUResult.", call. = FALSE)
    }
    x@episodes
}

#' @rdname result-accessors
#' @export
candidateTable <- function(x) {
    if (!methods::is(x, "TransientDTUResult")) {
        stop("'x' must be a TransientDTUResult.", call. = FALSE)
    }
    x@candidates
}

#' @rdname result-accessors
#' @export
diagnosticTable <- function(x) {
    if (!methods::is(x, "TransientDTUResult")) {
        stop("'x' must be a TransientDTUResult.", call. = FALSE)
    }
    x@diagnostics
}

#' @rdname result-accessors
#' @export
decisionParameters <- function(x) {
    if (!methods::is(x, "TransientDTUResult")) {
        stop("'x' must be a TransientDTUResult.", call. = FALSE)
    }
    x@parameters
}

#' Display a transient DTU result
#'
#' @param object A [TransientDTUResult-class] object.
#'
#' @return `object`, invisibly.
#'
#' @rdname TransientDTUResult-class
#' @importMethodsFrom methods show
#' @export
methods::setMethod("show", "TransientDTUResult", function(object) {
    cat("class: TransientDTUResult\n")
    cat("stage rows:", nrow(object@stageData), "\n")
    cat("episodes:", nrow(object@episodes), "\n")
    cat("candidate groups:", nrow(object@candidates), "\n")
    if (nrow(object@episodes)) {
        cat(
            "genes with episodes:",
            length(unique(as.character(object@episodes$gene_id))),
            "\n"
        )
    }
    invisible(object)
})
